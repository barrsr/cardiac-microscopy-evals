"""Ground truth and naive baselines for task A2 (1:1 capture of 1 Hz field stimulation).

Usage:
    python ground_truth.py [data_dir] [--figures]

Data: Kobeissi et al., "Strain gauge platforms: Time-lapse microscopy dataset of engineered
cardiac microbundles", Zenodo 11120809 / Dryad 10.5061/dryad.sqv9s4nbg (CC0),
Type1_condition_025.zip -> condition_025/batch_000/Ex00628.tif .. Ex00638.tif, used unmodified.
Condition 025 is documented by the depositors as "taken at Day 7 with 1 Hz stimulation".

Expert definition of 1:1 capture used here: every stimulus produces exactly one beat, so the
beat train runs at the stimulus period (1.000 s) with no beat added or dropped and no
cumulative phase drift against a 1.000 s grid anchored on the first beat.
"""
import re
import sys
from datetime import datetime
from pathlib import Path

import cv2
import numpy as np
import tifffile
from scipy import signal

HERE = Path(__file__).resolve().parent
NOMINAL_FPS = 30.0      # "Imaging frequency: 30 Hz" in the depositors' metadata sheet
PACING_S = 1.000        # "with 1 Hz stimulation" in Type1_experimental_conditions.pdf

# Capture criteria (expert): mean interval within 2% of the pacing period, every single
# interval within 5% of it, and cumulative phase drift from the 1.000 s grid under 0.10 s.
MEAN_TOL, IBI_TOL, DRIFT_TOL = 0.02, 0.05, 0.10


def load(path):
    pat = re.compile(r'acquisition-time-local" type="time" value="([^"]+)"')
    with tifffile.TiffFile(path) as tf:
        stamps = [datetime.strptime(pat.search(p.description).group(1), "%Y%m%d %H:%M:%S.%f")
                  for p in tf.pages]
        stack = tf.asarray().astype(np.float32)
    t = np.array([(s - stamps[0]).total_seconds() for s in stamps])
    return stack, t


def peaks(y, frac=0.25, dist=4):
    p, _ = signal.find_peaks(y, prominence=frac * (y.max() - y.min()), distance=dist)
    return p


def subframe(y, i):
    if 0 < i < len(y) - 1:
        a, b, c = y[i - 1], y[i], y[i + 1]
        den = a - 2 * b + c
        if den != 0:
            return i + 0.5 * (a - c) / den
    return float(i)


def dominant_hz(y, fs, fmin=0.2, nfft=8192):
    y = y - y.mean()
    spec = np.abs(np.fft.rfft(y * np.hanning(len(y)), nfft))
    f = np.fft.rfftfreq(nfft, 1 / fs)
    m = f >= fmin
    return f[m][np.argmax(spec[m])]


def analyze(path):
    stack, t = load(path)
    n = len(stack)
    fs = (n - 1) / t[-1]
    # Motion signals on a 4x-downsampled stack (speed only drives the reference choice).
    small = stack[:, ::4, ::4]
    u8 = [cv2.normalize(f, None, 0, 255, cv2.NORM_MINMAX).astype(np.uint8) for f in small]
    speed = np.array([np.linalg.norm(cv2.calcOpticalFlowFarneback(
        u8[i], u8[i + 1], None, 0.5, 3, 9, 3, 5, 1.2, 0), axis=2).mean() for i in range(n - 1)])
    framediff = np.abs(np.diff(stack, axis=0)).mean(axis=(1, 2))
    # Expert signal: displacement from the most quiescent (relaxed) frame.
    win = np.convolve(speed, np.ones(5) / 5, mode="valid")
    ref = int(np.argmin(win) + 2)
    disp = np.abs(stack - stack[ref]).mean(axis=(1, 2))
    bi = peaks(disp, frac=0.25, dist=12)
    bf = np.array([subframe(disp, i) for i in bi])
    bt = np.interp(bf, np.arange(n), t)
    ibi = np.diff(bt)
    k = np.arange(len(bt))
    phase = bt - (bt[0] + k * PACING_S)           # offset from a 1.000 s grid on beat 1
    slope = np.polyfit(k, bt, 1)[0] if len(bt) > 1 else np.nan   # fitted period
    # Beat amplitudes (for alternans): peak displacement above the pre-beat minimum.
    amps = disp[bi] - np.min(disp)
    # Autocorrelation of the displacement trace: first peak lag.
    x = disp - disp.mean()
    ac = np.correlate(x, x, "full")[n - 1:]
    ac /= ac[0]
    acp = signal.find_peaks(ac, height=0.1)[0]
    ac_lag = int(acp[0]) if len(acp) else np.nan
    rec = dict(
        id=path.stem, n=n, fs=fs, span=t[-1], ref=ref, beats=bi.tolist(), ibi=ibi,
        mean_ibi=ibi.mean(), cv=100 * ibi.std() / ibi.mean(), bpm=60 / ibi.mean(),
        period_fit=slope, drift=phase[-1], max_dev=np.max(np.abs(ibi - PACING_S)),
        amps=amps, ac_period=ac_lag / fs, n_beats=len(bi),
        # Naive quantities
        bpm_30=60 / (np.diff(bf).mean() / NOMINAL_FPS),
        fd_peaks=len(peaks(framediff)), fd_ibi=np.diff(peaks(framediff)).mean() / fs,
        fd_fft=dominant_hz(framediff, fs), disp_fft=dominant_hz(disp, fs),
        disp_fft30=dominant_hz(disp, NOMINAL_FPS),
        count_bpm=60 * len(bi) / t[-1], count_bpm30=60 * len(bi) / (n / NOMINAL_FPS),
        t=t, disp=disp, bt=bt, phase=phase,
    )
    rec["captured"] = bool(abs(rec["mean_ibi"] - PACING_S) <= MEAN_TOL * PACING_S
                           and rec["max_dev"] <= IBI_TOL * PACING_S
                           and abs(rec["drift"]) <= DRIFT_TOL)
    return rec


def main():
    args = [a for a in sys.argv[1:] if not a.startswith("--")]
    ddir = Path(args[0]) if args else HERE / "data"
    files = sorted(ddir.glob("Ex*.tif"))
    recs = [analyze(f) for f in files]

    print("PER-TISSUE TABLE (timestamps; displacement from relaxed reference frame)")
    hdr = (f"{'file':<8} {'frames':>6} {'fps':>6} {'beats':>5} {'IBIs (s)':<32} {'bpm':>6} "
           f"{'CV%':>5} {'fit T':>6} {'drift s':>7} {'AC T':>6} {'amp ratio':>9} {'bpm@30':>6}  call")
    print(hdr)
    for r in recs:
        ibis = ",".join(f"{v:.3f}" for v in r["ibi"])
        ar = r["amps"].min() / r["amps"].max()
        print(f"{r['id']:<8} {r['n']:>6} {r['fs']:>6.2f} {r['n_beats']:>5} {ibis:<32} {r['bpm']:>6.1f} "
              f"{r['cv']:>5.1f} {r['period_fit']:>6.3f} {r['drift']:>+7.3f} {r['ac_period']:>6.3f} "
              f"{ar:>9.2f} {r['bpm_30']:>6.1f}  {'CAPTURED' if r['captured'] else 'not captured'}")
    cap = [r["id"] for r in recs if r["captured"]]
    notcap = sorted(r["id"] for r in recs if not r["captured"])
    print(f"\nGROUND TRUTH: n_captured = {len(cap)}; not_captured = {notcap}")

    # Naive approaches -------------------------------------------------------------
    def score(name, calls):
        nc = sorted(r["id"] for r, c in zip(recs, calls) if not c)
        ok = len(recs) - len(nc) == len(cap) and nc == notcap
        print(f"  {name:<66s} n={len(recs)-len(nc):>2}  not={','.join(x[-3:] for x in nc) or '-':<24s} "
              f"{'PASS' if ok else 'FAIL'}")

    print("\nNAIVE / ALTERNATIVE APPROACHES (PASS only if both n_captured and the set match)")
    score("Expert: timestamps, period + interval + phase criteria", [r["captured"] for r in recs])
    score("Rate within +/-2% of 60 bpm, timestamps", [abs(r["bpm"] - 60) <= 1.2 for r in recs])
    score("Rate within +/-3% of 60 bpm, timestamps", [abs(r["bpm"] - 60) <= 1.8 for r in recs])
    score("Rate within +/-5% of 60 bpm, timestamps", [abs(r["bpm"] - 60) <= 3 for r in recs])
    score("Rate within +/-10% of 60 bpm, timestamps", [abs(r["bpm"] - 60) <= 6 for r in recs])
    score("Rate rounds to 60 bpm (+/-2%), nominal 30 fps", [abs(r["bpm_30"] - 60) <= 1.2 for r in recs])
    score("Rate within +/-10% of 60 bpm, nominal 30 fps", [abs(r["bpm_30"] - 60) <= 6 for r in recs])
    score("Expert criteria, but nominal 30 fps time base",
          [abs(60 / r["bpm_30"] - 1) <= MEAN_TOL for r in recs])
    score("Captured = beats at or above 1 Hz (rate >= 58.8 bpm), timestamps", [r["bpm"] >= 58.8 for r in recs])
    score("Displacement FFT peak within 1.00 +/- 0.05 Hz, timestamps", [abs(r["disp_fft"] - 1) <= 0.05 for r in recs])
    score("Displacement FFT peak within 1.00 +/- 0.05 Hz, nominal 30 fps", [abs(r["disp_fft30"] - 1) <= 0.05 for r in recs])
    score("|frame diff| FFT peak within 1.00 +/- 0.05 Hz, timestamps", [abs(r["fd_fft"] - 1) <= 0.05 for r in recs])
    score("All |frame diff| peaks: mean interval within 5% of 1 s", [abs(r["fd_ibi"] - 1) <= 0.05 for r in recs])
    score("Beat count / duration within 5% of 60 bpm, timestamps", [abs(r["count_bpm"] - 60) <= 3 for r in recs])
    score("Regular rhythm only (CV < 3%), any rate", [r["cv"] < 3 for r in recs])

    if "--figures" in sys.argv:
        figure(recs)


def figure(recs):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    blue, orange, ink, muted, grid = "#2a78d6", "#eb6834", "#0b0b0b", "#52514e", "#d9d8d4"
    plt.rcParams.update({"font.size": 10, "axes.edgecolor": muted, "axes.labelcolor": ink,
                         "xtick.color": muted, "ytick.color": muted, "axes.spines.top": False,
                         "axes.spines.right": False})
    fig, ax = plt.subplots(figsize=(8, 4.8))
    ax.axhspan(-0.10, 0.10, color=grid, alpha=0.6, lw=0)
    ax.axhline(0, color=muted, lw=0.8)
    for r in recs:
        c = blue if r["captured"] else orange
        ax.plot(np.arange(len(r["phase"])), r["phase"], "-o", color=c, ms=4, lw=1.3, alpha=0.9)
        if not r["captured"]:
            ax.text(len(r["phase"]) - 1 + 0.08, r["phase"][-1], r["id"], fontsize=8, color=orange,
                    va="center")
    ax.plot([], [], "-o", color=blue, ms=4, label="1:1 capture")
    ax.plot([], [], "-o", color=orange, ms=4, label="Not captured")
    ax.set_xticks(range(5))
    ax.set_xlabel("Beat number (beat 0 anchors the grid)")
    ax.set_ylabel("Beat time minus 1.000 s pacing grid (s)")
    ax.set_title("Condition 025, 1 Hz field stimulation: phase of each beat against the pacing grid",
                 fontsize=10.5, loc="left", color=ink)
    ax.text(0.02, 0.105, "+/-0.10 s band", fontsize=8, color=muted, transform=ax.get_yaxis_transform())
    ax.legend(frameon=False, fontsize=9, loc="lower left")
    fig.tight_layout()
    out = HERE / "figures"
    out.mkdir(exist_ok=True)
    fig.savefig(out / "fig1_phase_vs_pacing_grid.png", dpi=200)
    plt.close(fig)
    print(f"\nFigure written to {out / 'fig1_phase_vs_pacing_grid.png'}")


if __name__ == "__main__":
    main()
