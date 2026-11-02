"""Ground truth and naive baselines for task A3 (1:1 capture versus escape rhythm under 1 Hz pacing).

Usage:
    python ground_truth.py [data_dir] [--figures]

Data: Kobeissi et al., "Strain gauge platforms: Time-lapse microscopy dataset of engineered
cardiac microbundles", Zenodo 11120809 / Dryad 10.5061/dryad.sqv9s4nbg (CC0),
Type1_condition_023.zip, nineteen of its twenty-three TIFFs, used unmodified. Condition 023 is
documented by the depositors as L550 MYBPC3 patient line, "taken at Day 7 with 1 Hz stimulation".

Expert definition of 1:1 capture used here: every stimulus produces exactly one beat, so every
beat-to-beat interval equals the stimulus period (1.000 s, within the timing precision of about
one frame) and the beats keep a fixed phase to a 1.000 s grid anchored on the first beat. A beat
that arrives earlier than one pacing period after the previous one cannot have been evoked by the
stimulus, so a tissue whose intervals are all a few percent short is beating on its own rhythm
(escaping the pacing), however regular it is.

Beat times are taken at the contraction upstroke (half-maximum crossing of the displacement
signal, interpolated between frames), converted with the per-plane MetaMorph timestamps. The
displacement peak is used as a second, independent beat marker.
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

# Expert capture criterion: every interval within 2% of the pacing period (about two thirds of a
# frame) and the last beat within 0.05 s of the 1.000 s grid anchored on the first beat.
IBI_TOL, DRIFT_TOL = 0.02, 0.05


def load(path):
    pat = re.compile(r'acquisition-time-local" type="time" value="([^"]+)"')
    with tifffile.TiffFile(path) as tf:
        stamps = [datetime.strptime(pat.search(p.description).group(1), "%Y%m%d %H:%M:%S.%f")
                  for p in tf.pages]
        stack = tf.asarray().astype(np.float32)
    t = np.array([(s - stamps[0]).total_seconds() for s in stamps])
    return stack, t


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
    small = stack[:, ::4, ::4]
    u8 = [cv2.normalize(f, None, 0, 255, cv2.NORM_MINMAX).astype(np.uint8) for f in small]
    speed = np.array([np.linalg.norm(cv2.calcOpticalFlowFarneback(
        u8[i], u8[i + 1], None, 0.5, 3, 9, 3, 5, 1.2, 0), axis=2).mean() for i in range(n - 1)])
    # Displacement from the most quiescent (relaxed) frame: one event per beat.
    win = np.convolve(speed, np.ones(5) / 5, mode="valid")
    ref = int(np.argmin(win) + 2)
    disp = np.abs(stack - stack[ref]).mean(axis=(1, 2))
    pk, _ = signal.find_peaks(disp, prominence=0.25 * (disp.max() - disp.min()), distance=12)
    # Beat marker 1 (primary): upstroke, half-maximum crossing before each peak.
    lo = np.percentile(disp, 5)
    up = []
    for i in pk:
        h = lo + 0.5 * (disp[i] - lo)
        j = i
        while j > 0 and disp[j - 1] > h:
            j -= 1
        if j == 0:
            continue                      # beat already under way in frame 0
        up.append((j - 1) + (h - disp[j - 1]) / (disp[j] - disp[j - 1]))
    up = np.array(up)
    bt = np.interp(up, np.arange(n), t)
    ibi = np.diff(bt)
    k = np.arange(len(bt))
    phase = bt - (bt[0] + k * PACING_S)
    # Beat marker 2: displacement peak with parabolic sub-frame interpolation.
    pf = np.array([subframe(disp, i) for i in pk])
    pt = np.interp(pf, np.arange(n), t)
    pibi = np.diff(pt)
    pphase = pt - (pt[0] + np.arange(len(pt)) * PACING_S)
    rec = dict(
        id=path.stem, n=n, fs=fs, n_beats=len(bt), ibi=ibi, mean_ibi=ibi.mean(),
        bpm=60 / ibi.mean(), cv=100 * ibi.std() / ibi.mean(), drift=phase[-1],
        max_dev=np.max(np.abs(ibi - PACING_S)), min_ibi=ibi.min(), phase=phase,
        pibi=pibi, pdrift=pphase[-1], pmax_dev=np.max(np.abs(pibi - PACING_S)),
        pbpm=60 / pibi.mean(),
        bpm_30=60 / (np.diff(up).mean() / NOMINAL_FPS),
        disp_fft=dominant_hz(disp, fs), count_bpm=60 * len(pk) / t[-1],
        t=t, disp=disp, bt=bt,
    )
    rec["captured"] = bool(rec["max_dev"] <= IBI_TOL * PACING_S and abs(rec["drift"]) <= DRIFT_TOL)
    rec["captured_peaks"] = bool(rec["pmax_dev"] <= IBI_TOL * PACING_S and abs(rec["pdrift"]) <= DRIFT_TOL)
    return rec


def main():
    args = [a for a in sys.argv[1:] if not a.startswith("--")]
    ddir = Path(args[0]) if args else HERE / "data"
    recs = [analyze(f) for f in sorted(ddir.glob("Ex*.tif"))]

    print("PER-TISSUE TABLE (timestamps; upstroke beat times; peak-time intervals for comparison)")
    print(f"{'file':<8} {'fps':>6} {'beats':>5} {'upstroke IBIs (s)':<28} {'bpm':>6} {'CV%':>5} "
          f"{'maxdev':>6} {'drift s':>7} {'peak IBIs (s)':<28} {'peak drift':>10} {'bpm@30':>6}  call")
    for r in recs:
        ibis = ",".join(f"{v:.3f}" for v in r["ibi"])
        pibis = ",".join(f"{v:.3f}" for v in r["pibi"])
        print(f"{r['id']:<8} {r['fs']:>6.2f} {r['n_beats']:>5} {ibis:<28} {r['bpm']:>6.1f} {r['cv']:>5.1f} "
              f"{r['max_dev']:>6.3f} {r['drift']:>+7.3f} {pibis:<28} {r['pdrift']:>+10.3f} {r['bpm_30']:>6.1f}  "
              f"{'CAPTURED' if r['captured'] else 'not captured'}"
              f"{'' if r['captured'] == r['captured_peaks'] else '  (peak marker disagrees)'}")
    cap = [r["id"] for r in recs if r["captured"]]
    notcap = sorted(r["id"] for r in recs if not r["captured"])
    print(f"\nGROUND TRUTH: n_captured = {len(cap)}; not_captured = {notcap}")

    def score(name, calls):
        nc = sorted(r["id"] for r, c in zip(recs, calls) if not c)
        ok = nc == notcap
        print(f"  {name:<72s} n={len(recs) - len(nc):>2}  not={','.join(x[-3:] for x in nc) or '-':<52s} "
              f"{'PASS' if ok else 'FAIL'}")

    print("\nNAIVE / ALTERNATIVE APPROACHES (PASS only if the not-captured set matches)")
    score("Expert: every interval within 2%, drift under 0.05 s (upstroke)", [r["captured"] for r in recs])
    score("Expert criterion on peak times", [r["captured_peaks"] for r in recs])
    score("A2 criteria: mean within 2%, every interval within 5%, drift < 0.10 s",
          [abs(r["mean_ibi"] - 1) <= 0.02 and r["max_dev"] <= 0.05 and abs(r["drift"]) <= 0.10 for r in recs])
    score("Every interval within 5%, drift under 0.10 s (no mean test)",
          [r["max_dev"] <= 0.05 and abs(r["drift"]) <= 0.10 for r in recs])
    score("Every interval within 5% of 1.000 s", [r["max_dev"] <= 0.05 for r in recs])
    score("Mean rate within 2% of 60 bpm", [abs(r["bpm"] - 60) <= 1.2 for r in recs])
    score("Mean rate within 3% of 60 bpm", [abs(r["bpm"] - 60) <= 1.8 for r in recs])
    score("Mean rate within 5% of 60 bpm", [abs(r["bpm"] - 60) <= 3 for r in recs])
    score("Mean rate within 2% and CV under 3%", [abs(r["bpm"] - 60) <= 1.2 and r["cv"] < 3 for r in recs])
    score("Mean rate within 5% and CV under 3%", [abs(r["bpm"] - 60) <= 3 and r["cv"] < 3 for r in recs])
    score("Phase drift under 0.10 s only", [abs(r["drift"]) <= 0.10 for r in recs])
    score("Rate at or above 58.8 bpm (beats at least as fast as the stimulus)", [r["bpm"] >= 58.8 for r in recs])
    score("Displacement FFT peak within 1.00 +/- 0.05 Hz", [abs(r["disp_fft"] - 1) <= 0.05 for r in recs])
    score("Nominal 30 fps; mean rate within 2% of 60 bpm", [abs(r["bpm_30"] - 60) <= 1.2 for r in recs])
    score("Nominal 30 fps; mean rate within 10% of 60 bpm", [abs(r["bpm_30"] - 60) <= 6 for r in recs])
    score("Peak times; every interval within 2%", [r["pmax_dev"] <= 0.02 for r in recs])
    score("Peak times; every interval within 5%, drift under 0.10 s",
          [r["pmax_dev"] <= 0.05 and abs(r["pdrift"]) <= 0.10 for r in recs])

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
    order = sorted(recs, key=lambda r: (not r["captured"], r["max_dev"]))
    fig, ax = plt.subplots(figsize=(8, 6.2))
    ax.axvspan(PACING_S - 0.02, PACING_S + 0.02, color=grid, alpha=0.8, lw=0)
    ax.axvline(PACING_S, color=muted, lw=0.8)
    for y, r in enumerate(order):
        c = blue if r["captured"] else orange
        ax.plot(r["ibi"], np.full(len(r["ibi"]), y), "o", color=c, ms=5, alpha=0.9)
        ax.plot([r["ibi"].min(), r["ibi"].max()], [y, y], "-", color=c, lw=0.8, alpha=0.5)
    ax.set_yticks(range(len(order)))
    ax.set_yticklabels([r["id"] for r in order], fontsize=8)
    ax.invert_yaxis()
    ax.set_xlim(0.4, 1.65)
    ax.set_xlabel("Beat-to-beat interval (s), upstroke times on the per-frame timestamps")
    ax.set_title("Condition 023, 1 Hz field stimulation: every beat interval of every tissue",
                 fontsize=10.5, loc="left", color=ink)
    ax.text(PACING_S + 0.025, len(order) - 0.6, "1.000 s +/- 0.02 s", fontsize=8, color=muted)
    ax.plot([], [], "o", color=blue, ms=5, label="1:1 capture")
    ax.plot([], [], "o", color=orange, ms=5, label="Not captured")
    ax.legend(frameon=False, fontsize=9, loc="upper right")
    fig.tight_layout()
    out = HERE / "figures"
    out.mkdir(exist_ok=True)
    fig.savefig(out / "fig1_intervals_per_tissue.png", dpi=200)
    plt.close(fig)
    print(f"\nFigure written to {out / 'fig1_intervals_per_tissue.png'}")


if __name__ == "__main__":
    main()
