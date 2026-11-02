"""Ground truth and naive baselines for task A (beat rate of a paced hiPSC-CM microbundle).

Usage:
    python ground_truth.py [path/to/Ex00628.tif] [--figures]

Data: Kobeissi et al., "Strain gauge platforms: Time-lapse microscopy dataset of engineered
cardiac microbundles", Zenodo 11120809 / Dryad 10.5061/dryad.sqv9s4nbg (CC0), file
Type1_condition_025.zip -> condition_025/batch_000/Ex00628.tif, used unmodified.
Condition 025 is documented by the depositors as "taken at Day 7 with 1 Hz stimulation".
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
PACING_HZ = 1.0         # "with 1 Hz stimulation" in Type1_experimental_conditions.pdf


def load(path):
    """Return the stack (float32) and per-plane acquisition times in seconds from frame 0."""
    pat = re.compile(r'acquisition-time-local" type="time" value="([^"]+)"')
    with tifffile.TiffFile(path) as tf:
        stamps = [datetime.strptime(pat.search(p.description).group(1), "%Y%m%d %H:%M:%S.%f")
                  for p in tf.pages]
        exposure = re.search(r"Exposure: ([0-9.]+) msec", tf.pages[0].description).group(1)
        stack = tf.asarray().astype(np.float32)
    t = np.array([(s - stamps[0]).total_seconds() for s in stamps])
    return stack, t, float(exposure)


def peaks(y, frac=0.25, dist=4):
    p, _ = signal.find_peaks(y, prominence=frac * (y.max() - y.min()), distance=dist)
    return p


def subframe(y, i):
    """Parabolic interpolation of a peak position (in samples)."""
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
    return f[m][np.argmax(spec[m])], f, spec


def main():
    args = [a for a in sys.argv[1:] if not a.startswith("--")]
    path = Path(args[0]) if args else HERE / "data" / "Ex00628.tif"
    stack, t, exposure_ms = load(path)
    n = len(stack)
    dt = np.diff(t)
    fs_true = (n - 1) / t[-1]

    print(f"File: {path.name}  frames={n}  shape={stack.shape[1:]}  dtype=uint16")
    print(f"Exposure in MetaMorph header: {exposure_ms:g} ms")
    print(f"Per-plane timestamps: span {t[-1]:.3f} s, dt median {np.median(dt)*1e3:.1f} ms "
          f"(min {dt.min()*1e3:.1f}, max {dt.max()*1e3:.1f}), mean rate {fs_true:.3f} frames/s")
    print(f"Nominal rate in metadata sheet: {NOMINAL_FPS:g} Hz -> implied span {(n-1)/NOMINAL_FPS:.3f} s")

    # Motion signals ----------------------------------------------------------------
    u8 = [cv2.normalize(f, None, 0, 255, cv2.NORM_MINMAX).astype(np.uint8) for f in stack]
    speed = np.array([np.linalg.norm(cv2.calcOpticalFlowFarneback(
        u8[i], u8[i + 1], None, 0.5, 3, 15, 3, 5, 1.2, 0), axis=2).mean() for i in range(n - 1)])
    framediff = np.abs(np.diff(stack, axis=0)).mean(axis=(1, 2))   # consecutive |diff|
    t_mid = 0.5 * (t[:-1] + t[1:])
    mean_int = stack.mean(axis=(1, 2))

    # Expert method: displacement from a quiescent (relaxed) reference frame.
    # Reference = center of the 5-frame window with the least optical-flow speed.
    win = np.convolve(speed, np.ones(5) / 5, mode="valid")
    ref_idx = int(np.argmin(win) + 2)
    disp = np.abs(stack - stack[ref_idx]).mean(axis=(1, 2))
    beat_idx = peaks(disp, frac=0.25, dist=15)
    beat_t = np.interp([subframe(disp, i) for i in beat_idx], np.arange(n), t)
    ibi = np.diff(beat_t)
    bpm_gt = 60.0 / ibi.mean()

    # Pair each beat with its contraction (rising) and relaxation (falling) speed peaks.
    vpk = peaks(speed)
    pairs = []
    for b in beat_idx:
        before = vpk[vpk < b]
        after = vpk[vpk >= b]
        if len(before) and len(after):
            pairs.append((int(before[-1]), int(b), int(after[0])))

    # Autocorrelation of the speed trace with harmonic check.
    x = speed - speed.mean()
    ac = np.correlate(x, x, "full")[len(x) - 1:]
    ac /= ac[0]
    ac_lag = signal.find_peaks(ac, height=0.2)[0][0]
    ac_period = ac_lag * np.mean(dt)
    bpm_ac = 60.0 / ac_period

    print("\nEXPERT")
    print(f"  reference (relaxed) frame: {ref_idx}")
    print(f"  displacement peaks (frames): {beat_idx.tolist()}")
    print(f"  inter-beat intervals (s, from timestamps): {np.round(ibi, 4).tolist()}")
    print(f"  IBI mean {ibi.mean():.4f} s, CV {100*ibi.std()/ibi.mean():.2f}%")
    print(f"  contraction / peak / relaxation frame triplets: {pairs}")
    print(f"  autocorrelation first lag {ac_lag} frames -> {bpm_ac:.1f} bpm")
    print(f"  GROUND TRUTH (measured): {bpm_gt:.2f} bpm  | pacing: {60*PACING_HZ:.0f} bpm")

    # Naive approaches ---------------------------------------------------------------
    dur = t[-1]
    res = {}
    fd_pk = peaks(framediff)
    res["All peaks of consecutive |frame diff| / duration"] = 60 * len(fd_pk) / dur
    sp_pk = peaks(speed)
    res["All peaks of optical-flow speed / duration"] = 60 * len(sp_pk) / dur
    res["Mean IBI of all |frame diff| peaks"] = 60.0 / (np.diff(fd_pk).mean() / fs_true)
    res["FFT dominant peak of |frame diff| (true fs)"] = 60 * dominant_hz(framediff, fs_true)[0]
    res["FFT dominant peak of optical-flow speed (true fs)"] = 60 * dominant_hz(speed, fs_true)[0]
    res["Expert displacement method, nominal 30 fps"] = 60.0 / (np.diff(
        [subframe(disp, i) for i in beat_idx]).mean() / NOMINAL_FPS)
    res["Displacement peaks counted / duration (true time)"] = 60 * len(beat_idx) / dur
    res["Complete beats (contraction + relaxation) / duration"] = 60 * len(pairs) / dur
    res["Displacement peaks / duration (150 frames at 30 fps = 5 s)"] = 60 * len(beat_idx) / (n / NOMINAL_FPS)
    res["FFT dominant peak of |frame diff| (30 fps)"] = 60 * dominant_hz(framediff, NOMINAL_FPS)[0]
    res["Mean-intensity trace FFT (true fs)"] = 60 * dominant_hz(mean_int, fs_true)[0]
    res["Mean-intensity trace FFT (30 fps)"] = 60 * dominant_hz(mean_int, NOMINAL_FPS)[0]
    res["Autocorrelation of speed, nominal 30 fps"] = 60.0 / (ac_lag / NOMINAL_FPS)

    print("\nNAIVE / ALTERNATIVE APPROACHES (bpm; PASS if within 60 +/- 3)")
    for k, v in res.items():
        print(f"  {k:<58s} {v:7.1f}  {'PASS' if abs(v-60) <= 3 else 'FAIL'}")
    print(f"  {'Expert displacement method, timestamps':<58s} {bpm_gt:7.1f}  "
          f"{'PASS' if abs(bpm_gt-60) <= 3 else 'FAIL'}")

    if "--figures" in sys.argv:
        make_figures(t, t_mid, speed, framediff, disp, beat_idx, beat_t, pairs, fs_true, ref_idx)


def make_figures(t, t_mid, speed, framediff, disp, beat_idx, beat_t, pairs, fs_true, ref_idx):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    blue, orange, ink, muted = "#2a78d6", "#eb6834", "#0b0b0b", "#52514e"
    plt.rcParams.update({"font.size": 10, "axes.edgecolor": muted, "axes.labelcolor": ink,
                         "xtick.color": muted, "ytick.color": muted, "axes.spines.top": False,
                         "axes.spines.right": False})
    out = HERE / "figures"
    out.mkdir(exist_ok=True)

    # Figure 1: speed trace with contraction / relaxation peaks, displacement trace below.
    fig, (a1, a2) = plt.subplots(2, 1, figsize=(8, 6), sharex=True)
    a1.plot(t_mid, speed, color=ink, lw=1.4)
    c = [p[0] for p in pairs]
    r = [p[2] for p in pairs]
    a1.plot(t_mid[c], speed[c], "o", ms=8, color=blue, label="Contraction (shortening) peak")
    a1.plot(t_mid[r], speed[r], "o", ms=8, color=orange, label="Relaxation (re-lengthening) peak")
    other = sorted(set(peaks(speed)) - set(c) - set(r))
    if other:
        a1.plot(t_mid[other], speed[other], "o", ms=8, mfc="none", mec=muted,
                label="Relaxation peak of the beat under way at frame 0")
    a1.set_ylabel("Optical-flow speed\n(px/frame, mean)")
    a1.legend(frameon=False, fontsize=8, loc="upper center", bbox_to_anchor=(0.5, -0.04), ncol=3)
    a1.set_title("Ex00628: two speed peaks per beat; one displacement peak per beat",
                 fontsize=11, loc="left", color=ink)
    a2.plot(t, disp, color=ink, lw=1.4)
    a2.plot(t[beat_idx], disp[beat_idx], "v", ms=8, color=blue, label="Peak contraction")
    a2.axvline(t[ref_idx], ymax=0.8, color=muted, lw=1, ls=":")
    a2.text(t[ref_idx], disp.max() * 0.35, " reference\n (relaxed) frame", fontsize=8, color=muted)
    for i in range(len(beat_idx) - 1):
        tm = 0.5 * (beat_t[i] + beat_t[i + 1])
        a2.text(tm, disp[beat_idx].max() * 1.02, f"{beat_t[i+1]-beat_t[i]:.3f} s",
                ha="center", fontsize=8, color=muted)
    a2.set_ylabel("Mean |frame - reference|\n(gray levels)")
    a2.set_xlabel("Time from first frame (s, MetaMorph per-plane timestamps)")
    a2.legend(frameon=False, fontsize=8, loc="lower right")
    fig.tight_layout()
    fig.savefig(out / "fig1_motion_trace.png", dpi=200)
    plt.close(fig)

    # Figure 2: spectrum of the consecutive-frame-difference trace, true vs nominal time base.
    fig, ax = plt.subplots(figsize=(8, 4))
    _, f_true, s_true = dominant_hz(framediff, fs_true)
    _, f_nom, s_nom = dominant_hz(framediff, NOMINAL_FPS)
    ax.plot(f_true, s_true / s_true.max(), color=ink, lw=1.6, label=f"Timestamps ({fs_true:.2f} frames/s)")
    ax.plot(f_nom, s_nom / s_true.max(), color=orange, lw=1.2, ls="--",
            label="Nominal 30 frames/s from metadata sheet")
    for h in (1, 2, 3):
        ax.axvline(h * PACING_HZ, color=blue, lw=0.8, ls=":")
        ax.text(h * PACING_HZ, 1.04, f"{h} x 1 Hz" if h > 1 else "1 Hz pacing",
                ha="center", fontsize=8, color=blue)
    ax.set_xlim(0, 5)
    ax.set_ylim(0, 1.12)
    ax.set_xlabel("Frequency (Hz)")
    ax.set_ylabel("Amplitude (normalized)")
    ax.set_title("Spectrum of consecutive |frame difference|: tallest peak is the 3rd harmonic",
                 fontsize=11, loc="left", color=ink)
    ax.legend(frameon=False, fontsize=8, loc="upper center", bbox_to_anchor=(0.5, -0.2), ncol=2)
    fig.tight_layout()
    fig.savefig(out / "fig2_spectrum.png", dpi=200)
    plt.close(fig)
    print(f"\nFigures written to {out}")


if __name__ == "__main__":
    main()
