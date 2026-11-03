"""Looping side-by-side video of a captured and an escaping tissue for task A3.

Usage:
    MPLBACKEND=Agg python make_video.py [data_dir]

Left: Ex00567 (1:1 capture). Right: Ex00576 (escape rhythm). Beat times come from
ground_truth.analyze (contraction upstroke on the per-plane MetaMorph timestamps), so the
markers are exactly the beat times the ground truth uses. Frames play in real time from the
per-frame timestamps. Writes figures/beating_side_by_side.mp4 and a poster PNG of the last frame.
"""
import subprocess
import sys
from pathlib import Path

import cv2
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from ground_truth import PACING_S, analyze, load  # noqa: E402

FFMPEG = "/opt/homebrew/bin/ffmpeg"
OUT_FPS = 30
W, H = 960, 540
DPI = 100
ROWS = (106, 406)                     # crop: the tissue between the pillars, full width
INK, MUTED, FAINT, GRID = "#1a1a1a", "#5f5f5f", "#9a9a9a", "#d4d4d4"
PANELS = [
    dict(id="Ex00567", color="#1f8a5b", label="captured: every beat on the 1.000 s period"),
    dict(id="Ex00576", color="#c23b2a", label="escaping: beats every 0.96 s, ahead of the stimulus"),
]

plt.rcParams.update({"font.family": ["Arial", "Helvetica", "DejaVu Sans"], "font.size": 10,
                     "text.color": INK})


def prepare(ddir, p):
    path = ddir / f"{p['id']}.tif"
    rec = analyze(path)
    stack, t = load(path)
    crop = stack[:, ROWS[0]:ROWS[1], :]
    lo, hi = np.percentile(crop, [0.5, 99.5])            # one stretch for the whole recording
    u8 = (np.clip((crop - lo) / (hi - lo), 0, 1) * 255).astype(np.uint8)
    p.update(t=t, frames=u8, bt=rec["bt"], phase=rec["phase"], ibi=rec["ibi"])
    return p


def build(panels, dur):
    fig = plt.figure(figsize=(W / DPI, H / DPI), dpi=DPI, facecolor="white")
    px = lambda x, y, w, h: [x / W, 1 - (y + h) / H, w / W, h / H]   # top-left pixel rects
    fig.text(24 / W, 1 - 26 / H, "Engineered heart tissues under 1 Hz field stimulation",
             fontsize=13, color=INK, weight="bold", va="center")
    clock = fig.text(936 / W, 1 - 26 / H, "", fontsize=11, color=MUTED, ha="right", va="center",
                     family="DejaVu Sans Mono")
    pw, x0s = 444, (24, 492)
    ih = int(round((ROWS[1] - ROWS[0]) * pw / 512))
    dyn = []
    for p, x0 in zip(panels, x0s):
        fig.text(x0 / W, 1 - 58 / H, p["id"], fontsize=12, weight="bold", color=INK, va="center")
        fig.text(x0 / W, 1 - 78 / H, p["label"], fontsize=10.5, color=p["color"], va="center")
        ax = fig.add_axes(px(x0, 92, pw, ih))
        ax.set_axis_off()
        p["small"] = [cv2.resize(f, (pw, ih), interpolation=cv2.INTER_AREA) for f in p["frames"]]
        im = ax.imshow(p["small"][0], cmap="gray", vmin=0, vmax=255, interpolation="nearest",
                       aspect="auto")
        # Timeline strip
        tl = fig.add_axes(px(x0, 92 + ih + 16, pw, 74))
        tl.set_xlim(0, dur)
        tl.set_ylim(-1.15, 1.45)
        for s in ("top", "right", "left"):
            tl.spines[s].set_visible(False)
        tl.spines["bottom"].set_color(FAINT)
        tl.set_yticks([])
        grid = p["bt"][0] + PACING_S * np.arange(10)
        grid = grid[grid <= p["t"][-1]]
        for g in grid:
            tl.plot([g, g], [-0.55, 0.85], color=GRID, lw=1.2, zorder=1)
            tl.plot([g], [0.15], "o", ms=9, mfc="white", mec=FAINT, mew=1.3, zorder=2)
        tl.set_xticks(grid)
        tl.set_xticklabels([f"{g:.2f}" for g in grid], fontsize=8.5, color=MUTED)
        tl.tick_params(axis="x", length=0, pad=4)
        tl.text(0, 1.45, "gray rings: expected stimulus-locked beat times (s)", fontsize=8, color=FAINT,
                va="top")
        marks = []
        for b, ph in zip(p["bt"], p["phase"]):
            g = b - ph
            bar, = tl.plot([g, b], [0.15, 0.15], color=p["color"], lw=2, alpha=0.6, zorder=3)
            dot, = tl.plot([b], [0.15], "o", ms=7, color=p["color"], zorder=4)
            ms = int(round(ph * 1000))
            txt = tl.text(b, -0.35, f"{ms:+d} ms".replace("-", "−") if ms else "0 ms",
                          fontsize=8.5, color=p["color"], ha="center", va="top")
            marks.append((b, (bar, dot, txt)))
        cur = tl.axvline(0, ymin=0, ymax=0.8, color=INK, lw=1.2, zorder=5)
        dyn.append((p, im, marks, cur))
    cap = ["Brightfield, 4x, played in real time from the per-frame timestamps. Dots: contraction upstrokes "
           "from the task's ground-truth code.",
           "Gray rings: a 1.000 s grid anchored at each tissue's first upstroke. The stimulus itself is not "
           "recorded, so the grid is anchored to the first beat.",
           "A beat cut off by either end of the recording is not timed."]
    for k, line in enumerate(cap):
        fig.text(24 / W, 1 - (474 + 16 * k) / H, line, fontsize=8.5, color=MUTED, va="center")
    fig.text(24 / W, 1 - 526 / H, "Data: Kobeissi, Lejeune et al., Zenodo 11120809, CC0",
             fontsize=7.5, color=FAINT, va="center")
    return fig, clock, dyn


def render(fig, clock, dyn, tau):
    clock.set_text(f"t = {tau:4.2f} s")
    for p, im, marks, cur in dyn:
        i = int(np.searchsorted(p["t"], tau, side="right") - 1)
        im.set_data(p["small"][max(0, min(i, len(p["t"]) - 1))])
        for b, arts in marks:
            for a in arts:
                a.set_visible(tau >= b)
        cur.set_xdata([tau, tau])
    fig.canvas.draw()
    return np.asarray(fig.canvas.buffer_rgba())[:, :, :3].copy()


def main():
    ddir = Path(sys.argv[1]) if len(sys.argv) > 1 else HERE / "data"
    panels = [prepare(ddir, dict(p)) for p in PANELS]
    # Common loop length: the longer recording plus one frame; the shorter holds its last frame.
    dur = max(p["t"][-1] + np.median(np.diff(p["t"])) for p in panels)
    for p in panels:
        print(p["id"], "beat times (s):", ", ".join(f"{b:.3f}" for b in p["bt"]),
              "| intervals:", ", ".join(f"{v:.3f}" for v in p["ibi"]),
              "| offset from grid (s):", ", ".join(f"{v:+.3f}" for v in p["phase"]))
    fig, clock, dyn = build(panels, dur)
    out = HERE / "figures"
    out.mkdir(exist_ok=True)
    mp4 = out / "beating_side_by_side.mp4"
    n_out = int(np.ceil(dur * OUT_FPS))
    cmd = [FFMPEG, "-y", "-loglevel", "error", "-f", "rawvideo", "-pix_fmt", "rgb24",
           "-s", f"{W}x{H}", "-r", str(OUT_FPS), "-i", "-", "-an", "-c:v", "libx264",
           "-preset", "veryslow", "-crf", "18", "-pix_fmt", "yuv420p", "-movflags", "+faststart",
           str(mp4)]
    proc = subprocess.Popen(cmd, stdin=subprocess.PIPE)
    last = None
    for k in range(n_out):
        last = render(fig, clock, dyn, min(k / OUT_FPS, dur))
        proc.stdin.write(last.tobytes())
    proc.stdin.close()
    proc.wait()
    # Poster: the last frame, with every marker shown.
    last = render(fig, clock, dyn, dur)
    plt.imsave(out / "beating_side_by_side_poster.png", last)
    print(f"{mp4} ({mp4.stat().st_size / 1e6:.2f} MB), {n_out} frames at {OUT_FPS} fps, {dur:.3f} s")


if __name__ == "__main__":
    main()
