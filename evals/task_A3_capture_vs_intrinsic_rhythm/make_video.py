"""Looping side-by-side video of a captured tissue and one beating on its own rhythm, for task A3.

Usage:
    MPLBACKEND=Agg python make_video.py [data_dir]

Set FONT_DIR to a folder holding IBM Plex Sans, IBM Plex Mono and Instrument Serif TTFs to match the
results page; otherwise the system sans is used.

Left: Ex00567 (1:1 capture). Right: Ex00576 (intrinsic rhythm above the pacing rate). Beat times come from
ground_truth.analyze (contraction upstroke on the per-plane MetaMorph timestamps), so the
markers are exactly the beat times the ground truth uses. Frames play in real time from the
per-frame timestamps. Writes figures/beating_side_by_side.mp4 and a poster PNG of the last frame.
"""
import os
import shutil
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

FFMPEG = shutil.which("ffmpeg") or "ffmpeg"
OUT_FPS = 30
W, H = 960, 540                       # layout units; the video is rendered at 4/3 scale (1280x720)
DPI = 100 * 4 / 3
ROWS = (106, 406)                     # crop: the tissue between the pillars, full width
PAPER, INK, MUTED, FAINT, GRID = "#fcfcfa", "#0e0f12", "#6f737c", "#a3a7ae", "#dedcd6"
PANELS = [
    dict(id="Ex00567", color="#13774f", label="Follows the pulse",
         sub="A beat every 1.00 s, in step with the pulse"),
    dict(id="Ex00576", color="#7a1f2b", label="Keeps its own rhythm",
         sub="A beat every 0.96 s, pulling ahead of the pulse"),
]

SANS, SERIF, MONO = ["Arial", "Helvetica", "DejaVu Sans"], ["Georgia", "DejaVu Serif"], ["Menlo", "DejaVu Sans Mono"]
if os.environ.get("FONT_DIR"):
    from matplotlib import font_manager
    for f in Path(os.environ["FONT_DIR"]).glob("*.ttf"):
        font_manager.fontManager.addfont(str(f))
    SANS, SERIF, MONO = ["IBM Plex Sans"] + SANS, ["Instrument Serif"] + SERIF, ["IBM Plex Mono"] + MONO
plt.rcParams.update({"font.family": SANS, "font.size": 10, "text.color": INK})


def when(ms):
    """Offset from the expected beat time, in words."""
    if abs(ms) <= 10:
        return "on time"
    return f"{abs(ms)} ms {'early' if ms < 0 else 'late'}"


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
    fig = plt.figure(figsize=(W / 100, H / 100), dpi=DPI, facecolor=PAPER)
    px = lambda x, y, w, h: [x / W, 1 - (y + h) / H, w / W, h / H]   # top-left rects in layout units
    fig.text(24 / W, 1 - 28 / H, "Two engineered heart tissues, both paced once a second",
             fontsize=19, family=SERIF, color=INK, va="center")
    clock = fig.text(936 / W, 1 - 28 / H, "", fontsize=10, color=MUTED, ha="right", va="center",
                     family=MONO)
    pw, x0s = 444, (24, 492)
    ih = int(round((ROWS[1] - ROWS[0]) * pw / 512))
    dyn = []
    for p, x0 in zip(panels, x0s):
        fig.text(x0 / W, 1 - 64 / H, p["label"], fontsize=15, family=SERIF, color=p["color"], va="center")
        fig.text((x0 + pw) / W, 1 - 64 / H, p["id"], fontsize=8, family=MONO, color=FAINT, ha="right",
                 va="center")
        fig.text(x0 / W, 1 - 84 / H, p["sub"], fontsize=9, color=MUTED, va="center")
        ax = fig.add_axes(px(x0, 98, pw, ih))
        ax.set_axis_off()
        p["small"] = [cv2.resize(f, (pw, ih), interpolation=cv2.INTER_AREA) for f in p["frames"]]
        im = ax.imshow(p["small"][0], cmap="gray", vmin=0, vmax=255, interpolation="nearest",
                       aspect="auto")
        # Timeline: rings where a beat is due if the tissue follows the pulse, dots where it beat
        tl = fig.add_axes(px(x0, 98 + ih + 14, pw, 62))
        tl.set_facecolor(PAPER)
        tl.set_xlim(0, dur)
        tl.set_ylim(-1.15, 0.9)
        for sp in tl.spines.values():
            sp.set_visible(False)
        tl.set_xticks([])
        tl.set_yticks([])
        tl.axhline(0.15, color=GRID, lw=1, zorder=0)
        grid = p["bt"][0] + PACING_S * np.arange(10)
        grid = grid[grid <= p["t"][-1]]
        for g in grid:
            tl.plot([g], [0.15], "o", ms=10, mfc=PAPER, mec=FAINT, mew=1.3, zorder=2)
        marks = []
        for k, (b, ph) in enumerate(zip(p["bt"], p["phase"])):
            g = b - ph
            bar, = tl.plot([g, b], [0.15, 0.15], color=p["color"], lw=2.2, alpha=0.55, zorder=3)
            dot, = tl.plot([b], [0.15], "o", ms=7, color=p["color"], zorder=4)
            txt = tl.text(b, -0.4, when(int(round(ph * 1000))) if k else "first beat", fontsize=8, color=p["color"],
                          ha="center", va="top")
            marks.append((b, (bar, dot, txt)))
        cur = tl.axvline(0, ymin=0.3, ymax=0.85, color=INK, lw=1, zorder=5)
        dyn.append((p, im, marks, cur))
    # Key, drawn with the same marks as the timelines
    key = fig.add_axes(px(24, 466, 912, 18))
    key.set_axis_off()
    key.set_xlim(0, 912)
    key.set_ylim(-1, 1)
    key.plot([6], [0], "o", ms=9, mfc=PAPER, mec=FAINT, mew=1.3)
    key.text(18, 0, "when a beat is due if the tissue follows the pulse", fontsize=8.5, color=MUTED,
             va="center")
    key.plot([300], [0], "o", ms=7, color=INK)
    key.text(312, 0, "when the tissue actually beat", fontsize=8.5, color=MUTED, va="center")
    fig.text(24 / W, 1 - 500 / H,
             "Brightfield, 4x, shown in real time from the per-frame timestamps. Beat time: contraction upstroke, "
             "from the task's ground-truth code.", fontsize=7.5, color=FAINT, va="center")
    fig.text(24 / W, 1 - 514 / H,
             "Rings: a 1.000 s grid from each tissue's first beat; the pulse itself is not recorded. "
             "Data: Kobeissi, Lejeune et al., Zenodo 11120809, CC0.", fontsize=7.5, color=FAINT, va="center")
    return fig, clock, dyn


def render(fig, clock, dyn, tau):
    clock.set_text(f"{tau:4.2f} s")
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
           "-s", f"{round(W * DPI / 100)}x{round(H * DPI / 100)}", "-r", str(OUT_FPS), "-i", "-", "-an", "-c:v", "libx264",
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
