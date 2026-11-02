"""Ground truth and naive baselines for Task B2: which cells are striated, and their median sarcomere length.

Run:  python ground_truth.py            (per-cell table, GT, and every naive answer)
      python ground_truth.py --figures  (also writes figures/fig1_striated_vs_immature.png)

Inputs: data/cell_<X>_actn2.ome.tif (YX uint16, maximum projection of ACTN2-mEGFP) and
data/cell_<X>_mask.ome.tif (YX, 1 = target cell). Pixel size is read from each file's OME-XML.
Author annotations for the panel (expert scores, classifier fractions, published peak distance)
come from provenance/panel_author_annotations.csv, written by provenance/build_snapshot.py.

The spectral and profile methods below are the Task B v1 methods, unchanged.
"""
import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import tifffile
from numpy.fft import fft2, fftfreq, fftshift
from scipy import ndimage, signal

HERE = Path(__file__).resolve().parent
BAND = (1.2, 3.0)          # search band for the sarcomere fundamental, micrometers
CELLS = [f"cell_{c}" for c in "ABCDEFGH"]

# Striation rule from the authors' annotations (see eval notes, GROUND TRUTH):
#   striated     = both expert scores >= 4 and classifier organized z-disk fraction >= 0.5
#   not striated = both expert scores <= 2 and classifier organized z-disk fraction <= 0.05
# Score-3 cells (mixed puncta, fibers and z-disks) were excluded from the panel on purpose.
# Image confirmation: local spectral prominence of the in-band peak (power at f0 over the
# geometric mean of power at 0.75 f0 and 1.33 f0) >= 4 for striated, <= 2 for not striated.


# ----------------------------------------------------------------------------- io
def load(cid):
    with tifffile.TiffFile(HERE / "data" / f"{cid}_actn2.ome.tif") as t:
        im = t.asarray().astype(float)
        ome = t.ome_metadata
    px = float(re.search(r'PhysicalSizeX="([0-9.]+)"', ome).group(1))
    mask = tifffile.imread(HERE / "data" / f"{cid}_mask.ome.tif") > 0
    return im, mask, px


def highpass(im, px, sigma_um=3.0):
    return im - ndimage.gaussian_filter(im, sigma_um / px)


# ------------------------------------------------- GT method 1: 2D power spectrum
def radial_power(im, mask, px, pad=2048, bgsub=True, soft=True):
    x = highpass(im, px) if bgsub else im.copy()
    x = (x - x[mask].mean()) * mask
    if soft:
        x = x * ndimage.gaussian_filter(mask.astype(float), 5)
    P = np.abs(fftshift(fft2(x, s=(pad, pad)))) ** 2
    f = fftshift(fftfreq(pad, px))
    FX, FY = np.meshgrid(f, f)
    R = np.hypot(FX, FY)
    df = f[1] - f[0]
    edges = np.arange(0, R.max(), df)
    idx = np.digitize(R.ravel(), edges)
    n = len(edges)
    prof = np.bincount(idx, P.ravel(), minlength=n + 1)[1:n] / np.maximum(np.bincount(idx, minlength=n + 1)[1:n], 1)
    fc = edges[:-1] + df / 2
    return fc, prof, P, FX, FY


def band_peak(fc, prof, band=BAND):
    sel = np.where((fc > 1 / band[1]) & (fc < 1 / band[0]))[0]
    k = sel[np.argmax(prof[sel])]
    a, b, c = np.log(prof[k - 1: k + 2])
    d = 0.5 * (a - c) / (a - 2 * b + c)  # parabolic refinement on log power
    f0 = fc[k] + d * (fc[1] - fc[0])
    return 1 / f0, k


def tile_spectra(im, mask, px, size=128, step=64, pad=1024, band=BAND):
    win = np.outer(np.hanning(size), np.hanning(size))
    f = fftshift(fftfreq(pad, px))
    FX, FY = np.meshgrid(f, f)
    R = np.hypot(FX, FY)
    ann = (R > 1 / band[1]) & (R < 1 / band[0])
    out = []
    H, W = mask.shape
    for y in range(0, H - size + 1, step):
        for x in range(0, W - size + 1, step):
            if mask[y:y + size, x:x + size].mean() < 0.95:
                continue
            t = im[y:y + size, x:x + size]
            P = np.abs(fftshift(fft2((t - t.mean()) * win, s=(pad, pad)))) ** 2
            i = np.argmax(np.where(ann, P, 0))
            out.append(dict(y=y, x=x, period=1 / R.flat[i], prom=P.flat[i] / np.median(P[ann]),
                            ang=np.degrees(np.arctan2(FY.flat[i], FX.flat[i])) % 180))
    return out


# ------------------------------- GT method 2: orientation-aligned z-line profiles
def structure_tensor_axis(im, px, sigma_int_um=4.0):
    g = ndimage.gaussian_filter(im, 1.0)
    gy, gx = ndimage.sobel(g, 0), ndimage.sobel(g, 1)
    s = sigma_int_um / px
    Jxx, Jyy, Jxy = (ndimage.gaussian_filter(v, s) for v in (gx * gx, gy * gy, gx * gy))
    theta = 0.5 * np.arctan2(2 * Jxy, Jxx - Jyy)  # dominant gradient = across z-lines = along the myofibril
    coh = np.sqrt((Jxx - Jyy) ** 2 + 4 * Jxy ** 2) / (Jxx + Jyy + 1e-12)
    return theta, coh


def aligned_profiles(im, mask, px, L_um=12.0, n_lines=9, sep_um=0.6, grid=40, sub=4):
    theta, coh = structure_tensor_axis(im, px)
    sm = ndimage.gaussian_filter(im, 0.8)
    L = int(L_um / px)
    t = np.arange(-L // 2, L // 2 + 1, 1 / sub)
    thr = np.percentile(coh[mask], 50)  # keep the more coherent (organized) half of the cell
    per_window, all_d, angs, windows = [], [], [], []
    for y in range(grid, im.shape[0] - grid, grid):
        for x in range(grid, im.shape[1] - grid, grid):
            if not mask[y, x] or coh[y, x] < thr:
                continue
            th = theta[y, x]
            ux, uy = np.cos(th), np.sin(th)
            d_win = []
            nprof = 0
            for k in range(n_lines):
                off = (k - n_lines // 2) * sep_um / px
                xs, ys = x - off * uy + t * ux, y + off * ux + t * uy
                if xs.min() < 0 or ys.min() < 0 or xs.max() > im.shape[1] - 1 or ys.max() > im.shape[0] - 1:
                    continue
                if not mask[ys.astype(int), xs.astype(int)].all():
                    continue
                p = ndimage.map_coordinates(sm, [ys, xs], order=1)
                p = p - ndimage.uniform_filter1d(p, int(3.0 / px * sub))
                pk, _ = signal.find_peaks(p, distance=int(1.0 / px * sub), prominence=0.5 * p.std())
                pos = []
                for q in pk:
                    if 0 < q < len(p) - 1:
                        a, b, c = p[q - 1], p[q], p[q + 1]
                        den = a - 2 * b + c
                        pos.append(q + (0.5 * (a - c) / den if den != 0 else 0.0))
                d_win += list(np.diff(pos) / sub * px)
                nprof += 1
            if nprof < n_lines // 2 or len(d_win) < 5:
                continue
            d_win = np.array(d_win)
            per_window.append(np.median(d_win))
            all_d += d_win.tolist()
            angs.append(np.degrees(th) % 180)
            windows.append((y, x, th))
    return np.array(per_window), np.array(all_d), np.array(angs), windows




def local_prominence(fc, prof, f0):
    pw = lambda f: prof[np.argmin(abs(fc - f))]
    return pw(f0) / np.sqrt(pw(0.75 * f0) * pw(1.33 * f0))


def measure(im, mask, px):
    """All v1 measurements for one cell, inside its mask."""
    fc, prof, *_ = radial_power(im, mask, px)
    p1, _ = band_peak(fc, prof)
    f0 = 1 / p1
    pw = lambda m: prof[np.argmin(abs(fc - m * f0))]
    tl = tile_spectra(im, mask, px)
    tp = np.array([t["period"] for t in tl])
    tpr = np.array([t["prom"] for t in tl])
    tile = np.median(tp[tpr >= np.median(tpr)]) if len(tl) else np.nan
    pwm, d_all, _, _ = aligned_profiles(im, mask, px)
    return dict(fft=p1, prom=local_prominence(fc, prof, f0), h2=pw(2) / pw(1), tile=tile, ntile=len(tl),
                prof=np.median(pwm) if len(pwm) else np.nan, nwin=len(pwm), nint=len(d_all))


def main(figures=False):
    ann = pd.read_csv(HERE / "provenance" / "panel_author_annotations.csv").set_index("cell_id")
    rows = []
    for cid in CELLS:
        im, mask, px = load(cid)
        r = measure(im, mask, px)
        a = ann.loc[cid]
        e_lo, e_hi = min(a.expert_mh, a.expert_kg), max(a.expert_mh, a.expert_kg)
        if e_lo >= 4 and a.frac_organized_zdisks >= 0.5:
            author_call = "striated"
        elif e_hi <= 2 and a.frac_organized_zdisks <= 0.05:
            author_call = "not striated"
        else:
            author_call = "AMBIGUOUS"
        image_call = "striated" if r["prom"] >= 4 else ("not striated" if r["prom"] <= 2 else "AMBIGUOUS")
        # unmasked crop, for the neighbor-contamination baseline
        fcu, profu, *_ = radial_power(im, np.ones_like(mask), px)
        pu, _ = band_peak(fcu, profu)
        r.update(cell=cid, px=px, mh=a.expert_mh, kg=a.expert_kg, zd=a.frac_organized_zdisks,
                 op=a.frac_organized_puncta, dp=a.frac_disorganized_puncta, fib=a.frac_fibers,
                 author_pk=a.author_peak_distance_um, author_call=author_call, image_call=image_call,
                 unmasked_fft=pu, unmasked_prom=local_prominence(fcu, profu, 1 / pu))
        r["gt_sl"] = np.mean([r["fft"], r["prof"]]) if author_call == "striated" else np.nan
        rows.append(r)
    d = pd.DataFrame(rows).set_index("cell")
    assert (d.author_call == d.image_call).all(), "author annotation and image check disagree"
    assert not (d.author_call == "AMBIGUOUS").any()

    print(f"pixel size from OME-XML: {sorted(set(d.px))} um")
    print("\nPER-CELL TABLE (methods inside the cell mask; lengths in um)")
    cols = ["mh", "kg", "zd", "op", "fib", "dp", "author_pk", "fft", "prom", "h2", "tile", "prof", "nwin",
            "unmasked_fft", "unmasked_prom", "author_call", "gt_sl"]
    with pd.option_context("display.width", 250, "display.max_columns", 30):
        print(d[cols].round(3).to_string())

    S = d[d.author_call == "striated"]
    ids = sorted(S.index)
    med = float(np.median(S.gt_sl))
    print(f"\nGROUND TRUTH: striated cells {ids} (n = {len(ids)}); per-cell SL = mean(method 1a, method 2): "
          + ", ".join(f"{c} {v:.3f}" for c, v in S.gt_sl.items()))
    print(f"  median sarcomere length across striated cells: {med:.3f} um -> reported as {round(med, 2):.2f}")
    print(f"  method 1a only: median {np.median(S.fft):.3f}; method 1b only: {np.median(S.tile):.3f}; "
          f"method 2 only: {np.median(S.prof):.3f}")
    print(f"  separation: striated prominence {S.prom.min():.2f}-{S.prom.max():.2f}; "
          f"non-striated {d[d.author_call != 'striated'].prom.min():.2f}-{d[d.author_call != 'striated'].prom.max():.2f}")

    # ---------------- naive approaches
    N = d[d.author_call != "striated"]
    naive = [
        ("2D FFT band peak (1.2-3.0 um) on every cell, all called striated", len(d), np.median(d.fft)),
        ("tile 2D FFT on every cell, all called striated", len(d), np.median(d.tile)),
        ("aligned z-line profiles on every cell, all called striated", len(d), np.median(d.prof)),
        ("mean of 1a and 2 on every cell (GT method without the striation call)", len(d), np.median((d.fft + d.prof) / 2)),
        ("2D FFT on the unmasked crops, all called striated", len(d), np.median(d.unmasked_fft)),
        ("authors' published GLCM peak distance, all cells", len(d), np.median(d.author_pk)),
        ("authors' GLCM peak distance, striated cells only", len(S), np.median(S.author_pk)),
        ("striated = FFT band peak lands inside 1.5-2.2 um (masked)", int(d.fft.between(1.5, 2.2).sum()),
         np.median(d.fft[d.fft.between(1.5, 2.2)])),
        ("striated = FFT band peak inside 1.5-2.2 um, unmasked crops", int(d.unmasked_fft.between(1.5, 2.2).sum()),
         np.median(d.unmasked_fft[d.unmasked_fft.between(1.5, 2.2)])),
        ("striated = aligned-profile spacing inside 1.5-2.5 um", int(d.prof.between(1.5, 2.5).sum()),
         np.median(d.prof[d.prof.between(1.5, 2.5)])),
    ]
    print("\nNAIVE APPROACHES  (n_striated, median um)")
    for name, n, v in naive:
        print(f"  n={n}  {v:6.3f}   {name}")
    print("\nNon-striated cells, what each method returns: "
          + "; ".join(f"{c}: FFT {r.fft:.2f}, tiles {r.tile:.2f}, profiles {r.prof:.2f}, unmasked {r.unmasked_fft:.2f}, GLCM {r.author_pk:.2f}"
                      for c, r in N.iterrows()))

    if figures:
        make_figure(d)


def make_figure(d):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    order = list(d[d.author_call == "striated"].index) + list(d[d.author_call != "striated"].index)
    fig, ax = plt.subplots(2, len(order), figsize=(2.6 * len(order), 5.6), gridspec_kw=dict(height_ratios=[1, 0.8]))
    W = 15.0  # um window
    for j, cid in enumerate(order):
        im, mask, px = load(cid)
        h = int(W / px / 2)
        # window centered on the most in-mask region: maximum of the smoothed mask
        cy, cx = np.unravel_index(np.argmax(ndimage.uniform_filter(mask.astype(float), 2 * h)), mask.shape)
        cy, cx = max(cy, h), max(cx, h)
        crop = im[cy - h:cy + h, cx - h:cx + h]
        lo, hi = np.percentile(im[mask], [1, 99.7])
        a = ax[0, j]
        a.imshow(np.clip((crop - lo) / (hi - lo), 0, 1), cmap="gray")
        a.contour(mask[cy - h:cy + h, cx - h:cx + h], [0.5], colors="#e6c229", linewidths=0.8)
        a.plot([8, 8 + 5 / px], [2 * h - 10] * 2, color="w", lw=2.5)
        r = d.loc[cid]
        col = "#2a9d8f" if r.author_call == "striated" else "#e76f51"
        a.set_title(f"{cid}  experts {r.mh:.0f}/{r.kg:.0f}\nz-disk class {r.zd:.2f}", fontsize=9, color=col)
        a.set_axis_off()
        fc, prof, *_ = radial_power(im, mask, px)
        s = (fc > 0.2) & (fc < 1.5)
        b = ax[1, j]
        b.semilogy(fc[s], prof[s] / np.median(prof[s]), color=col, lw=1.2)
        b.axvspan(1 / BAND[1], 1 / BAND[0], color="0.85", zorder=0)
        b.axvline(1 / r.fft, color="k", ls=":", lw=1)
        b.set_title(f"band peak {r.fft:.2f} µm\nprominence {r.prom:.1f}", fontsize=8)
        b.set_xlabel("cycles/µm", fontsize=8)
        b.tick_params(labelsize=7)
        b.spines[["top", "right"]].set_visible(False)
        if j == 0:
            b.set_ylabel("radial power (rel.)", fontsize=8)
    fig.suptitle("ACTN2-mEGFP, day 32: striated cells (teal) vs puncta/fiber cells without z-disks (orange). "
                 "Scale bar 5 µm; gray band = 1.2-3.0 µm search band; dotted = band maximum.", fontsize=9)
    fig.tight_layout(rect=(0, 0, 1, 0.95))
    (HERE / "figures").mkdir(exist_ok=True)
    fig.savefig(HERE / "figures" / "fig1_striated_vs_immature.png", dpi=150)
    plt.close(fig)
    print("figure written to", HERE / "figures")


if __name__ == "__main__":
    main(figures="--figures" in sys.argv)
