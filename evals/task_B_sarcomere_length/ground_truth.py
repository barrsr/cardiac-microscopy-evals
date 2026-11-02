"""Ground truth and naive baselines for Task B: sarcomere length from ACTN2-mEGFP z-lines.

Run:  python ground_truth.py            (prints GT and every naive answer)
      python ground_truth.py --figures  (also writes figures/*.png)

Input: data/actn2_zstack.ome.tif (Z,Y,X uint16) and data/cell_mask.ome.tif (Y,X, 1 = target cell).
Pixel size is read from the OME-XML, never assumed.
"""
import re
import sys
from pathlib import Path

import numpy as np
import tifffile
from numpy.fft import fft2, fftfreq, fftshift, rfftfreq
from scipy import ndimage, signal, stats

HERE = Path(__file__).resolve().parent
BAND = (1.2, 3.0)  # physiological search band for the fundamental, micrometers


# ----------------------------------------------------------------------------- io
def load():
    with tifffile.TiffFile(HERE / "data" / "actn2_zstack.ome.tif") as t:
        stack = t.asarray().astype(float)
        ome = t.ome_metadata
    px = float(re.search(r'PhysicalSizeX="([0-9.]+)"', ome).group(1))
    pz = float(re.search(r'PhysicalSizeZ="([0-9.]+)"', ome).group(1))
    mask = tifffile.imread(HERE / "data" / "cell_mask.ome.tif") > 0
    return stack, mask, px, pz


def focus_plane(stack, mask):
    """In-focus plane = max variance of the Laplacian-of-Gaussian inside the cell."""
    s = [ndimage.gaussian_laplace(p, 1.0)[mask].var() for p in stack]
    return int(np.argmax(s)), np.array(s)


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


# -------------------------------------------------------------- naive approaches
def axis_spectrum_period(im, mask, px, axis, lo_um=0.6, hi_um=6.0):
    x = highpass(im, px)
    xm = np.where(mask, x - x[mask].mean(), 0)
    P = (np.abs(np.fft.rfft(xm, n=2048, axis=axis)) ** 2).mean(axis=1 - axis)
    f = rfftfreq(2048, px)
    sel = (f > 1 / hi_um) & (f < 1 / lo_um)
    return 1 / f[np.argmax(np.where(sel, P, 0))]


def axis_acf_first_peak(im, mask, px, axis):
    x = highpass(im, px)
    xm = np.where(mask, x - x[mask].mean(), 0)
    F = np.fft.rfft(xm, n=2048, axis=axis)
    ac = np.fft.irfft(np.abs(F) ** 2, n=2048, axis=axis).mean(axis=1 - axis)[:200]
    pk, _ = signal.find_peaks(ac / ac[0])
    return pk[0] * px


def axis_line_scan(im, mask, px, axis):
    cy, cx = (int(v) for v in ndimage.center_of_mass(mask))
    sm = ndimage.gaussian_filter(im, 0.8)
    line, m = (sm[cy], mask[cy]) if axis == 1 else (sm[:, cx], mask[:, cx])
    seg = line[m]
    seg = seg - ndimage.uniform_filter1d(seg, 25)
    pk, _ = signal.find_peaks(seg, distance=int(1.0 / px), prominence=0.5 * seg.std())
    return np.median(np.diff(pk)) * px, len(pk)


def raw_spectrum_global_max(im, mask, px, fmin):
    fc, prof, *_ = radial_power(im, mask, px, bgsub=False, soft=False)
    s = fc > fmin
    return 1 / fc[s][np.argmax(prof[s])]


# ------------------------------------------------------------------------- main
def main(figures=False):
    stack, mask, px, pz = load()
    zf, fscore = focus_plane(stack, mask)
    mip = stack.max(0)
    print(f"pixel size from OME-XML: {px} um (xy), {pz} um (z); stack {stack.shape}")
    print(f"in-focus plane (max LoG variance in mask): z={zf}; focus score z={zf}: {fscore[zf]:.0f}, "
          f"middle plane z={stack.shape[0] // 2}: {fscore[stack.shape[0] // 2]:.0f}")

    res = {}
    for name, im in [("MIP", mip), (f"plane z={zf}", stack[zf])]:
        fc, prof, *_ = radial_power(im, mask, px)
        p1, k = band_peak(fc, prof)
        f0 = 1 / p1
        pw = lambda m: prof[np.argmin(abs(fc - m * f0))]
        tl = tile_spectra(im, mask, px)
        tp = np.array([t["period"] for t in tl]); tpr = np.array([t["prom"] for t in tl])
        ta = np.array([t["ang"] for t in tl])
        org = tpr >= np.median(tpr)
        pw_med, d_all, angs, _ = aligned_profiles(im, mask, px)
        core = d_all[(d_all > 1.2) & (d_all < 3.0)]
        kde = stats.gaussian_kde(core, bw_method=0.05)
        xs = np.linspace(1.2, 3.0, 1801)
        mode = xs[np.argmax(kde(xs))]
        res[name] = dict(fft=p1, tile=np.median(tp[org]), prof=np.median(pw_med), mode=mode)
        print(f"\n[{name}]")
        print(f"  GT method 1a  whole-cell 2D power spectrum, radial peak in {BAND} um band: {p1:.3f} um ({p1 / px:.2f} px)")
        print(f"      harmonic check: P(f0/2)/P(f0)={pw(0.5) / pw(1):.3f}, P(2f0)/P(f0)={pw(2) / pw(1):.3f}, P(3f0)/P(f0)={pw(3) / pw(1):.4f}")
        print(f"  GT method 1b  128-px tiles, median of organized half ({org.sum()}/{len(tl)} tiles): {np.median(tp[org]):.3f} um; "
              f"IQR {np.percentile(tp[org], 25):.3f}-{np.percentile(tp[org], 75):.3f}; myofibril axis median {np.median(ta[org]):.1f} deg from +x (y down)")
        print(f"  GT method 2   structure-tensor-aligned z-line profiles: median of {len(pw_med)} window medians {np.median(pw_med):.3f} um "
              f"(IQR {np.percentile(pw_med, 25):.3f}-{np.percentile(pw_med, 75):.3f}); mode of {len(d_all)} intervals {mode:.3f} um; "
              f"axis median {np.median(angs):.1f} deg")

    gt = np.mean([res["MIP"]["fft"], res["MIP"]["prof"]])
    print(f"\nGROUND TRUTH (mean of method 1a and method 2 on the MIP): {gt:.3f} um -> reported as {round(gt, 2):.2f} um")

    # ---------------- naive approaches
    full = np.ones_like(mask)
    zm = stack.shape[0] // 2
    fcm, profm, *_ = radial_power(stack[zm], mask, px)
    fcu, profu, *_ = radial_power(mip, full, px)
    ax_angle = np.median([t["ang"] for t in tile_spectra(mip, mask, px)])
    naive = [
        ("1D spectrum along image x (rows), band 0.6-6 um", axis_spectrum_period(mip, mask, px, 1)),
        ("1D spectrum along image y (columns), band 0.6-6 um", axis_spectrum_period(mip, mask, px, 0)),
        ("1D autocorrelation along x, first peak", axis_acf_first_peak(mip, mask, px, 1)),
        ("1D autocorrelation along y, first peak", axis_acf_first_peak(mip, mask, px, 0)),
        ("single horizontal line scan through centroid, median peak spacing", axis_line_scan(mip, mask, px, 1)[0]),
        ("single vertical line scan through centroid, median peak spacing", axis_line_scan(mip, mask, px, 0)[0]),
        ("radial spectrum global max, no high-pass, f > 0", raw_spectrum_global_max(mip, mask, px, 0.0)),
        ("radial spectrum global max, no high-pass, f > 0.05 cyc/um", raw_spectrum_global_max(mip, mask, px, 0.05)),
        ("radial spectrum, fundamental band peak, reported in pixels", res["MIP"]["fft"] / px),
        ("second harmonic taken as the answer (radial peak in 0.6-1.2 um)", band_peak(*radial_power(mip, mask, px)[:2], band=(0.6, 1.2))[0]),
        (f"middle plane z={zm} only, radial band peak", band_peak(fcm, profm)[0]),
        ("whole crop with no cell mask, radial band peak", band_peak(fcu, profu)[0]),
        ("radial band peak with a nominal pixel size of 6.5 um / 52.29x = 0.1243 um", res["MIP"]["fft"] / px * 0.1243),
    ]
    print(f"\nexpected inflation if measured along x: GT/cos({ax_angle:.1f} deg) = {gt / np.cos(np.radians(ax_angle)):.3f} um")
    print("\nNAIVE APPROACHES")
    for name, v in naive:
        print(f"  {v:8.3f}   {name}")

    if figures:
        make_figures(stack, mask, px, zf, gt, res)


def make_figures(stack, mask, px, zf, gt, res):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.patches import Rectangle

    mip = stack.max(0)
    fig_dir = HERE / "figures"
    fig_dir.mkdir(exist_ok=True)
    lo, hi = np.percentile(mip, 1), np.percentile(mip[mask], 99.7)
    disp = np.clip((mip - lo) / (hi - lo), 0, 1)

    # ---- Figure 1: crop, aligned profile vs x-axis profile
    pw, d_all, angs, wins = aligned_profiles(mip, mask, px)
    L = int(12 / px)
    t = np.arange(-L // 2, L // 2 + 1, 0.25)
    sm = ndimage.gaussian_filter(mip, 0.8)

    def prof_stats(p):
        q = p - ndimage.uniform_filter1d(p, int(3 / px * 4))
        pk, _ = signal.find_peaks(q, distance=int(1 / px * 4), prominence=0.5 * q.std())
        dd = np.diff(pk) * 0.25 * px
        return (len(pk), dd.std() / dd.mean()) if len(dd) > 2 else (0, np.inf)

    # illustrate at the analysis window whose aligned profile is most regular (>= 6 z-lines, lowest CV)
    best = None
    for (y, x, th) in wins:
        p = ndimage.map_coordinates(sm, [y + t * np.sin(th), x + t * np.cos(th)], order=1)
        n, cv = prof_stats(p)
        if n >= 6 and (best is None or cv < best[0]):
            best = (cv, y, x, th)
    _, y, x, th = best
    p_al = ndimage.map_coordinates(sm, [y + t * np.sin(th), x + t * np.cos(th)], order=1)
    p_x = ndimage.map_coordinates(sm, [np.full_like(t, y), x + t], order=1)
    fig = plt.figure(figsize=(12, 5.2))
    ax0 = fig.add_axes([0.02, 0.08, 0.42, 0.86])
    ax0.imshow(disp, cmap="gray")
    ax0.contour(mask, [0.5], colors="#e6c229", linewidths=0.8)
    ax0.plot(x + t * np.cos(th), y + t * np.sin(th), color="#2a9d8f", lw=2, label="along myofibril")
    ax0.plot(x + t, np.full_like(t, y), color="#e76f51", lw=2, ls="--", label="along image x")
    ax0.plot([40, 40 + 10 / px], [960, 960], color="w", lw=3)
    ax0.text(40, 945, "10 µm", color="w", fontsize=9)
    ax0.set_axis_off()
    ax0.legend(loc="upper right", fontsize=8, framealpha=0.8)
    ax0.set_title("ACTN2-mEGFP, max projection; yellow = target cell", fontsize=10)
    ax1 = fig.add_axes([0.52, 0.56, 0.45, 0.36])
    ax2 = fig.add_axes([0.52, 0.10, 0.45, 0.36])
    for ax, p, col, lab in [(ax1, p_al, "#2a9d8f", "along myofibril"), (ax2, p_x, "#e76f51", "along image x")]:
        q = p - ndimage.uniform_filter1d(p, int(3 / px * 4))
        pk, _ = signal.find_peaks(q, distance=int(1 / px * 4), prominence=0.5 * q.std())
        ax.plot(t * px, p, color=col)
        ax.plot(t[pk] * px, p[pk], "k.", ms=6)
        sp = np.median(np.diff(t[pk])) * px
        ax.set_title(f"{lab}: median z-line spacing {sp:.2f} µm", fontsize=10)
        ax.set_ylabel("intensity")
        ax.spines[["top", "right"]].set_visible(False)
    ax2.set_xlabel("distance along line (µm)")
    fig.savefig(fig_dir / "fig1_profile_axis_vs_myofibril.png", dpi=150)
    plt.close(fig)

    # ---- Figure 2: 2D power spectrum and radial profile
    fc, prof, P, FX, FY = radial_power(mip, mask, px)
    p1, _ = band_peak(fc, prof)
    f0 = 1 / p1
    fig, (a, b) = plt.subplots(1, 2, figsize=(12, 5))
    lim = 1.4
    c = (np.abs(FX[0]) <= lim)
    sub = np.log10(P[np.ix_(c, c)] + 1)
    a.imshow(sub, extent=[-lim, lim, lim, -lim], cmap="magma", vmin=np.percentile(sub, 50), vmax=np.percentile(sub, 99.99))
    for m, col, lab in [(1, "#2a9d8f", f"fundamental {p1:.2f} µm"), (2, "#8ecae6", f"2nd harmonic {p1 / 2:.2f} µm")]:
        a.add_patch(plt.Circle((0, 0), m * f0, fill=False, ec=col, ls="--", lw=1.2, label=lab))
    a.legend(loc="lower left", fontsize=8)
    a.set_xlabel("fx (cycles/µm)")
    a.set_ylabel("fy (cycles/µm)")
    a.set_title("2D power spectrum of the target cell (log)", fontsize=10)
    s = (fc > 0.02) & (fc < 2.5)
    b.semilogy(fc[s], prof[s], color="#333")
    b.axvline(f0, color="#2a9d8f", ls="--", label=f"fundamental f0 = 1/{p1:.2f} µm")
    b.axvline(2 * f0, color="#8ecae6", ls="--", label=f"2 f0 = 1/{p1 / 2:.2f} µm")
    b.axvspan(1 / BAND[1], 1 / BAND[0], color="#2a9d8f", alpha=0.08, label="1.2-3.0 µm search band")
    b.set_xlabel("spatial frequency (cycles/µm)")
    b.set_ylabel("radially averaged power")
    b.set_title("Radial profile: fundamental vs harmonic", fontsize=10)
    b.legend(fontsize=8)
    b.spines[["top", "right"]].set_visible(False)
    fig.tight_layout()
    fig.savefig(fig_dir / "fig2_power_spectrum.png", dpi=150)
    plt.close(fig)
    print("figures written to", fig_dir)


if __name__ == "__main__":
    main(figures="--figures" in sys.argv)
