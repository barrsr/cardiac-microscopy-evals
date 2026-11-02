#!/usr/bin/env python3
"""Fetch and rebuild the task B2 data snapshot (which cells are striated, and their sarcomere length).

Source: Gerbin KA, et al. "Cell states beyond transcriptomics: Integrating structural
organization and gene expression in hiPSC-derived cardiomyocytes." Cell Systems 12(6):670-687
(2021), doi:10.1016/j.cels.2021.05.001. Data from the Allen Institute for Cell Science public
bucket s3://allencell/aics/integrated_transcriptomics_structural_organization_hipsc_cm/,
read anonymously over plain HTTPS, so no AWS tools or credentials are needed.

License: Allen Institute for Cell Science terms of use
(https://www.allencell.org/terms-of-use.html), research and noncommercial use with
attribution. This repository does not redistribute the data; this script fetches it from the
Allen Institute directly.

What is rebuilt. Eight day-32 cells from one imaging session (20190816_M02_001, plate
5500000014, wells B2, B4 and B5). For each cell the script downloads the authors' 2D field
(2d_segmented_fields/2D_fov_tiff_path/<...>_annotations_corrected.tiff, 10 channels) and takes
  channel 1 ("488"): the authors' maximum projection of ACTN2-mEGFP, uint16
  channel 9 ("cell"): the authors' cell segmentation
(channel map in 2d_segmented_fields/supporting_files/channel_defs.json). The target cell is the
segmentation label listed in PANEL. The crop is the cell's bounding box grown by an 8 um margin
(round(8 / 0.1195257) = 67 px) on each side and clipped to the field. Pixel values are unchanged.
  data/cell_<X>_actn2.ome.tif  the crop of channel 1, uint16, YX
  data/cell_<X>_mask.ome.tif   the same crop of channel 9 with the target label set to 1, uint8
Both carry PhysicalSizeX/Y = 0.1195257 um, the SlideBook calibration (tag 65007 and TIFF
XResolution) of the session's raw files. The panel table, crop windows and author annotations
are those of the original build (provenance/panel_author_annotations.csv);
the script checks that the recomputed crop windows and mask sizes match them.

Intermediate downloads (eight fields, about 136 MB) go to the cache directory, default
~/.cache/cardiac-microscopy-evals, overridable with the CME_CACHE environment variable.
Pass --clean-cache to delete them after a successful run.

Verification. The OME-TIFFs are written with zlib compression and with the OME UUID and
Creator string of the original snapshot fixed, so with tifffile 2026.10.15 and the same zlib
the output is byte-identical to the original and is checked against the file MD5s below.
Other tifffile or zlib versions can lay out or compress the same pixels into different bytes.
When the file MD5 does not match, the script falls back to checking the decoded pixel array
(shape, dtype and MD5 of the raw pixel bytes) plus the physical pixel size in the OME-XML, and
reports the file as "pixel-identical". Either result means the data the model sees is the same.

The script is idempotent: files that already verify are left alone, and only the fields of
cells with a missing or wrong file are downloaded.

usage: python fetch_data.py [--clean-cache]

Requires: numpy, tifffile, requests (plus the standard library).
"""
import hashlib
import os
import re
import sys
import time
import urllib.parse
from pathlib import Path

import numpy as np
import requests
import tifffile

HERE = Path(__file__).resolve().parent
DATA = HERE / "data"
CACHE = Path(os.environ.get("CME_CACHE", "~/.cache/cardiac-microscopy-evals")).expanduser() / "allencell"
BUCKET = "https://allencell.s3.amazonaws.com/aics/integrated_transcriptomics_structural_organization_hipsc_cm/"
FIELD_DIR = "2d_segmented_fields/2D_fov_tiff_path/"

PX = 0.1195257  # um/px, SlideBook calibration of session 20190816_M02_001
MARGIN_UM = 8.0
CH_IMAGE, CH_CELL = 1, 9  # channel_defs.json: '488' = 1, 'cell' = 9
CREATOR = "tifffile.py 2026.10.15"

# cell id -> (FOV, 2D field file, field bytes, segmentation label, crop rows, crop cols, mask pixels)
# Panel ids were assigned with numpy default_rng(20261102).permutation(8) in the original build.
PANEL = {
    "cell_A": ("20190816_M02_001_Capture 3 - Capture 3 - B4 [144] Montage - 13",
               "c204b1dc_5500000014_63X_20190816_S1_P13_B4_annotations_corrected.tiff", 16_381_693, 1, (0, 775), (97, 795), 169_424),
    "cell_B": ("20190816_M02_001_Capture 3 - Capture 3 - B2 [144] Montage - 17",
               "c395b4de_5500000014_63X_20190816_S1_P17_B2_annotations_corrected.tiff", 16_399_078, 4, (45, 943), (509, 1331), 195_096),
    "cell_C": ("20190816_M02_001_Capture 3 - Capture 3 - B5 [144] Montage - 2",
               "f0235e31_5500000014_63X_20190816_S1_P2_B5_annotations_corrected.tiff", 17_104_143, 5, (1175, 1689), (231, 1107), 198_496),
    "cell_D": ("20190816_M02_001_Capture 3 - Capture 3 - B4 [144] Montage - 14",
               "32c432b4_5500000014_63X_20190816_S1_P14_B4_annotations_corrected.tiff", 17_590_873, 5, (1061, 1736), (245, 1143), 308_996),
    "cell_E": ("20190816_M02_001_Capture 3 - Capture 3 - B2 [144] Montage - 11",
               "b0344f39_5500000014_63X_20190816_S1_P11_B2_annotations_corrected.tiff", 16_728_186, 5, (165, 747), (63, 583), 110_252),
    "cell_F": ("20190816_M02_001_Capture 3 - Capture 3 - B5 [144] Montage - 9",
               "ccafc394_5500000014_63X_20190816_S1_P9_B5_annotations_corrected.tiff", 17_660_948, 5, (765, 1337), (0, 665), 144_864),
    "cell_G": ("20190816_M02_001_Capture 3 - Capture 3 - B5 [144] Montage - 15",
               "c8aa2331_5500000014_63X_20190816_S1_P15_B5_annotations_corrected.tiff", 17_936_091, 4, (315, 807), (707, 1583), 122_048),
    "cell_H": ("20190816_M02_001_Capture 3 - Capture 3 - B4 [144] Montage - 27",
               "e9264983_5500000014_63X_20190816_S1_P27_B4_annotations_corrected.tiff", 16_391_151, 4, (183, 915), (465, 1377), 305_148),
}

# file -> (file MD5, shape, dtype, pixel-bytes MD5, OME UUID)
OUTPUTS = {
    "cell_A_actn2.ome.tif": ("ad7537038e676a34318382ccd98a9b68", (775, 698), "uint16", "9ff36034a5b9ffbb2f01ecf3615662e3", "8c4dafc2-d6f3-11f1-b6b0-ae699ee1431c"),
    "cell_A_mask.ome.tif": ("c3a56d8b6d327e5168856097b180c3f3", (775, 698), "uint8", "6d551e9d9981fad4e8278a1315ff5a9b", "8c503a12-d6f3-11f1-b6b0-ae699ee1431c"),
    "cell_B_actn2.ome.tif": ("beb4cc0ec2ace9aba66f97cf9046ef97", (898, 822), "uint16", "5cef7482318c73ba88b05be0cf50c67d", "8d22de18-d6f3-11f1-b6b0-ae699ee1431c"),
    "cell_B_mask.ome.tif": ("4f52a53f277a9893def27d5ee989cdd8", (898, 822), "uint8", "565723a70a03f4d7bdbbe60497c0346d", "8d267898-d6f3-11f1-b6b0-ae699ee1431c"),
    "cell_C_actn2.ome.tif": ("055817ddffcead40078e650b22deb37e", (514, 876), "uint16", "3e77e0a0ce23423ba4973d4a8a1bb630", "8dfa4ae2-d6f3-11f1-b6b0-ae699ee1431c"),
    "cell_C_mask.ome.tif": ("35c500f6eba7bd0b1bf0fa047092b448", (514, 876), "uint8", "13d48eda1bbc9df19515336bf19aaf64", "8dfca2d8-d6f3-11f1-b6b0-ae699ee1431c"),
    "cell_D_actn2.ome.tif": ("5e4d2c8c3df3b65f3bee23f9563f9267", (675, 898), "uint16", "9b3cc1965868f6eddc0401d7cb9dc9a8", "90c68ede-d6f3-11f1-b6b0-ae699ee1431c"),
    "cell_D_mask.ome.tif": ("1606b497712cf7058422ad4b1f228114", (675, 898), "uint8", "b795c9114bbb4b5c013021c57582f58a", "90c9728e-d6f3-11f1-b6b0-ae699ee1431c"),
    "cell_E_actn2.ome.tif": ("d35107b8bbfe75b2a352ad80af80e68e", (582, 520), "uint16", "53922f2553860b395eb1c24c18053ac3", "9185f5c6-d6f3-11f1-b6b0-ae699ee1431c"),
    "cell_E_mask.ome.tif": ("0e3b59e1578d07d9d6a1790d992df24e", (582, 520), "uint8", "012cb67dc10200e8af6dbc923b2d2a3c", "9187eade-d6f3-11f1-b6b0-ae699ee1431c"),
    "cell_F_actn2.ome.tif": ("839179f7fed257a11d7347fc8a1b2b7d", (572, 665), "uint16", "3c86d4d8e4fd31bcd5236a43a992c0e4", "92638d00-d6f3-11f1-b6b0-ae699ee1431c"),
    "cell_F_mask.ome.tif": ("fb2ac4b1c0dd4dcc2d2d8590e7c8658d", (572, 665), "uint8", "954e0012122cf42dfb0d29a4334eb115", "92658074-d6f3-11f1-b6b0-ae699ee1431c"),
    "cell_G_actn2.ome.tif": ("13733937329313f3d7a564621ce60f93", (492, 876), "uint16", "0a0dd5b50b16353a0948eae0fa5c8c15", "93509672-d6f3-11f1-b6b0-ae699ee1431c"),
    "cell_G_mask.ome.tif": ("2153954e62f7d4305fab201df840f7d7", (492, 876), "uint8", "fca2e8b2fcec714b4b21a6062efb5646", "9352c082-d6f3-11f1-b6b0-ae699ee1431c"),
    "cell_H_actn2.ome.tif": ("a9eda9ff11a4cf026470ca5ebb24986d", (732, 912), "uint16", "096228d74c9061217e810557366c55f0", "96133b1c-d6f3-11f1-b6b0-ae699ee1431c"),
    "cell_H_mask.ome.tif": ("5f700b93b2e7d7a7dc3e1473fda7ba30", (732, 912), "uint8", "55140e86fd8762568a7fae5189ac2a84", "9616fcca-d6f3-11f1-b6b0-ae699ee1431c"),
}


def md5sum(path, chunk=1 << 22):
    h = hashlib.md5()
    with open(path, "rb") as f:
        for block in iter(lambda: f.read(chunk), b""):
            h.update(block)
    return h.hexdigest()


def download(key, size):
    dst = CACHE / Path(key).name
    if dst.exists() and dst.stat().st_size == size:
        print(f"  cached: {dst.name}")
        return dst
    CACHE.mkdir(parents=True, exist_ok=True)
    url = BUCKET + urllib.parse.quote(key)
    part = dst.with_name(dst.name + ".part")
    print(f"  downloading {url}")
    with requests.get(url, stream=True, timeout=60) as r:
        r.raise_for_status()
        t0, done = time.time(), 0
        with open(part, "wb") as f:
            for block in r.iter_content(1 << 20):
                f.write(block)
                done += len(block)
                rate = done / max(time.time() - t0, 1e-6) / 1e6
                sys.stderr.write(f"\r  {dst.name[:40]}: {done / 1e6:6.1f} / {size / 1e6:.1f} MB "
                                 f"({100 * done / size:5.1f}%) {rate:6.1f} MB/s")
                sys.stderr.flush()
    sys.stderr.write("\n")
    if part.stat().st_size != size:
        raise RuntimeError(f"{dst.name}: expected {size} bytes, got {part.stat().st_size}; please rerun")
    part.replace(dst)
    return dst


def verify(name):
    """Return 'byte-identical', 'pixel-identical' or None."""
    path = DATA / name
    if not path.exists():
        return None
    fmd5, shape, dtype, pmd5, _ = OUTPUTS[name]
    if md5sum(path) == fmd5:
        return "byte-identical"
    try:
        with tifffile.TiffFile(path) as t:
            a = t.asarray()
            ome = t.ome_metadata or ""
    except Exception:
        return None
    if a.shape != shape or str(a.dtype) != dtype:
        return None
    if hashlib.md5(np.ascontiguousarray(a).tobytes()).hexdigest() != pmd5:
        return None
    for k in ("PhysicalSizeX", "PhysicalSizeY"):
        m = re.search(rf'{k}="([0-9.]+)"', ome)
        if not m or abs(float(m.group(1)) - PX) > 1e-9:
            return None
    return "pixel-identical"


def write(name, arr, title):
    md = dict(axes="YX", UUID=OUTPUTS[name][4], Creator=CREATOR, Name=title,
              PhysicalSizeX=PX, PhysicalSizeXUnit="µm", PhysicalSizeY=PX, PhysicalSizeYUnit="µm")
    tmp = DATA / (name + ".part")
    tifffile.imwrite(tmp, arr, ome=True, compression="zlib", metadata=md)
    tmp.replace(DATA / name)


def build_cell(cid):
    fov, field_file, size, label, rows, cols, npx = PANEL[cid]
    a = tifffile.imread(download(FIELD_DIR + field_file, size))
    if a.ndim != 3 or a.shape[0] != 10:
        raise RuntimeError(f"{field_file}: unexpected shape {a.shape}")
    img, m = a[CH_IMAGE], a[CH_CELL] == label
    ys, xs = np.where(m)
    M = int(round(MARGIN_UM / PX))
    H, W = m.shape
    y0, y1 = max(ys.min() - M, 0), min(ys.max() + M + 1, H)
    x0, x1 = max(xs.min() - M, 0), min(xs.max() + M + 1, W)
    if (y0, y1) != rows or (x0, x1) != cols:
        raise RuntimeError(f"{cid}: crop {y0}:{y1}, {x0}:{x1} differs from the recorded {rows}, {cols}")
    crop, mk = img[y0:y1, x0:x1], m[y0:y1, x0:x1].astype(np.uint8)
    if int(mk.sum()) != npx:
        raise RuntimeError(f"{cid}: mask holds {int(mk.sum())} pixels, expected {npx}")
    write(f"{cid}_actn2.ome.tif", crop, f"{cid} ACTN2-mEGFP maximum projection")
    write(f"{cid}_mask.ome.tif", mk, f"{cid} mask (1 = target cell, 0 = other)")
    print(f"  {cid}: {fov}, label {label}, crop rows {y0}:{y1} cols {x0}:{x1}, mask {npx} px")


def main():
    t_start = time.time()
    DATA.mkdir(exist_ok=True)
    todo = [c for c in PANEL if not (verify(f"{c}_actn2.ome.tif") and verify(f"{c}_mask.ome.tif"))]
    if todo:
        print(f"[task B2] rebuilding {len(todo)} cell(s) from the Allen Institute bucket")
        for cid in todo:
            build_cell(cid)
    else:
        print("[task B2] all files present")
    failed = []
    for name in OUTPUTS:
        res = verify(name)
        print(f"  {'OK  ' if res else 'FAIL'} {name}  ({res or 'does not match'})")
        if not res:
            failed.append(name)
    if failed:
        sys.exit(f"[task B2] verification failed: {', '.join(failed)}")
    if "--clean-cache" in sys.argv:
        for _, field_file, *_ in PANEL.values():
            p = CACHE / field_file
            if p.exists():
                p.unlink()
                print(f"  removed {p.name}")
    print(f"[task B2] {len(OUTPUTS)} files verified in {time.time() - t_start:.1f} s")


if __name__ == "__main__":
    main()
