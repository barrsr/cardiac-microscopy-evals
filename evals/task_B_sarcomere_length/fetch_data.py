#!/usr/bin/env python3
"""Fetch and rebuild the task B data snapshot (sarcomere length from ACTN2-mEGFP z-lines).

Source: Gerbin KA, et al. "Cell states beyond transcriptomics: Integrating structural
organization and gene expression in hiPSC-derived cardiomyocytes." Cell Systems 12(6):670-687
(2021), doi:10.1016/j.cels.2021.05.001. Data from the Allen Institute for Cell Science public
bucket s3://allencell/aics/integrated_transcriptomics_structural_organization_hipsc_cm/,
read anonymously over plain HTTPS, so no AWS tools or credentials are needed.

License: Allen Institute for Cell Science terms of use
(https://www.allencell.org/terms-of-use.html), research and noncommercial use with
attribution. This repository does not redistribute the data; this script fetches it from the
Allen Institute directly.

What is rebuilt (plate 5500000014, well B4, field Montage-23, cell 3, day 32, fixed):
  data/actn2_zstack.ome.tif  rows 250:1250 and columns 600:1630 of the 1736 x 1776 raw
                             ACTN2-mEGFP stack (channel C5, 488 nm / 525 nm), all 50 z planes,
                             pixel values unchanged, uint16, ZYX.
  data/cell_mask.ome.tif     the same crop of the authors' 2D cell segmentation (channel 9,
                             "cell"), label 3 set to 1 and everything else 0, uint8.
                             It holds 501,844 pixels, the authors' Total_Area for this cell.
Both files carry PhysicalSizeX/Y = 0.1195257 um, the SlideBook calibration (tag 65007 and
TIFF XResolution) of the raw file; the stack also carries PhysicalSizeZ = 0.33 um. The
source OME-XML rounds X/Y to 0.120 um.

Intermediate downloads (about 325 MB) go to the cache directory, default
~/.cache/cardiac-microscopy-evals, overridable with the CME_CACHE environment variable.
Pass --clean-cache to delete them after a successful run.

Verification. The OME-TIFFs are written with zlib compression and with the OME UUID and
Creator string of the original snapshot fixed, so with tifffile 2026.10.15 and the same zlib
the output is byte-identical to the original and is checked against the file MD5s below.
Other tifffile or zlib versions can lay out or compress the same pixels into different bytes.
When the file MD5 does not match, the script falls back to checking the decoded pixel array
(shape, dtype and MD5 of the raw pixel bytes) plus the physical pixel size in the OME-XML, and
reports the file as "pixel-identical". Either result means the data the model sees is the same.

The script is idempotent: files that already verify are left alone.

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

PX_XY = 0.1195257  # um/px, SlideBook calibration of session 20190816_M02_001
PX_Z = 0.33        # um
ROWS, COLS = slice(250, 1250), slice(600, 1630)
LABEL = 3
MASK_PIXELS = 501_844

# key in the bucket -> (local cache name, bytes, MD5 of the object as downloaded on 2 Nov 2026)
INPUTS = {
    "raw_images/FISH/5500000014/20190816_M02_001_Capture 3 - Capture 3 - B4 [144] Montage - 23_C5.tif":
        ("B4_Montage-23_C5.tif", 308_456_737, "f5fb21a45863957476f7eb1a36b17b6f"),
    "2d_segmented_fields/2D_fov_tiff_path/3dc13480_5500000014_63X_20190816_S1_P23_B4_annotations_corrected.tiff":
        ("3dc13480_5500000014_63X_20190816_S1_P23_B4_annotations_corrected.tiff", 16_837_885,
         "869b9921b022729ce66289180bf20d8f"),
}
RAW_KEY, FIELD_KEY = list(INPUTS)

CREATOR = "tifffile.py 2026.10.15"
# file -> (file MD5, shape, dtype, pixel-bytes MD5, OME UUID, OME metadata)
OUTPUTS = {
    "actn2_zstack.ome.tif": (
        "4073aba0e6c5aa30ec4c660d0eab3566", (50, 1000, 1030), "uint16", "fffa65ee65b63e60db0d13be84226590",
        "c81fe86a-d6f0-11f1-84b1-ae699ee1431c",
        dict(axes="ZYX",
             Name="hiPSC-CM ACTN2-mEGFP, Allen Institute FISH plate 5500000014, day 32, field B4 Montage-23, cell 3 region",
             Channel={"Name": "ACTN2-mEGFP (488 nm ex, 525 nm em)"},
             PhysicalSizeZ=PX_Z, PhysicalSizeZUnit="µm")),
    "cell_mask.ome.tif": (
        "bf96b430a00d877c40b8f283a603bb07", (1000, 1030), "uint8", "d827a084f048caab907ef726b6c7e6c3",
        "c91f0368-d6f0-11f1-84b1-ae699ee1431c",
        dict(axes="YX", Name="target cell mask (1 = target cell, 0 = other)")),
}


def md5sum(path, chunk=1 << 22):
    h = hashlib.md5()
    with open(path, "rb") as f:
        for block in iter(lambda: f.read(chunk), b""):
            h.update(block)
    return h.hexdigest()


def download(key, name, size, md5):
    dst = CACHE / name
    if dst.exists() and dst.stat().st_size == size:
        print(f"  cached: {dst}")
        return dst
    CACHE.mkdir(parents=True, exist_ok=True)
    url = BUCKET + urllib.parse.quote(key)
    part = dst.with_name(dst.name + ".part")
    have = part.stat().st_size if part.exists() else 0
    headers = {"Range": f"bytes={have}-"} if 0 < have < size else {}
    print(f"  downloading {url}")
    with requests.get(url, stream=True, headers=headers, timeout=60) as r:
        r.raise_for_status()
        if have and r.status_code != 206:
            have = 0
        t0, done = time.time(), have
        with open(part, "ab" if have else "wb") as f:
            for block in r.iter_content(1 << 20):
                f.write(block)
                done += len(block)
                rate = (done - have) / max(time.time() - t0, 1e-6) / 1e6
                sys.stderr.write(f"\r  {name[:40]}: {done / 1e6:7.1f} / {size / 1e6:.1f} MB "
                                 f"({100 * done / size:5.1f}%) {rate:6.1f} MB/s")
                sys.stderr.flush()
    sys.stderr.write("\n")
    if part.stat().st_size != size:
        raise RuntimeError(f"{name}: expected {size} bytes, got {part.stat().st_size}; please rerun")
    got = md5sum(part)
    if got != md5:
        # The outputs are verified independently, so a changed upstream object is only a warning here.
        print(f"  warning: {name} MD5 {got} differs from the recorded {md5}")
    part.replace(dst)
    return dst


def verify(name):
    """Return 'byte-identical', 'pixel-identical' or None."""
    path = DATA / name
    if not path.exists():
        return None
    fmd5, shape, dtype, pmd5, _, meta = OUTPUTS[name]
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
    sizes = {k: re.search(rf'{k}="([0-9.]+)"', ome) for k in ("PhysicalSizeX", "PhysicalSizeY", "PhysicalSizeZ")}
    if not all(sizes[k] and abs(float(sizes[k].group(1)) - PX_XY) < 1e-9 for k in ("PhysicalSizeX", "PhysicalSizeY")):
        return None
    if "PhysicalSizeZ" in meta and not (sizes["PhysicalSizeZ"] and abs(float(sizes["PhysicalSizeZ"].group(1)) - PX_Z) < 1e-9):
        return None
    return "pixel-identical"


def write(name, arr):
    _, _, _, _, uuid, meta = OUTPUTS[name]
    md = dict(meta, UUID=uuid, Creator=CREATOR, PhysicalSizeX=PX_XY, PhysicalSizeXUnit="µm",
              PhysicalSizeY=PX_XY, PhysicalSizeYUnit="µm")
    tmp = DATA / (name + ".part")
    tifffile.imwrite(tmp, arr, ome=True, compression="zlib", metadata=md)
    tmp.replace(DATA / name)


def build():
    raw = download(RAW_KEY, *INPUTS[RAW_KEY])
    field = download(FIELD_KEY, *INPUTS[FIELD_KEY])
    print("  cropping the raw ACTN2-mEGFP stack")
    # Read the 50 TIFF pages directly; the embedded OME-XML points at sibling channel files.
    with tifffile.TiffFile(raw, is_ome=False) as t:
        stack = np.stack([p.asarray()[ROWS, COLS] for p in t.pages])
    if stack.shape != (50, 1000, 1030) or stack.dtype != np.uint16:
        raise RuntimeError(f"unexpected raw stack {stack.shape} {stack.dtype}")
    seg = tifffile.imread(field)
    if seg.shape[0] != 10:
        raise RuntimeError(f"unexpected 2D field shape {seg.shape}")
    mask = (seg[9][ROWS, COLS] == LABEL).astype(np.uint8)  # channel_defs.json: 'cell' = 9
    if int(mask.sum()) != MASK_PIXELS:
        raise RuntimeError(f"mask holds {int(mask.sum())} pixels, expected {MASK_PIXELS}")
    write("actn2_zstack.ome.tif", stack)
    write("cell_mask.ome.tif", mask)


def main():
    t_start = time.time()
    DATA.mkdir(exist_ok=True)
    if all(verify(n) for n in OUTPUTS):
        print("[task B] all files present")
    else:
        print("[task B] rebuilding from the Allen Institute bucket")
        build()
    failed = []
    for name in OUTPUTS:
        res = verify(name)
        print(f"  {'OK  ' if res else 'FAIL'} {name}  ({res or 'does not match'})")
        if not res:
            failed.append(name)
    if failed:
        sys.exit(f"[task B] verification failed: {', '.join(failed)}")
    if "--clean-cache" in sys.argv:
        for name, *_ in INPUTS.values():
            p = CACHE / name
            if p.exists():
                p.unlink()
                print(f"  removed {p}")
    print(f"[task B] {len(OUTPUTS)} files verified in {time.time() - t_start:.1f} s")


if __name__ == "__main__":
    main()
