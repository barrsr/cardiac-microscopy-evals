#!/usr/bin/env python3
"""Fetch the task A3 data snapshot (1:1 capture versus escape rhythm under 1 Hz pacing).

Source: Kobeissi H, et al. "Strain gauge platforms: Time-lapse microscopy dataset of
engineered cardiac microbundles." Dryad, doi:10.5061/dryad.sqv9s4nbg, mirrored as
Zenodo record 11120809. License: CC0 1.0.

The snapshot is nineteen of the twenty-three TIFFs of Type1_condition_023.zip, taken
unmodified, with the three batch subfolders (condition_023/batch_000/ to batch_002/) flattened
into data/. Ex00568, Ex00573, Ex00577 and Ex00578 are left out; eval.json gives the reasons.

The archive (1,258,653,935 bytes) is downloaded once from Zenodo into the cache directory
(default ~/.cache/cardiac-microscopy-evals, override with the CME_CACHE environment
variable) and checked against its published MD5. Each TIFF is streamed out of the archive and
verified byte for byte against the MD5 table below.

The script is idempotent: files already present with the right MD5 are left alone, and the
archive is only downloaded if at least one file is missing or wrong.

usage: python fetch_data.py [--clean-cache]
    --clean-cache   delete the cached archive after a successful run

Requires: requests (plus the standard library).
"""
import hashlib
import os
import shutil
import sys
import time
import zipfile
from pathlib import Path

import requests

HERE = Path(__file__).resolve().parent
DATA = HERE / "data"
CACHE = Path(os.environ.get("CME_CACHE", "~/.cache/cardiac-microscopy-evals")).expanduser()

ZIP_NAME = "Type1_condition_023.zip"
ZIP_URL = f"https://zenodo.org/records/11120809/files/{ZIP_NAME}?download=1"
ZIP_SIZE = 1_258_653_935
ZIP_MD5 = "1dc07af0f95876443ef555d8968a6e12"
ZIP_PATH = CACHE / "zenodo_11120809" / ZIP_NAME

# output file -> (member path inside the archive, bytes, MD5)
FILES = {
    "Ex00556.tif": ("condition_023/batch_000/Ex00556.tif", 79_444_634, "f852aba118a6e3635ed722c78e267aba"),
    "Ex00557.tif": ("condition_023/batch_000/Ex00557.tif", 79_444_334, "eba14f60b2e1199fe329cf90cb00e9bd"),
    "Ex00558.tif": ("condition_023/batch_000/Ex00558.tif", 79_444_634, "1534e95e0591974b8ffa9db9aaf70bb5"),
    "Ex00559.tif": ("condition_023/batch_000/Ex00559.tif", 79_444_634, "39c5dccfcfd0024e502acf234a23a088"),
    "Ex00560.tif": ("condition_023/batch_001/Ex00560.tif", 79_446_134, "60d56436f70cfbc3e25652cdaa535b45"),
    "Ex00561.tif": ("condition_023/batch_001/Ex00561.tif", 79_446_434, "7a1b7af2f2c3ab9a2d716d63dbc5e08b"),
    "Ex00562.tif": ("condition_023/batch_001/Ex00562.tif", 79_446_134, "2c63c39958af34fe9dd233d656609c3c"),
    "Ex00563.tif": ("condition_023/batch_001/Ex00563.tif", 79_446_134, "ecc30082b78797a7195438feec775d5c"),
    "Ex00564.tif": ("condition_023/batch_001/Ex00564.tif", 79_446_134, "f9114e1401a774da5ca889a5a422dc8d"),
    "Ex00565.tif": ("condition_023/batch_001/Ex00565.tif", 79_446_134, "b21ee8c8f7c7b2e6d9c0e130ee08084e"),
    "Ex00566.tif": ("condition_023/batch_002/Ex00566.tif", 79_444_634, "54f191740d22346676f9b4faf6051b52"),
    "Ex00567.tif": ("condition_023/batch_002/Ex00567.tif", 79_445_234, "7a6231ed2e58f531856e737be24af90c"),
    "Ex00569.tif": ("condition_023/batch_002/Ex00569.tif", 79_445_234, "f33becfea170485932f670972963374d"),
    "Ex00570.tif": ("condition_023/batch_002/Ex00570.tif", 79_444_634, "a8322af5513e517a39b240d7b99e602e"),
    "Ex00571.tif": ("condition_023/batch_002/Ex00571.tif", 79_444_634, "58978aa85ada6ff5f478a6dff1dbf74b"),
    "Ex00572.tif": ("condition_023/batch_002/Ex00572.tif", 79_444_934, "9bd82c9ebaff21e54f7241e00254e505"),
    "Ex00574.tif": ("condition_023/batch_002/Ex00574.tif", 79_445_234, "0a680ffc9ef8a5e8dead40efccbd4c35"),
    "Ex00575.tif": ("condition_023/batch_002/Ex00575.tif", 79_444_634, "ed71d33e8c0e1b95e3f270d0a11611fe"),
    "Ex00576.tif": ("condition_023/batch_002/Ex00576.tif", 79_444_634, "46c5bc14dcae12305459d5721e174591"),
}


def md5sum(path, chunk=1 << 22):
    h = hashlib.md5()
    with open(path, "rb") as f:
        for block in iter(lambda: f.read(chunk), b""):
            h.update(block)
    return h.hexdigest()


def progress(label, done, total, t0):
    rate = done / max(time.time() - t0, 1e-6) / 1e6
    pct = 100.0 * done / total if total else 0.0
    sys.stderr.write(f"\r  {label}: {done / 1e6:8.1f} / {total / 1e6:.1f} MB ({pct:5.1f}%) {rate:6.1f} MB/s")
    sys.stderr.flush()


def download(url, dst, size, md5):
    """Stream url to dst, resuming a partial download, then check size and MD5."""
    if dst.exists() and dst.stat().st_size == size and md5sum(dst) == md5:
        print(f"  cached: {dst}")
        return
    dst.parent.mkdir(parents=True, exist_ok=True)
    part = dst.with_name(dst.name + ".part")
    have = part.stat().st_size if part.exists() else 0
    headers = {"Range": f"bytes={have}-"} if 0 < have < size else {}
    if have >= size:
        part.unlink()
        have = 0
    print(f"  downloading {url}")
    with requests.get(url, stream=True, headers=headers, timeout=60) as r:
        r.raise_for_status()
        if have and r.status_code != 206:
            have = 0  # server ignored the range request; start over
        t0, done = time.time(), have
        with open(part, "ab" if have else "wb") as f:
            for block in r.iter_content(1 << 20):
                f.write(block)
                done += len(block)
                progress(dst.name, done, size, t0)
    sys.stderr.write("\n")
    if part.stat().st_size != size or md5sum(part) != md5:
        part.unlink()
        raise RuntimeError(f"{dst.name}: size or MD5 mismatch after download; please rerun")
    part.replace(dst)


def valid(path, size, md5):
    return path.exists() and path.stat().st_size == size and md5sum(path) == md5


def main():
    t_start = time.time()
    DATA.mkdir(exist_ok=True)
    todo = {k: v for k, v in FILES.items() if not valid(DATA / k, v[1], v[2])}
    if todo:
        print(f"[task A3] fetching {len(todo)} file(s) from Zenodo 11120809")
        download(ZIP_URL, ZIP_PATH, ZIP_SIZE, ZIP_MD5)
        with zipfile.ZipFile(ZIP_PATH) as z:
            for name, (member, size, md5) in todo.items():
                out = DATA / name
                tmp = out.with_name(out.name + ".part")
                with z.open(member) as src, open(tmp, "wb") as dst:
                    shutil.copyfileobj(src, dst, 1 << 22)
                tmp.replace(out)
                print(f"  extracted {member} -> data/{name}")
    else:
        print("[task A3] all files present")

    failed = []
    for name, (_, size, md5) in FILES.items():
        ok = valid(DATA / name, size, md5)
        print(f"  {'OK  ' if ok else 'FAIL'} {name}  md5 {md5}  (byte-identical)")
        if not ok:
            failed.append(name)
    if failed:
        sys.exit(f"[task A3] verification failed: {', '.join(failed)}")
    if "--clean-cache" in sys.argv and ZIP_PATH.exists():
        ZIP_PATH.unlink()
        print(f"  removed {ZIP_PATH}")
    print(f"[task A3] {len(FILES)} file(s) verified in {time.time() - t_start:.1f} s")


if __name__ == "__main__":
    main()
