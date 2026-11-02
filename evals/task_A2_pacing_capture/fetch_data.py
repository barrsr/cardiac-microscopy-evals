#!/usr/bin/env python3
"""Fetch the task A2 data snapshot (1:1 capture of 1 Hz field stimulation).

Source: Kobeissi H, et al. "Strain gauge platforms: Time-lapse microscopy dataset of
engineered cardiac microbundles." Dryad, doi:10.5061/dryad.sqv9s4nbg, mirrored as
Zenodo record 11120809. License: CC0 1.0.

The snapshot is all eleven TIFFs of Type1_condition_025.zip, taken unmodified, with the two
batch subfolders (condition_025/batch_000/ and condition_025/batch_001/) flattened into data/.

The archive (654,914,176 bytes) is downloaded once from Zenodo into the cache directory
(default ~/.cache/cardiac-microscopy-evals, override with the CME_CACHE environment
variable) and checked against its published MD5. Task A uses the same archive and the same
cache path, so running both scripts downloads it only once. Each TIFF is streamed out of the
archive and verified byte for byte against the MD5 table below.

The script is idempotent: files already present with the right MD5 are left alone, and the
archive is only downloaded if at least one file is missing or wrong.

usage: python fetch_data.py [--clean-cache]
    --clean-cache   delete the cached archive after a successful run
                    (leave it if you will also fetch task A)

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

ZIP_NAME = "Type1_condition_025.zip"
ZIP_URL = f"https://zenodo.org/records/11120809/files/{ZIP_NAME}?download=1"
ZIP_SIZE = 654_914_176
ZIP_MD5 = "012aca6b7861e1685fc3f38bc0585d20"
ZIP_PATH = CACHE / "zenodo_11120809" / ZIP_NAME

# output file -> (member path inside the archive, bytes, MD5)
FILES = {
    "Ex00628.tif": ("condition_025/batch_000/Ex00628.tif", 79_444_634, "a1d8dac26038f1c5c6798d3c5e32b709"),
    "Ex00629.tif": ("condition_025/batch_000/Ex00629.tif", 79_444_634, "6297d51e5222bfe3836fe26a2758b697"),
    "Ex00630.tif": ("condition_025/batch_000/Ex00630.tif", 79_444_934, "0f7f1525f1ce74bc21bed80c5eba4695"),
    "Ex00631.tif": ("condition_025/batch_000/Ex00631.tif", 79_444_934, "a9a15a39e2b52ab76d200e0b0bf66fd5"),
    "Ex00632.tif": ("condition_025/batch_000/Ex00632.tif", 79_444_634, "08da3e8f86812fcddf5a6122b3be7a5b"),
    "Ex00633.tif": ("condition_025/batch_000/Ex00633.tif", 79_444_934, "f82adafd815314c96b43f6f706597751"),
    "Ex00634.tif": ("condition_025/batch_000/Ex00634.tif", 79_444_934, "ef6d71faf5b59e3c5d5874d7ed0686e2"),
    "Ex00635.tif": ("condition_025/batch_001/Ex00635.tif", 79_446_434, "9be9371ba75065ddc888285048950393"),
    "Ex00636.tif": ("condition_025/batch_001/Ex00636.tif", 79_446_434, "a4ab40c931291f06577abb4c5f13074f"),
    "Ex00637.tif": ("condition_025/batch_001/Ex00637.tif", 79_446_134, "51f847c6244c3f8b66f592a2a7ed57a9"),
    "Ex00638.tif": ("condition_025/batch_001/Ex00638.tif", 79_446_134, "442af8713d7a78ff4a073f6cf2bfa233"),
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
        print(f"[task A2] fetching {len(todo)} file(s) from Zenodo 11120809")
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
        print("[task A2] all files present")

    failed = []
    for name, (_, size, md5) in FILES.items():
        ok = valid(DATA / name, size, md5)
        print(f"  {'OK  ' if ok else 'FAIL'} {name}  md5 {md5}  (byte-identical)")
        if not ok:
            failed.append(name)
    if failed:
        sys.exit(f"[task A2] verification failed: {', '.join(failed)}")
    if "--clean-cache" in sys.argv and ZIP_PATH.exists():
        ZIP_PATH.unlink()
        print(f"  removed {ZIP_PATH}")
    print(f"[task A2] {len(FILES)} file(s) verified in {time.time() - t_start:.1f} s")


if __name__ == "__main__":
    main()
