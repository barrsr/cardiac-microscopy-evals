#!/usr/bin/env python3
"""Fetch and verify the data for all four evals.

Runs evals/<task>/fetch_data.py for each task in turn with the same Python interpreter.
Tasks A and A2 share one Zenodo archive in the cache, so it is downloaded once. Each script is
idempotent, so rerunning this only fills in whatever is missing or does not verify.

Cache directory: ~/.cache/cardiac-microscopy-evals, or the CME_CACHE environment variable.

usage: python scripts/fetch_all.py [--clean-cache] [task ...]
    task           one or more task folder names (default: all four)
    --clean-cache  after every task has verified, delete the cached downloads
"""
import subprocess
import sys
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
TASKS = [
    "task_A_beat_rate",
    "task_A2_pacing_capture",
    "task_B_sarcomere_length",
    "task_B2_sarcomere_maturity",
]


def main():
    args = sys.argv[1:]
    clean = "--clean-cache" in args
    chosen = [a for a in args if not a.startswith("--")] or TASKS
    unknown = [t for t in chosen if t not in TASKS]
    if unknown:
        sys.exit(f"unknown task(s): {', '.join(unknown)}; choose from {', '.join(TASKS)}")
    results = {}
    for i, task in enumerate(chosen):
        script = ROOT / "evals" / task / "fetch_data.py"
        # Clean the shared Zenodo archive only after the last task that needs it.
        later = chosen[i + 1:]
        cmd = [sys.executable, str(script)]
        if clean and not (task.startswith("task_A") and any(t.startswith("task_A") for t in later)):
            cmd.append("--clean-cache")
        t0 = time.time()
        rc = subprocess.call(cmd, cwd=script.parent)
        results[task] = (rc, time.time() - t0)
        print()
    print("summary")
    for task, (rc, dt) in results.items():
        print(f"  {'OK  ' if rc == 0 else 'FAIL'} {task}  ({dt:.1f} s)")
    sys.exit(0 if all(rc == 0 for rc, _ in results.values()) else 1)


if __name__ == "__main__":
    main()
