#!/usr/bin/env python3
"""Regenerate the deposit-facing DNLS paper figure assets."""

from __future__ import annotations

import subprocess
import sys
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[2]
SCRIPTS = [
    REPO_ROOT / "d2_natural_lengths.py",
    REPO_ROOT / "fss_analyze.py",
]


def main() -> int:
    for script in SCRIPTS:
        print(f"Running {script.name} ...")
        subprocess.run([sys.executable, str(script)], cwd=REPO_ROOT, check=True)
    print("Done.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
