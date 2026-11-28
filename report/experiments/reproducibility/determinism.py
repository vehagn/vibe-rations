#!/usr/bin/env python3
"""Compares the trace files of the same run made three ways, byte for byte.

    python3 report/experiments/reproducibility/determinism.py <raw directory> <data directory>
"""

import hashlib
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
import reportlib  # noqa: E402

raw, data = Path(sys.argv[1]), Path(sys.argv[2])
runs = ["threads-all", "threads-one", "plain-loops"]
files = ["time.npy", "vx.npy", "vz.npy"]

rows = []
identical = True
for name in files:
    hashes = [hashlib.sha256((raw / run / name).read_bytes()).hexdigest() for run in runs]
    same = len(set(hashes)) == 1
    identical &= same
    rows.append([name, *(h[:16] for h in hashes), "yes" if same else "no"])
reportlib.write_csv(data / "hashes.csv", ["file", *runs, "identical"], rows)

kernels = [json.loads((raw / run / "run.json").read_text())["kernels"] for run in runs]
reportlib.write_values(
    data,
    {
        "identical": "identical" if identical else "different",
        "files": len(files),
        "kernels_default": kernels[0],
        "kernels_off": kernels[2],
    },
)
if not identical:
    sys.exit("the three runs' traces differ; see hashes.csv")
