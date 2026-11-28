#!/usr/bin/env python3
"""Reduces the explosion example's outputs to what the tutorial shows: a few thinned traces,
the gather and three snapshots the solver drew, and the run's facts.

    python3 report/experiments/usage/explosion.py <the run's output directory> <data directory>
"""

import json
import shutil
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
import reportlib  # noqa: E402

run, data = Path(sys.argv[1]), Path(sys.argv[2])
facts = json.loads((run / "run.json").read_text())

# Traces at the epicentre and at three offsets along the line, thinned for plotting.
_, time = reportlib.load_npy(run / "time.npy")
positions = [
    float(line.split(",")[1]) for line in (run / "receivers.csv").read_text().splitlines()[1:]
]
source_x = facts["sources"][0]["position"][0]
offsets = [0.0, 120.0, 240.0, 360.0]
chosen = [positions.index(source_x + offset) for offset in offsets]
header = ["t"]
columns = [time]
for component in facts["components"]:
    (receivers, samples), values = reportlib.load_npy(run / f"{component}.npy")
    for offset, receiver in zip(offsets, chosen):
        header.append(f"{component}_{offset:.0f}")
        columns.append(values[receiver * samples : (receiver + 1) * samples])
keep = reportlib.thin(len(time), 300)
reportlib.write_csv(data / "traces.csv", header, ([c[i] for c in columns] for i in keep))

# The solver's own pictures: the gather of vz, and the section through the source at three
# times.
shutil.copyfile(run / "vz_gather.png", data / "vz_gather.png")
dt = facts["dt"]
snapshots = {}
for label, seconds in [("a", 0.15), ("b", 0.25), ("c", 0.35)]:
    stills = sorted((run / "snapshots").glob("vz_xz*_*.png"))
    best = min(stills, key=lambda p: abs(int(p.stem.rsplit("_", 1)[1]) * dt - seconds))
    shutil.copyfile(best, data / f"snapshot_{label}.png")
    snapshots[label] = int(best.stem.rsplit("_", 1)[1]) * dt

reportlib.write_values(
    data,
    {
        "nx": facts["points"][0],
        "ny": facts["points"][1],
        "nz": facts["points"][2],
        "spacing": facts["spacing"][0],
        "extent_x": facts["extent"][0],
        "extent_z": facts["extent"][2],
        "dt_ms": facts["dt"] * 1e3,
        "dt_stable_ms": facts["dt_stable"] * 1e3,
        "steps": facts["steps"],
        "duration": facts["duration"],
        "order": facts["order"],
        "points_per_wavelength": facts["points_per_wavelength"],
        "max_frequency": facts["max_frequency"],
        "memory_mib": facts["memory_bytes"] / 2**20,
        "receivers": facts["receivers"],
        "receiver_first_x": positions[0],
        "receiver_last_x": positions[-1],
        "source_x": source_x,
        "source_z": facts["sources"][0]["position"][2],
        "snapshot_a_time": snapshots["a"],
        "snapshot_b_time": snapshots["b"],
        "snapshot_c_time": snapshots["c"],
    },
)
