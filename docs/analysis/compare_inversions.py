#!/usr/bin/env python3
"""Compares inversions of one configuration run on the CPU and on a GPU.

Usage: compare_inversions.py NAME CONFIG TRUTH_DIR RUN_DIR... (the first run is the reference)

Each RUN_DIR is what docs/analysis/run_timed.sh leaves: `log` (every line timestamped, UTC),
`iterations/<n>/{vp,vs}.npy` and `inversion.json` (the model and summary after iteration n;
0 is the start) and `out/inversion.json`. TRUTH_DIR holds the true vp.npy and vs.npy,
sections (nz, nx) or grids (nz, ny, nx). Writes NAME.md and NAME.json beside the runs.
"""
import datetime
import json
import os
import re
import sys
import tomllib

import numpy as np

name, config, truth_dir, *runs = sys.argv[1:]
cfg = tomllib.load(open(config, "rb"))
nx, ny, nz = cfg["grid"]["points"]
h = cfg["grid"]["spacing"]
hx, hy, hz = (h, h, h) if not isinstance(h, list) else h
pml = cfg.get("boundary", {}).get("pml_width", 10)
fixed = cfg["inversion"].get("fixed_above", 0.0)
lines = cfg["receivers"].get("line", [])
y = lines[0]["from"][1] if lines else cfg["receivers"]["point"][0]["position"][1]
jr = int(round(y / hy))
ACROSS = 4  # nodes either side of the receivers' plane, as tests/marmousi2_study.rs judges

k0 = int(np.ceil(fixed / hz - 1e-9))
region = (slice(k0, nz - pml), slice(max(jr - ACROSS, 0), min(jr + ACROSS + 1, ny)), slice(pml, nx - pml))
if ny == 1:
    region = (region[0], slice(0, 1), region[2])
bands = np.array_split(np.arange(k0, nz - pml), 4)


def grid(a):
    a = np.asarray(a, dtype=np.float64)
    return a.reshape(nz, 1, nx).repeat(ny, axis=1) if a.ndim == 2 else a.reshape(nz, ny, nx)


truth = {p: grid(np.load(os.path.join(truth_dir, p + ".npy"))) for p in ("vp", "vs")}


def stamp(line):
    return datetime.datetime.strptime(line[:23], "%Y-%m-%dT%H:%M:%S.%f")


def timing(run):
    t = {"iterations": []}
    with open(os.path.join(run, "log")) as f:
        for line in f:
            text = line[25:].strip()
            if text.startswith("start:"):
                t["start"] = stamp(line)
            elif "unknowns;" in text:
                t["ready"] = stamp(line)
            elif re.match(r"iteration +\d+:", text):
                t["iterations"].append(stamp(line))
            elif text.startswith("exit"):
                t["end"], t["exit"] = stamp(line), text.split()[1]
            elif "Maximum resident set size" in text:
                t["rss_gib"] = int(text.split()[-1]) / 2**20
    out = {"exit": t.get("exit"), "rss_gib": t.get("rss_gib")}
    if "ready" in t:
        out["start_s"] = (t["ready"] - t["start"]).total_seconds()
        marks = [t["ready"]] + t["iterations"]
        out["iteration_s"] = [(b - a).total_seconds() for a, b in zip(marks, marks[1:])]
    if "end" in t:
        out["total_s"] = (t["end"] - t["start"]).total_seconds()
    return out


def iterations(run):
    d = os.path.join(run, "iterations")
    return sorted(int(n) for n in os.listdir(d) if n.isdigit()) if os.path.isdir(d) else []


def model(run, n, p):
    return np.load(os.path.join(run, "iterations", str(n), p + ".npy"))


def summary(run):
    path = os.path.join(run, "out", "inversion.json")
    return json.load(open(path)) if os.path.exists(path) else {}


def errors(run, n):
    """For vp and vs: |m - true| over |start - true| in the region, whole and by depth band."""
    out = {}
    for p in ("vp", "vs"):
        m, m0, t = grid(model(run, n, p)), grid(model(run, 0, p)), truth[p]
        e, e0 = (m - t)[region], (m0 - t)[region]
        whole = np.linalg.norm(e) / np.linalg.norm(e0)
        by_band = []
        for b in bands:
            rows = slice(b[0] - k0, b[-1] - k0 + 1)
            by_band.append(float(np.linalg.norm(e[rows]) / np.linalg.norm(e0[rows])))
        out[p] = {"whole": float(whole), "bands": by_band}
    return out


result = {"name": name, "config": config, "runs": {}}
for run in runs:
    s = summary(run)
    its = iterations(run)
    misfits = [s.get("starting_misfit")] + [i["misfit"] for i in s.get("iterations", [])]
    r = {
        "device": s.get("device") or "CPU",
        "timing": timing(run),
        "misfits": misfits,
        "trials": [i["trials"] for i in s.get("iterations", [])],
        "memory_bytes": s.get("memory_bytes"),
        "stopped": s.get("stopped"),
        "error": {n: errors(run, n) for n in its} if its else {},
    }
    result["runs"][os.path.basename(run.rstrip("/"))] = r

ref = runs[0]
for run in runs[1:]:
    diff = []
    for n in sorted(set(iterations(ref)) & set(iterations(run))):
        row = {"iteration": n}
        for p in ("vp", "vs"):
            a, b = model(ref, n, p), model(run, n, p)
            differ = int(np.count_nonzero(a.view(np.uint32) != b.view(np.uint32)))
            rel = np.abs(b.astype(np.float64) - a) / np.maximum(np.abs(a), 1e-30)
            row[p] = {"differ": differ, "of": int(a.size), "max_rel": float(rel.max()),
                      "rms_rel": float(np.sqrt((rel**2).mean()))}
        diff.append(row)
    result["runs"][os.path.basename(run.rstrip("/"))]["against_" + os.path.basename(ref.rstrip("/"))] = diff

here = os.path.dirname(runs[0].rstrip("/"))
json.dump(result, open(os.path.join(here, name + ".json"), "w"), indent=1, default=str)

# A page to read.
names = list(result["runs"])
md = [f"# {name}: CPU against GPU\n", f"Configuration `{config}`; error judged at depths "
      f"{k0 * hz:.0f}–{(nz - pml - 1) * hz:.0f} m, x outside the layers, within {ACROSS} nodes of the receivers' plane.\n"]
md.append("| run | device | start (s) | mean iteration (s) | iterations | total (h) | max RSS (GiB) | exit |")
md.append("|---|---|---|---|---|---|---|---|")
for n in names:
    r = result["runs"][n]; t = r["timing"]; it = t.get("iteration_s", [])
    md.append(f"| {n} | {r['device'][:60]} | {t.get('start_s', float('nan')):.0f} | "
              f"{(sum(it) / len(it)) if it else float('nan'):.0f} | {len(it)} | "
              f"{t.get('total_s', float('nan')) / 3600:.2f} | {t.get('rss_gib') or float('nan'):.2f} | {t.get('exit')} |")
md.append("\n## Misfit by iteration (fraction of the start)\n")
md.append("| iteration | " + " | ".join(names) + " |")
md.append("|---" * (len(names) + 1) + "|")
longest = max(len(result["runs"][n]["misfits"]) for n in names)
for i in range(longest):
    cells = []
    for n in names:
        m = result["runs"][n]["misfits"]
        cells.append(f"{m[i] / m[0]:.6f}" if i < len(m) and m[0] else "")
    md.append(f"| {i} | " + " | ".join(cells) + " |")
for n in names[1:]:
    d = result["runs"][n].get("against_" + names[0], [])
    same = [row["iteration"] for row in d if row["vp"]["differ"] == 0 and row["vs"]["differ"] == 0]
    md.append(f"\n## {n} against {names[0]}\n")
    md.append(f"Models identical bit for bit at {len(same)} of {len(d)} iterations compared.\n")
    md.append("| iteration | vp values differing | vp max rel | vs values differing | vs max rel |")
    md.append("|---|---|---|---|---|")
    for row in d:
        md.append(f"| {row['iteration']} | {row['vp']['differ']} | {row['vp']['max_rel']:.2e} | "
                  f"{row['vs']['differ']} | {row['vs']['max_rel']:.2e} |")
md.append("\n## Error against the truth, final iteration (|m - true| / |start - true|)\n")
labels = [f"{b[0] * hz:.0f}–{b[-1] * hz:.0f} m" for b in bands]
md.append("| run | parameter | whole | " + " | ".join(labels) + " |")
md.append("|---" * (len(labels) + 3) + "|")
for n in names:
    e = result["runs"][n]["error"]
    if not e:
        continue
    last = e[max(e)]
    for p in ("vp", "vs"):
        md.append(f"| {n} (it {max(e)}) | {p} | {last[p]['whole']:.4f} | " +
                  " | ".join(f"{x:.4f}" for x in last[p]["bands"]) + " |")
open(os.path.join(here, name + ".md"), "w").write("\n".join(md) + "\n")
print("\n".join(md))
