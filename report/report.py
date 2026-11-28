#!/usr/bin/env python3
"""The report's generator: runs the experiments behind its tables and figures, stamps each
with the commit and commands that made it, and builds the PDF.

    report/report.py list                     every experiment and whether its data are current
    report/report.py run usage/explosion      run experiments (names, --tier, --stale, --missing)
    report/report.py render                   write build/gen/*.tex from data/ and git
    report/report.py pdf                      render, then latexmk: build/vibe-rations-report.pdf
    report/report.py check                    fail unless every experiment's data are current

An experiment is a manifest, report/experiments/<section>/<name>.toml: a title, a tier and
the shell commands that make its data, run from the repository's root; its scripts and
inputs sit beside it as <name>.py and <name>/. Its results, reduced
to small text files, live in report/data/<section>/<name>/ and are committed, so the PDF
builds without rerunning anything; `run` clears that directory, runs the commands and writes
provenance.json beside the results. report/README.md has the conventions.

Python's standard library only (3.11 or later, for tomllib).
"""

from __future__ import annotations

import argparse
import datetime
import hashlib
import json
import os
import platform
import re
import shutil
import subprocess
import sys
import time
import tomllib
from pathlib import Path

REPORT = Path(__file__).resolve().parent
ROOT = REPORT.parent
EXPERIMENTS = REPORT / "experiments"
DATA = REPORT / "data"
BUILD = REPORT / "build"
GEN = BUILD / "gen"
RAW = BUILD / "raw"
PDF = BUILD / "vibe-rations-report.pdf"

TIERS = ("quick", "study", "heavy")
# What every experiment's results depend on, besides its own manifest and scripts and what
# its manifest lists under `depends`: the solver and how it is built.
ALWAYS = ("src", "Cargo.toml", "Cargo.lock", "rust-toolchain.toml")
# Committed data are meant to be small text; larger files are reported.
FILE_LIMIT = 512 * 1024
DIRECTORY_LIMIT = 2 * 1024 * 1024


def git(*args: str) -> str:
    """git's output, or an empty string where there is no repository or no such object."""
    try:
        done = subprocess.run(
            ["git", *args], cwd=ROOT, capture_output=True, text=True, check=False
        )
    except FileNotFoundError:
        return ""
    return done.stdout.strip() if done.returncode == 0 else ""


def object_hash(path: str) -> str:
    """The hash of a path's content in the commit checked out: a tree's or a file's."""
    return git("rev-parse", "--verify", "--quiet", f"HEAD:{path}")


class Experiment:
    def __init__(self, manifest: Path):
        self.manifest = manifest
        self.id = manifest.relative_to(EXPERIMENTS).with_suffix("").as_posix()
        spec = tomllib.loads(manifest.read_text())
        unknown = set(spec) - {"title", "tier", "depends", "command", "needs"}
        if unknown:
            raise SystemExit(f"{self.id}: unknown keys {sorted(unknown)}")
        self.title: str = spec["title"]
        self.tier: str = spec.get("tier", "quick")
        if self.tier not in TIERS:
            raise SystemExit(f"{self.id}: tier must be one of {TIERS}")
        self.needs: str = spec.get("needs", "")
        self.commands: list[dict] = spec.get("command", [])
        for command in self.commands:
            if set(command) - {"run", "output"} or "run" not in command:
                raise SystemExit(f"{self.id}: a command has `run` and optionally `output`")
        # The manifest and whatever sits beside it under the same name: its scripts
        # (<name>.py) and a directory of its inputs (<name>/, say the solver's
        # configurations).
        beside = sorted(
            p.relative_to(ROOT).as_posix()
            for p in [*manifest.parent.glob(manifest.stem + ".*"), manifest.with_suffix("")]
            if p.exists()
        )
        self.depends: list[str] = [*ALWAYS, *beside, *spec.get("depends", [])]
        self.data = DATA / self.id
        self.raw = RAW / self.id

    def shown(self, command: str) -> str:
        """A command as the report prints it: paths relative to the repository's root."""
        return command.format(
            data=self.data.relative_to(ROOT).as_posix(),
            raw=self.raw.relative_to(ROOT).as_posix(),
        )

    def provenance(self) -> dict | None:
        path = self.data / "provenance.json"
        return json.loads(path.read_text()) if path.exists() else None

    def status(self) -> tuple[str, list[str]]:
        """missing, dirty (made with uncommitted changes), stale (what it depends on has
        changed since) or current; and for stale, which paths changed."""
        made = self.provenance()
        if made is None:
            return "missing", []
        changed = [
            path
            for path in self.depends
            if made["sources"].get(path, "") != object_hash(path)
        ]
        if changed:
            return "stale", changed
        if made["dirty"]:
            return "dirty", made["dirty"]
        return "current", []


def experiments() -> list[Experiment]:
    # experiments/<section>/<name>.toml exactly: TOML files further down are inputs.
    return [Experiment(path) for path in sorted(EXPERIMENTS.glob("*/*.toml"))]


def machine() -> dict:
    """The machine a run was made on: its processor, threads, memory and system. No host
    name: this goes into a public repository. VIBE_RATIONS_REPORT_MACHINE gives it a label."""
    cpu = ""
    try:
        if sys.platform == "darwin":
            cpu = subprocess.run(
                ["sysctl", "-n", "machdep.cpu.brand_string"],
                capture_output=True,
                text=True,
                check=False,
            ).stdout.strip()
        else:
            for line in Path("/proc/cpuinfo").read_text().splitlines():
                if line.startswith("model name"):
                    cpu = line.split(":", 1)[1].strip()
                    break
    except OSError:
        pass
    memory = 0
    try:
        memory = os.sysconf("SC_PHYS_PAGES") * os.sysconf("SC_PAGE_SIZE")
    except (ValueError, OSError, AttributeError):
        pass
    return {
        "label": os.environ.get("VIBE_RATIONS_REPORT_MACHINE", ""),
        "cpu": cpu or platform.processor() or platform.machine(),
        "architecture": platform.machine(),
        "threads": os.cpu_count() or 0,
        "memory_gib": round(memory / 2**30, 1),
        "system": f"{platform.system()} {platform.release()}",
    }


def machine_text(m: dict) -> str:
    what = f"{m['cpu']}, {m['threads']} threads, {m['memory_gib']:g} GiB, {m['system']}"
    return f"{m['label']} ({what})" if m.get("label") else what


def build_solver() -> Path:
    """Builds the command-line program as a user would and returns its directory, which
    goes first on the commands' PATH so they read `vibe-rations run ...`."""
    print("building: cargo build --release --bin vibe-rations", flush=True)
    subprocess.run(
        ["cargo", "build", "--release", "--bin", "vibe-rations"], cwd=ROOT, check=True
    )
    target = os.environ.get("CARGO_TARGET_DIR")
    if not target:
        meta = subprocess.run(
            ["cargo", "metadata", "--format-version", "1", "--no-deps"],
            cwd=ROOT,
            capture_output=True,
            text=True,
            check=True,
        )
        target = json.loads(meta.stdout)["target_directory"]
    return Path(target) / "release"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as f:
        for block in iter(lambda: f.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def run_experiment(e: Experiment, solver: Path | None, toolchain: str) -> bool:
    print(f"\n== {e.id}: {e.title}", flush=True)
    dirty = git("status", "--porcelain", "--", *e.depends).splitlines()
    if dirty:
        print("   uncommitted changes in what it depends on; its data will say so:")
        for line in dirty:
            print(f"     {line}")
    # The directory holds exactly what this run makes; a run that fails leaves no
    # provenance, so its data count as missing (`git checkout` brings the old ones back).
    shutil.rmtree(e.data, ignore_errors=True)
    e.data.mkdir(parents=True)
    e.raw.mkdir(parents=True, exist_ok=True)
    env = dict(os.environ)
    env["VIBE_RATIONS_REPORT_DATA"] = str(e.data)
    env["VIBE_RATIONS_REPORT_RAW"] = str(e.raw)
    if solver:
        env["PATH"] = f"{solver}{os.pathsep}{env['PATH']}"
    ran = []
    for number, command in enumerate(e.commands, 1):
        shown = e.shown(command["run"])
        print(f"   $ {shown}", flush=True)
        log = e.raw / f"command-{number}.log"
        start = time.monotonic()
        # A command with `output` is one whose words the report shows: everything it
        # prints, to either stream, is kept in that file.
        keep = "output" in command
        with log.open("w") as out:
            done = subprocess.run(
                command["run"].format(data=e.data, raw=e.raw),
                shell=True,
                executable="/bin/bash",
                cwd=ROOT,
                env=env,
                stdout=subprocess.PIPE,
                stderr=subprocess.STDOUT if keep else out,
                text=True,
                check=False,
            )
            out.write("\n--- stdout ---\n" + done.stdout)
        seconds = time.monotonic() - start
        if done.returncode != 0:
            tail = "\n".join(log.read_text().splitlines()[-25:])
            print(f"   failed (exit {done.returncode}); the end of {log}:\n{tail}")
            return False
        entry = {"run": shown, "seconds": round(seconds, 1)}
        if keep:
            (e.data / command["output"]).write_text(done.stdout)
            entry["output"] = command["output"]
        ran.append(entry)
        print(f"     {seconds:.1f} s", flush=True)
    files = {}
    total = 0
    for path in sorted(p for p in e.data.rglob("*") if p.is_file()):
        size = path.stat().st_size
        total += size
        files[path.relative_to(e.data).as_posix()] = {"bytes": size, "sha256": sha256(path)}
        if size > FILE_LIMIT:
            print(f"   warning: {path.relative_to(ROOT)} is {size // 1024} KiB; thin it")
    if total > DIRECTORY_LIMIT:
        print(f"   warning: {e.data.relative_to(ROOT)} holds {total // 1024} KiB")
    provenance = {
        "experiment": e.id,
        "title": e.title,
        "tier": e.tier,
        "commit": git("rev-parse", "HEAD"),
        "commit_date": git("log", "-1", "--format=%cI"),
        "dirty": dirty,
        "sources": {path: object_hash(path) for path in e.depends},
        "date": datetime.datetime.now(datetime.timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"),
        "machine": machine(),
        "toolchain": toolchain,
        "commands": ran,
        "seconds": round(sum(c["seconds"] for c in ran), 1),
        "files": files,
    }
    (e.data / "provenance.json").write_text(json.dumps(provenance, indent=2) + "\n")
    return True


def cmd_run(args: argparse.Namespace) -> int:
    known = {e.id: e for e in experiments()}
    chosen: list[Experiment] = []
    for name in args.names:
        matches = [e for e in known.values() if e.id == name or e.id.startswith(name + "/")]
        if not matches:
            raise SystemExit(f"no experiment {name!r}; `report.py list` shows them")
        chosen += matches
    if not args.names:
        tiers = args.tier or ["quick"]
        chosen = [e for e in known.values() if e.tier in tiers]
    elif args.tier:
        chosen = [e for e in chosen if e.tier in args.tier]
    if args.stale or args.missing:
        wanted = {"missing"} | ({"stale", "dirty"} if args.stale else set())
        chosen = [e for e in chosen if e.status()[0] in wanted]
    if not chosen:
        print("nothing to run")
        return 0
    if args.dry_run:
        for e in chosen:
            print(f"{e.id} ({e.tier}, {e.status()[0]})")
            for command in e.commands:
                print(f"  $ {e.shown(command['run'])}")
        return 0
    solver = None if args.no_build else build_solver()
    toolchain = subprocess.run(
        ["rustc", "--version"], cwd=ROOT, capture_output=True, text=True, check=False
    ).stdout.strip()
    failed = []
    for e in chosen:
        if not run_experiment(e, solver, toolchain):
            failed.append(e.id)
            if not args.keep_going:
                break
    if failed:
        print(f"\nfailed: {', '.join(failed)}")
        return 1
    return 0


def cmd_list(_args: argparse.Namespace) -> int:
    rows = []
    for e in experiments():
        status, why = e.status()
        made = e.provenance() or {}
        rows.append(
            (
                e.id,
                e.tier,
                status,
                made.get("commit", "")[:7],
                made.get("date", "")[:10],
                f"{made['seconds']:.0f} s" if made else "",
                ", ".join(why) if status == "stale" else "",
            )
        )
    head = ("experiment", "tier", "status", "commit", "made", "took", "changed since")
    widths = [max(len(str(r[i])) for r in [head, *rows]) for i in range(len(head))]
    for row in [head, *rows]:
        print("  ".join(str(cell).ljust(w) for cell, w in zip(row, widths)).rstrip())
    return 0


def cmd_check(_args: argparse.Namespace) -> int:
    bad = [(e.id, e.status()[0]) for e in experiments() if e.status()[0] != "current"]
    for name, status in bad:
        print(f"{name}: {status}")
    return 1 if bad else 0


# ----------------------------------------------------------------------------------------
# What LaTeX reads: build/gen


def tex(text: object) -> str:
    """Plain text made safe for LaTeX."""
    replace = {
        "\\": r"\textbackslash{}",
        "&": r"\&",
        "%": r"\%",
        "$": r"\$",
        "#": r"\#",
        "_": r"\_",
        "{": r"\{",
        "}": r"\}",
        "~": r"\textasciitilde{}",
        "^": r"\textasciicircum{}",
    }
    return "".join(replace.get(ch, ch) for ch in str(text))


def value_text(value: object) -> str:
    """A value as \\val gives it: numbers as written (siunitx formats them), text escaped."""
    if isinstance(value, bool):
        return "true" if value else "false"
    if isinstance(value, (int, float)):
        return repr(value)
    return tex(value)


def duration(seconds: float) -> str:
    if seconds < 90:
        return f"{seconds:.0f} s"
    if seconds < 5400:
        return f"{seconds / 60:.0f} min"
    return f"{seconds / 3600:.1f} h"


def without_notes(bib: str) -> str:
    """The bibliography without its `note` fields: in docs/references.bib they record how
    each entry was checked, which is for the repository's maintainers, not for readers."""
    out = []
    i = 0
    while True:
        match = re.search(r"^[ \t]*note[ \t]*=[ \t]*\{", bib[i:], flags=re.MULTILINE)
        if not match:
            return "".join(out) + bib[i:]
        out.append(bib[i : i + match.start()])
        j = i + match.end()
        depth = 1
        while depth:
            depth += {"{": 1, "}": -1}.get(bib[j], 0)
            j += 1
        while bib[j] in ", \t":
            j += 1
        i = j + 1 if bib[j] == "\n" else j


def cmd_render(_args: argparse.Namespace) -> int:
    shutil.rmtree(GEN, ignore_errors=True)
    (GEN / "commands").mkdir(parents=True)
    bib = (ROOT / "docs" / "references.bib").read_text()
    (GEN / "references.bib").write_text(without_notes(bib))
    commit = git("rev-parse", "HEAD")
    changes = git("status", "--porcelain")
    build = [
        "% Written by report/report.py render: do not edit.",
        rf"\renewcommand{{\buildcommit}}{{{commit or 'unknown'}}}",
        rf"\renewcommand{{\buildcommitshort}}{{{commit[:7] or 'unknown'}}}",
        rf"\renewcommand{{\builddate}}{{{git('log', '-1', '--format=%cs') or 'unknown'}}}",
        rf"\builddirty{'true' if changes else 'false'}",
    ]
    (GEN / "build.tex").write_text("\n".join(build) + "\n")

    values = ["% Written by report/report.py render: do not edit."]
    defs = ["% Written by report/report.py render: do not edit."]
    appendix = ["% Written by report/report.py render: do not edit."]
    rows = []
    section = None
    for e in experiments():
        status, why = e.status()
        made = e.provenance() or {}
        short = made.get("commit", "")[:7]
        defs.append(
            rf"\experiment{{{e.id}}}{{{tex(e.title)}}}{{{e.tier}}}{{{status}}}"
            rf"{{{short}}}{{{made.get('date', '')[:10]}}}"
        )
        lines = [e.shown(c["run"]) for c in e.commands]
        name = e.id.replace("/", "--")
        (GEN / "commands" / f"{name}.sh").write_text("\n".join(lines) + "\n")
        path = e.data / "values.json"
        if path.exists():
            for key, value in json.loads(path.read_text()).items():
                values.append(rf"\setval{{{e.id}/{key}}}{{{value_text(value)}}}")

        if e.id.split("/")[0] != section:
            section = e.id.split("/")[0]
            appendix.append(rf"\subsection*{{Experiments of \code{{{section}}}}}")
        appendix.append(rf"\experimententry{{{e.id}}}{{{tex(e.title)}}}")
        facts = [f"Tier \\code{{{e.tier}}}"]
        if e.needs:
            facts.append(f"needs {tex(e.needs)}")
        if made:
            facts.append(
                f"made on {made['date'][:10]} at commit \\commit{{{made['commit']}}} "
                f"in {duration(made['seconds'])}, on {tex(machine_text(made['machine']))}, "
                f"with {tex(made['toolchain'])}"
            )
        if status == "missing":
            facts.append(r"\statusmark{missing}: no data yet")
        elif status == "stale":
            listed = ", ".join(rf"\code{{{p}}}" for p in why)
            facts.append(rf"\statusmark{{stale}}: {listed} changed since")
        elif status == "dirty":
            facts.append(r"\statusmark{dirty}: made with uncommitted changes")
        appendix.append("; ".join(facts) + ".")
        appendix.append(rf"\lstinputlisting[style=shell]{{\gen/commands/{name}.sh}}")
        if made and made["files"]:
            listed = ", ".join(rf"\code{{{f}}}" for f in made["files"])
            where = e.data.relative_to(ROOT).as_posix()
            appendix.append(rf"\noindent Files in \code{{{where}/}}: {listed}.")
        made_at = rf"\commit{{{made['commit']}}}" if made else "--"
        took = duration(made["seconds"]) if made else "--"
        rows.append(
            rf"\hyperref[exp:{e.id}]{{\code{{{e.id}}}}} & {e.tier} & {made_at} & "
            rf"{made.get('date', '--')[:10]} & {took} & \statusmark{{{status}}} \\"
        )
    (GEN / "values.tex").write_text("\n".join(values) + "\n")
    (GEN / "experiments.tex").write_text("\n".join(defs) + "\n")
    (GEN / "appendix.tex").write_text("\n\n".join(appendix) + "\n")
    (GEN / "table.tex").write_text("\n".join(rows) + "\n")
    print(f"wrote {GEN.relative_to(ROOT)}/")
    return 0


def cmd_pdf(args: argparse.Namespace) -> int:
    cmd_render(args)
    done = subprocess.run(["latexmk", "main.tex"], cwd=REPORT, check=False)
    built = BUILD / "latex" / "main.pdf"
    if done.returncode != 0 or not built.exists():
        print(f"latexmk failed; see {BUILD.relative_to(ROOT)}/latex/main.log")
        return 1
    shutil.copyfile(built, PDF)
    print(f"wrote {PDF.relative_to(ROOT)}")
    return 0


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    sub = parser.add_subparsers(dest="command", required=True)
    sub.add_parser("list", help="every experiment and whether its data are current")
    run = sub.add_parser("run", help="run experiments and record their provenance")
    run.add_argument("names", nargs="*", help="experiments, or sections of them")
    run.add_argument(
        "--tier",
        action="append",
        choices=TIERS,
        help="only these tiers (several allowed); with no names, the default is quick",
    )
    run.add_argument("--stale", action="store_true", help="only stale, dirty or missing ones")
    run.add_argument("--missing", action="store_true", help="only ones without data")
    run.add_argument("--keep-going", action="store_true", help="carry on after a failure")
    run.add_argument("--dry-run", action="store_true", help="print the commands only")
    run.add_argument(
        "--no-build", action="store_true", help="don't build the solver or change PATH"
    )
    sub.add_parser("render", help="write build/gen/*.tex from data/ and git")
    sub.add_parser("pdf", help="render, then build the PDF with latexmk")
    sub.add_parser("check", help="fail unless every experiment's data are current")
    args = parser.parse_args()
    return {
        "list": cmd_list,
        "run": cmd_run,
        "render": cmd_render,
        "pdf": cmd_pdf,
        "check": cmd_check,
    }[args.command](args)


if __name__ == "__main__":
    sys.exit(main())
