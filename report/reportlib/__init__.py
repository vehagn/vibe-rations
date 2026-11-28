"""What the report's scripts share: reading the solver's .npy files, writing the small CSV
and JSON files the report's LaTeX reads, and drawing a 2D field as a PNG.

Standard library only, so the report needs nothing beyond the dev container's Python. A
script under report/experiments/ finds this package with

    sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
    import reportlib
"""

from __future__ import annotations

import array
import ast
import json
import os
import struct
import sys
import zlib
from pathlib import Path
from typing import Iterable, Sequence

ROOT = Path(__file__).resolve().parents[2]


def data_dir(experiment: str) -> Path:
    """Where an experiment's committed data go: where `report.py run` says, and otherwise
    under out/report/, so that a script run by hand never touches the committed data."""
    given = os.environ.get("VIBE_RATIONS_REPORT_DATA")
    path = Path(given) if given else ROOT / "out" / "report" / "data" / experiment
    path.mkdir(parents=True, exist_ok=True)
    return path


def raw_dir(experiment: str) -> Path:
    """Where an experiment's bulky, uncommitted outputs go."""
    given = os.environ.get("VIBE_RATIONS_REPORT_RAW")
    path = Path(given) if given else ROOT / "out" / "report" / "raw" / experiment
    path.mkdir(parents=True, exist_ok=True)
    return path


def load_npy(path: Path | str) -> tuple[tuple[int, ...], array.array]:
    """An .npy file as the solver writes them (little-endian f4 or f8, C order): its shape
    and its values, row-major."""
    blob = Path(path).read_bytes()
    if blob[:6] != b"\x93NUMPY":
        raise ValueError(f"{path}: not an .npy file")
    major = blob[6]
    if major == 1:
        (length,) = struct.unpack("<H", blob[8:10])
        start = 10
    else:
        (length,) = struct.unpack("<I", blob[8:12])
        start = 12
    header = ast.literal_eval(blob[start : start + length].decode("latin-1"))
    if header["fortran_order"]:
        raise ValueError(f"{path}: Fortran order isn't supported")
    code = {"<f4": "f", "<f8": "d", "<i4": "i", "<i8": "q"}.get(header["descr"])
    if code is None:
        raise ValueError(f"{path}: unsupported dtype {header['descr']}")
    values = array.array(code)
    values.frombytes(blob[start + length :])
    if sys.byteorder != "little":
        values.byteswap()
    return tuple(header["shape"]), values


def thin(count: int, limit: int) -> range:
    """Indices of at most about `limit` of `count` samples, evenly spaced from the first:
    what a plotted curve needs, and what keeps the committed files small."""
    return range(0, count, max(1, -(-count // limit)))


def number(value: float, digits: int = 6) -> str:
    """A number with `digits` significant digits, as short as it gets."""
    return f"{value:.{digits}g}"


def write_csv(
    path: Path | str,
    header: Sequence[str],
    rows: Iterable[Sequence[object]],
    digits: int = 6,
) -> None:
    """A comma-separated table: a header and rows, floats to `digits` significant digits."""
    with Path(path).open("w") as f:
        f.write(",".join(header) + "\n")
        for row in rows:
            cells = [number(c, digits) if isinstance(c, float) else str(c) for c in row]
            f.write(",".join(cells) + "\n")


def write_values(directory: Path | str, values: dict) -> None:
    """values.json: the numbers and words the text quotes as \\val{<experiment>/<key>}."""
    path = Path(directory) / "values.json"
    path.write_text(json.dumps(values, indent=2, ensure_ascii=False) + "\n")


# The report's colour maps, as report/preamble.tex defines them for pgfplots (`vrdiverging`,
# `vrsequential`), so a PNG drawn here matches a colour bar drawn there.
COLORMAPS = {
    "diverging": [(0x1B, 0x64, 0xBE), (0xEF, 0xEF, 0xEF), (0xD6, 0x3A, 0x3B)],
    "sequential": [(0xF4, 0xF8, 0xFE), (0x9D, 0xC8, 0xF7), (0x25, 0x7D, 0xE0), (0x0E, 0x3D, 0x7E)],
}


def _colour(stops: list[tuple[int, int, int]], x: float) -> tuple[int, int, int]:
    x = min(max(x, 0.0), 1.0) * (len(stops) - 1)
    i = min(int(x), len(stops) - 2)
    t = x - i
    return tuple(round(a + (b - a) * t) for a, b in zip(stops[i], stops[i + 1]))


def write_png(
    path: Path | str,
    rows: Sequence[Sequence[float]],
    lo: float,
    hi: float,
    colormap: str = "sequential",
    scale: int = 1,
) -> None:
    """A 2D field as an RGB PNG, the first row on top, coloured from `lo` to `hi` with one
    of the report's colour maps; `scale` repeats each value over a square of pixels."""
    stops = COLORMAPS[colormap]
    span = (hi - lo) or 1.0
    height, width = len(rows), len(rows[0])
    lines = bytearray()
    for row in rows:
        line = bytearray([0])
        for value in row:
            line += bytes(_colour(stops, (value - lo) / span)) * scale
        lines += line * scale

    def chunk(kind: bytes, body: bytes) -> bytes:
        return (
            struct.pack(">I", len(body))
            + kind
            + body
            + struct.pack(">I", zlib.crc32(kind + body))
        )

    header = struct.pack(">IIBBBBB", width * scale, height * scale, 8, 2, 0, 0, 0)
    Path(path).write_bytes(
        b"\x89PNG\r\n\x1a\n"
        + chunk(b"IHDR", header)
        + chunk(b"IDAT", zlib.compress(bytes(lines), 9))
        + chunk(b"IEND", b"")
    )
