"""What removing leapfrog's time error would buy: von Neumann analysis of the time-dispersion
transform (TDT), and of staggered-to-node interpolation.

Leapfrog turns the semi-discrete (time-continuous) frequency w' into the discrete one by
w' dt / 2 = sin(w dt / 2), whatever the medium and the spatial operator. A transform pair
(pre-warping the source wavelet, un-warping the traces; Wang & Xu 2015, Koene et al. 2018)
undoes that, leaving only the spatial error, for which

    w' dt = 2 r |K|,   K_i = sum_n c_n sin((2n - 1) k_i h / 2),   r = v dt / h.

Prints, for S waves at vp/vs = sqrt(3) and the worst of axis and diagonal:
  1. the points per S wavelength (ppw) each scheme needs for a group-velocity error of 0.5%
     and 0.1%, and its work relative to today's default (order 8 at C = 0.3), with work per
     unit volume and time ~ ppw^4 sum|c_n| / (C * speed), speeds measured in Mcell/s;
  2. the same for least-squares operators, judged without time error;
  3. the amplitude that midpoint interpolation from staggered positions to a node keeps.
Run: python3 docs/analysis/time_dispersion.py (pure Python, about a minute).
Results: docs/research-log.md, 2026-11-20, "A research brief, checked".
"""

import math

from dispersion import AXIS, DIAGONAL, K, least_squares, taylor

VPVS = math.sqrt(3)
# Measured Mcell/s on the container's four cores, by half-length (orders 4, 8, 12, 16).
SPEED = {2: 220.0, 4: 160.0, 6: 130.0, 8: 105.0}


def wdt(c, kh, direction, r, transform):
    s = math.sqrt(sum((K(c, kh * d) / 2) ** 2 for d in direction))
    x = r * s
    if transform:
        return 2 * x  # the semi-discrete relation, which the transform recovers
    return 2 * math.asin(x) if x <= 1 else None


def group_error(c, kh, direction, r, transform):
    e = 1e-6
    hi, lo = wdt(c, kh + e, direction, r, transform), wdt(c, kh - e, direction, r, transform)
    if hi is None or lo is None:
        return math.inf
    return (hi - lo) / (2 * e) / r - 1


def needed_ppw(c, r, tol, transform):
    for ppw10 in range(20, 400):
        ppw = ppw10 / 10
        khmax = 2 * math.pi / ppw
        if all(
            abs(group_error(c, khmax * s / 60, d, r, transform)) <= tol
            for s in range(1, 61)
            for d in (AXIS, DIAGONAL)
        ):
            return ppw
    return math.inf


def r_s(c, courant):
    """S-wave v dt / h at a fraction `courant` of the P-wave stability limit."""
    return courant / (math.sqrt(3) * sum(abs(x) for x in c)) / VPVS


def work(ppw, c, courant, L):
    return ppw**4 * sum(abs(x) for x in c) / (courant * SPEED[L])


def midpoint_weights(L):
    """Lagrange weights at 0 from the 2L staggered positions +-(n - 1/2), n = 1..L."""
    points = [n - 0.5 for n in range(1, L + 1)]
    every = points + [-x for x in points]
    out = []
    for x in points:
        a = 1.0
        for y in every:
            if y != x:
                a *= -y / (x - y)
        out.append(a)
    return out


if __name__ == "__main__":
    tols = (0.005, 0.001)
    c8 = taylor(4)
    ref = {t: needed_ppw(c8, r_s(c8, 0.3), t, False) for t in tols}
    ref_work = {t: work(ref[t], c8, 0.3, 4) for t in tols}

    print("1. S group-velocity error: ppw needed, and work relative to order 8 at C = 0.3")
    print(f"   {'scheme':28} {'ppw 0.5%':>9} {'work':>6} {'ppw 0.1%':>9} {'work':>6}")
    for L in (2, 4, 6, 8):
        c = taylor(L)
        for courant, transform in ((0.85, False), (0.3, False), (0.85, True), (0.95, True)):
            label = f"order {2 * L}, C = {courant}" + (" + TDT" if transform else "")
            cells = []
            for t in tols:
                p = needed_ppw(c, r_s(c, courant), t, transform)
                cells.append(f"{p:9.1f} {work(p, c, courant, L) / ref_work[t]:6.2f}")
            print(f"   {label:28} " + " ".join(cells))

    print("\n2. Least-squares operators, judged without time error (C = 0.85 + TDT)")
    for L in (4, 8):
        for t in tols:
            best = None
            for b100 in range(40, 300, 5):
                c = least_squares(L, b100 / 100)
                p = needed_ppw(c, r_s(c, 0.85), t, True)
                if best is None or p < best[0]:
                    best = (p, b100 / 100, c)
            p, band, c = best
            print(
                f"   order {2 * L:2d}, {100 * t:.1f}%: {p:.1f} ppw (band {band:.2f}, "
                f"sum|c| {sum(abs(x) for x in c):.4f}), work {work(p, c, 0.85, L) / ref_work[t]:.2f}"
            )

    print("\n3. Interpolating staggered values onto a node: amplitude lost (receiver; twice")
    print("   that for a source and receiver), at kh = 2 pi / ppw")
    for L in (1, 2, 3, 4):
        w = midpoint_weights(L)
        cells = []
        for ppw in (3.5, 4.5, 4.5 * 2.5):
            kept = sum(2 * wn * math.cos((n - 0.5) * 2 * math.pi / ppw) for n, wn in enumerate(w, 1))
            cells.append(f"{ppw:4.1f} ppw: {100 * (1 - kept):5.2f}%")
        print(f"   {2 * L}-point: " + ", ".join(cells))
