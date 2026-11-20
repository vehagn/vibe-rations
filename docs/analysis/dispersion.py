"""Von Neumann dispersion analysis of the 3D staggered-grid velocity-stress scheme.

For a homogeneous isotropic medium the staggered scheme with derivative coefficients c_n
(half-length L) and leapfrog time stepping has, for P and S alike (speed v),

    sin(w dt / 2) = r * |K|,   K_i = sum_n c_n sin((2n - 1) k_i h / 2),   r = v dt / h,

so the numerical phase velocity is w/|k| and the group velocity dw/d|k| along a direction.
Run: python3 docs/analysis/dispersion.py (pure Python, about 15 s). Prints the points per
wavelength each operator needs for a given S group-velocity error, Taylor and least-squares
optimised. Results: docs/research-log.md, 2026-11-20, "Stencil order and time step".
"""

import math
from fractions import Fraction


def taylor(L):
    """Staggered Taylor coefficients: sum c_n (2n-1)^(2m-1) = [m == 1], m = 1..L."""
    A = [[Fraction((2 * n - 1) ** (2 * m - 1)) for n in range(1, L + 1)] for m in range(1, L + 1)]
    b = [Fraction(1 if m == 0 else 0) for m in range(L)]
    return [float(x) for x in solve(A, b)]


def solve(A, b):
    n = len(b)
    M = [row[:] + [b[i]] for i, row in enumerate(A)]
    for c in range(n):
        p = max(range(c, n), key=lambda r: abs(M[r][c]))
        M[c], M[p] = M[p], M[c]
        for r in range(n):
            if r != c and M[r][c] != 0:
                f = M[r][c] / M[c][c]
                M[r] = [x - f * y for x, y in zip(M[r], M[c])]
    return [M[i][n] / M[i][i] for i in range(n)]


def least_squares(L, beta, samples=400):
    """Coefficients minimising sum over xi in (0, beta] of (2 sum c_n sin((2n-1) xi/2) - xi)^2,
    weighted to relative error, with the consistency condition sum c_n (2n-1) = 1 imposed by
    eliminating c_1."""
    # Unknowns c_2..c_L; c_1 = 1 - sum_{n>=2} c_n (2n-1).
    rows, rhs = [], []
    for s in range(1, samples + 1):
        xi = beta * s / samples
        w = 1.0 / xi  # relative error
        base = 2 * math.sin(xi / 2)  # from c_1 = 1
        rows.append([w * (2 * math.sin((2 * n - 1) * xi / 2) - (2 * n - 1) * 2 * math.sin(xi / 2)) for n in range(2, L + 1)])
        rhs.append(w * (xi - base))
    m = L - 1
    AtA = [[sum(r[i] * r[j] for r in rows) for j in range(m)] for i in range(m)]
    Atb = [sum(r[i] * y for r, y in zip(rows, rhs)) for i in range(m)]
    rest = solve(AtA, Atb)
    c1 = 1 - sum(c * (2 * n - 1) for n, c in zip(range(2, L + 1), rest))
    return [c1] + rest


def K(c, xi):
    return sum(cn * math.sin((2 * n - 1) * xi / 2) for n, cn in enumerate(c, start=1)) * 2


def omega(c, kh, direction, r):
    """w dt for |k| h = kh along `direction` (unit vector), Courant number r = v dt / h.
    Returns None when the scheme is unstable for this mode."""
    s = math.sqrt(sum((K(c, kh * d) / 2) ** 2 for d in direction))
    x = r * s
    if x > 1:
        return None
    return 2 * math.asin(x)


def errors(c, kh, direction, r):
    """(phase, group) velocity errors relative to the true speed at |k| h = kh."""
    w = omega(c, kh, direction, r)
    if w is None:
        return math.inf, math.inf
    phase = w / (r * kh) - 1
    e = 1e-6
    group = (omega(c, kh + e, direction, r) - omega(c, kh - e, direction, r)) / (2 * e) / r - 1
    return phase, group


AXIS = (1.0, 0.0, 0.0)
DIAGONAL = (1 / math.sqrt(3),) * 3


def needed_ppw(c, r, tol, which=1):
    """Fewest points per wavelength such that every wavenumber up to 2 pi / ppw has
    |error| <= tol along the axis and the diagonal."""
    lo = 2.0
    for ppw10 in range(20, 400):
        ppw = ppw10 / 10
        khmax = 2 * math.pi / ppw
        ok = True
        for s in range(1, 61):
            kh = khmax * s / 60
            for d in (AXIS, DIAGONAL):
                if abs(errors(c, kh, d, r)[which]) > tol:
                    ok = False
                    break
            if not ok:
                break
        if ok:
            return ppw
    return math.inf


if __name__ == "__main__":
    vpvs = math.sqrt(3)  # Poisson solid; S is the slow wave that sets the grid
    courant = 0.85  # our default fraction of the stability limit
    schemes = []
    for L in (1, 2, 3, 4, 6, 8):
        schemes.append((f"Taylor {2*L}th, {2*L}-point", taylor(L)))
    # Optimised: pick the band that minimises the ppw for 0.5% group error.
    for L in (2, 4, 6, 8):
        best = None
        for b100 in range(30, 300, 5):
            beta = b100 / 100
            c = least_squares(L, beta)
            sc = sum(abs(x) for x in c)
            r_s = courant / (math.sqrt(3) * sc) / vpvs
            ppw = needed_ppw(c, r_s, 0.005, which=1)
            if best is None or ppw < best[0]:
                best = (ppw, beta, c)
        schemes.append((f"Optimised {2*L}-point (band {best[1]:.2f})", best[2]))

    print(f"S waves at vp/vs = {vpvs:.2f}, dt at {courant} of the stability limit\n")
    print(f"{'operator':34} {'sum|c|':>7} {'dt':>6} | ppw for S group-velocity error  | phase 0.5%")
    print(f"{'':34} {'':>7} {'':>6} |   1%   0.5%   0.1%")
    for name, c in schemes:
        sc = sum(abs(x) for x in c)
        r_p = courant / (math.sqrt(3) * sc)  # P waves at the stability fraction
        r_s = r_p / vpvs
        g = [needed_ppw(c, r_s, t, which=1) for t in (0.01, 0.005, 0.001)]
        p = needed_ppw(c, r_s, 0.005, which=0)
        dt = 1 / sc / (1 / (9 / 8 + 1 / 24))  # relative to our 4th-order dt
        print(f"{name:34} {sc:7.4f} {dt:6.3f} | {g[0]:5.1f} {g[1]:6.1f} {g[2]:6.1f} | {p:6.1f}")
    print()
    for name, c in schemes:
        print(f"{name:34} " + ", ".join(f"{x:.10g}" for x in c))
