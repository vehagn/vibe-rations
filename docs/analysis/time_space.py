"""Is the time step what limits long stencils? And coefficients fitted to the time-space
dispersion relation (after Liu & Sen 2011). Run: python3 docs/analysis/time_space.py
(pure Python, about 10 s). Results: docs/research-log.md, 2026-11-20, "Stencil order and
time step"."""
import math
from dispersion import taylor, solve, errors, needed_ppw, AXIS, DIAGONAL

VPVS = math.sqrt(3)

def r_s(c, courant):
    return courant / (math.sqrt(3) * sum(abs(x) for x in c)) / VPVS

def time_space(L, r, beta, samples=300):
    """Least squares on sum_n c_n sin((2n-1) xi_i/2) so that along the axis and the diagonal
    sin(w dt/2) = r |K| gives w = v k exactly, for xi = k h in (0, beta]; relative weighting."""
    rows, rhs = [], []
    for s in range(1, samples + 1):
        xi = beta * s / samples
        target = math.sin(r * xi / 2) / r
        w = 1 / target
        # axis: S = sum c_n sin((2n-1) xi/2)
        rows.append([w * math.sin((2*n-1) * xi / 2) for n in range(1, L + 1)]); rhs.append(w * target)
        # diagonal: S = sqrt(3) sum c_n sin((2n-1) xi/(2 sqrt 3))
        rows.append([w * math.sqrt(3) * math.sin((2*n-1) * xi / (2*math.sqrt(3))) for n in range(1, L + 1)]); rhs.append(w * target)
    AtA = [[sum(r_[i]*r_[j] for r_ in rows) for j in range(L)] for i in range(L)]
    Atb = [sum(r_[i]*y for r_, y in zip(rows, rhs)) for i in range(L)]
    return solve(AtA, Atb)

print("1) Taylor operators at smaller time steps: ppw for S group error 0.5% / 0.1%")
print(f"{'operator':12}" + "".join(f"  courant {c:4}" for c in (0.85, 0.5, 0.25)))
for L in (2, 4, 8):
    c = taylor(L)
    cells = []
    for courant in (0.85, 0.5, 0.25):
        r = r_s(c, courant)
        cells.append(f"{needed_ppw(c, r, 0.005):5.1f}/{needed_ppw(c, r, 0.001):4.1f}")
    print(f"{'Taylor '+str(2*L)+'pt':12}" + "".join(f"  {x:>12}" for x in cells))

print("\n2) Time-space optimised for S at courant 0.85 (band chosen for 0.5% group error)")
print(f"{'operator':12} {'band':>5} {'sum|c|':>7} | S ppw 0.5%/0.1% | at half the S speed | P on the same grid")
for L in (2, 4, 8):
    best = None
    for b100 in range(40, 300, 5):
        beta = b100 / 100
        c = time_space(L, 0.2, beta)  # first guess of r
        for _ in range(3):  # r depends on sum|c|: iterate to consistency
            c = time_space(L, r_s(c, 0.85), beta)
        r = r_s(c, 0.85)
        if any(math.isinf(errors(c, 0.1, d, r)[1]) for d in (AXIS, DIAGONAL)):
            continue
        ppw = needed_ppw(c, r, 0.005)
        if best is None or ppw < best[0]:
            best = (ppw, beta, c, r)
    ppw, beta, c, r = best
    slow = needed_ppw(c, r / 2, 0.005)  # a slower S wave: half the Courant number
    # P travels sqrt(3) faster: its Courant number is sqrt(3) r and it has sqrt(3) more ppw.
    p_ok = needed_ppw(c, r * VPVS, 0.005)
    print(f"{'TS '+str(2*L)+'pt':12} {beta:5.2f} {sum(abs(x) for x in c):7.4f} | {ppw:5.1f} / {needed_ppw(c, r, 0.001):4.1f}  | {slow:5.1f} ppw needed   | needs {p_ok:4.1f} (has {ppw*VPVS:4.1f})")
    print("   c =", ", ".join(f"{x:.10g}" for x in c))
