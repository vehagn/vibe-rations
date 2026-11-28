//! Staggered-grid first derivatives along one axis of a padded array, a row of nodes at a time.
//!
//! A [`Stencil`] of half-length L approximates the derivative half a node away from a node as
//! (1/h) Σₙ cₙ [f(x + (n - ½)h) - f(x - (n - ½)h)], n = 1..L.
//!
//! `start` is the padded index of the row's first node and `stride` the distance in memory
//! between neighbours along the axis being differentiated. Both methods fill `out` with the
//! derivative for `out.len()` consecutive nodes (x fastest), each placed half a node forward or
//! backward of the node it is computed for. They read L nodes either way, so the grid's halo
//! must be at least L wide.

#![allow(clippy::needless_range_loop)]

use std::f64::consts::PI;
use std::sync::Mutex;

use anyhow::{Result, ensure};

use crate::simd::Arithmetic;

/// The longest stencil offered: 16 points.
pub const MAX_ORDER: usize = 16;

/// A staggered first-derivative operator, with the weights that interpolate staggered values
/// onto a node over the same span.
#[derive(Clone, Debug, PartialEq)]
pub struct Stencil {
    coefficients: Vec<f64>,
    /// The same, in single precision for the kernels.
    single: Vec<f32>,
    /// Weights w₁..w_L of the values at ±(n - ½) nodes that interpolate onto the node.
    interpolation: Vec<f64>,
    /// For optimised coefficients, the largest k·h they were designed for.
    band: Option<f64>,
}

impl Stencil {
    /// The Taylor operator of the given order (2, 4, ..., 16): exact for polynomials up to
    /// degree `order`, with half-length `order / 2`. Order 4 is Levander's (1988) (9/8, -1/24).
    pub fn taylor(order: usize) -> Result<Self> {
        ensure!(
            (2..=MAX_ORDER).contains(&order) && order.is_multiple_of(2),
            "the stencil order must be even and from 2 to {MAX_ORDER}, not {order}"
        );
        let l = order / 2;
        // Σ cₙ (2n-1)^(2k-1) = [k = 1] for k = 1..L. With xₙ = (2n-1)², aₙ = cₙ (2n-1) solves a
        // Vandermonde system whose answer is the Lagrange basis at zero: aₙ = Π_{m≠n} xₘ/(xₘ-xₙ).
        let x = |n: usize| ((2 * n - 1) * (2 * n - 1)) as f64;
        let coefficients: Vec<f64> = (1..=l)
            .map(|n| {
                let a: f64 = (1..=l)
                    .filter(|&m| m != n)
                    .map(|m| x(m) / (x(m) - x(n)))
                    .product();
                a / (2 * n - 1) as f64
            })
            .collect();
        Ok(Self::from_parts(coefficients, midpoint_weights(l), None))
    }

    /// Coefficients chosen for wavenumbers up to `band` = k·h < π rather than for the
    /// longest waves: those that minimise the largest phase error per node,
    /// | |K| - |k|·h |, over every direction of propagation in 3D and every |k|·h up to the
    /// band, weighted by (|k|·h/band)^-2 so that longer waves get less of it, with the
    /// derivative still exact for linear functions. With leapfrog's time error removed that
    /// error is all there is: over a distance D it turns a wave's phase by D/h times as much.
    ///
    /// Spreading the error over a band too wide for the stencil would make every wave worse,
    /// the longest included. So the (weighted) error is capped at `tolerance` per node: if the
    /// band can't be had within it, the band shrinks to the widest that can ([`Stencil::band`]
    /// says which). The interpolation weights are chosen the same way over that band, for the
    /// largest amplitude error.
    ///
    /// Designed by Lawson's algorithm (Rice & Usow 1968) on the linearised error, in the
    /// manner of Holberg (1987), who bounded the group-velocity error the same way. Order 2
    /// has nothing to choose and is Taylor's.
    ///
    /// The design takes a good fraction of a second, and an inversion sets up many runs with
    /// the same operator, so designs are remembered (a few, the latest kept).
    pub fn optimised(order: usize, band: f64, tolerance: f64) -> Result<Self> {
        /// Designs remembered: (order, band and tolerance bits) and the stencil.
        type Designs = Vec<((usize, u64, u64), Stencil)>;
        static DESIGNS: Mutex<Designs> = Mutex::new(Vec::new());
        const REMEMBERED: usize = 8;
        let key = (order, band.to_bits(), tolerance.to_bits());
        if let Some((_, s)) = DESIGNS.lock().unwrap().iter().find(|(k, _)| *k == key) {
            return Ok(s.clone());
        }
        let s = Self::design(order, band, tolerance)?;
        let mut designs = DESIGNS.lock().unwrap();
        if designs.len() == REMEMBERED {
            designs.remove(0);
        }
        designs.push((key, s.clone()));
        Ok(s)
    }

    /// [`Stencil::optimised`], designed afresh.
    fn design(order: usize, band: f64, tolerance: f64) -> Result<Self> {
        let taylor = Self::taylor(order)?;
        ensure!(
            band > 0.0 && band < PI,
            "the band, k·h, must be in (0, π), not {band}"
        );
        ensure!(tolerance > 0.0, "the tolerance must be positive");
        let l = order / 2;
        if l == 1 {
            return Ok(taylor);
        }
        // Unknowns c₂..c_L; c₁ = 1 - Σ (2n-1) cₙ keeps the derivative exact for linear functions.
        let design = |b: f64| minimax(&derivative_rows(l, b));
        let (mut rest, worst) = design(band);
        let mut band = band;
        if worst > tolerance {
            // The widest band within the tolerance, by bisection; Taylor's if none is.
            let (mut lo, mut hi) = (0.0, band);
            for _ in 0..16 {
                let mid = 0.5 * (lo + hi);
                let (x, e) = design(mid);
                if e <= tolerance {
                    (lo, rest) = (mid, x);
                } else {
                    hi = mid;
                }
            }
            if lo == 0.0 {
                return Ok(taylor);
            }
            band = lo;
        }
        let mut coefficients = vec![
            1.0 - (2..=l)
                .zip(&rest)
                .map(|(n, c)| (2 * n - 1) as f64 * c)
                .sum::<f64>(),
        ];
        coefficients.extend(rest);
        // Unknowns w₂..w_L; w₁ = ½ - Σ wₙ keeps a constant exact.
        let (rest, _) = minimax(&interpolation_rows(l, band));
        let mut interpolation = vec![0.5 - rest.iter().sum::<f64>()];
        interpolation.extend(rest);
        Ok(Self::from_parts(coefficients, interpolation, Some(band)))
    }

    fn from_parts(coefficients: Vec<f64>, interpolation: Vec<f64>, band: Option<f64>) -> Self {
        let single = coefficients.iter().map(|&c| c as f32).collect();
        Self {
            coefficients,
            single,
            interpolation,
            band,
        }
    }

    /// For optimised coefficients, the largest k·h they were designed for.
    pub fn band(&self) -> Option<f64> {
        self.band
    }

    /// Weights w₁..w_L for the values at ±(n - ½) nodes that interpolate onto the node:
    /// Lagrange's for Taylor coefficients, optimised over the band for optimised ones.
    pub fn interpolation(&self) -> &[f64] {
        &self.interpolation
    }

    /// Weights for the 2L values at offsets j - δ (j = 1 - L, ..., L, in nodes) from a point
    /// that lies δ ∈ (0, 1) past the nearest value below it: how a point between two values of
    /// a lattice is read from the 2L values around it. Lagrange's for Taylor coefficients; for
    /// optimised ones chosen as their interpolation onto a node is, minimising the largest
    /// complex error |Σ wⱼ e^{i x (j - δ)} - 1| (its real and imaginary parts) for k·h up to
    /// the band, weighted the same way, and exact for constants and linear functions. At δ = ½
    /// they are [`Stencil::interpolation`]'s. Designs are remembered.
    pub fn interpolation_at(&self, delta: f64) -> Vec<f64> {
        let l = self.half_length();
        let points: Vec<f64> = (0..2 * l)
            .map(|k| k as f64 - (l - 1) as f64 - delta)
            .collect();
        if delta == 0.5 {
            // p = ±(n - ½): weight wₙ.
            return points
                .iter()
                .map(|p| self.interpolation[(p.abs() - 0.5).round() as usize])
                .collect();
        }
        let Some(band) = self.band else {
            return lagrange_at_zero(&points);
        };
        type Designs = Vec<((usize, u64, u64), Vec<f64>)>;
        static DESIGNS: Mutex<Designs> = Mutex::new(Vec::new());
        const REMEMBERED: usize = 256;
        let key = (l, band.to_bits(), delta.to_bits());
        if let Some((_, w)) = DESIGNS.lock().unwrap().iter().find(|(k, _)| *k == key) {
            return w.clone();
        }
        let w = interpolation_design(&points, band);
        let mut designs = DESIGNS.lock().unwrap();
        if designs.len() == REMEMBERED {
            designs.remove(0);
        }
        designs.push((key, w.clone()));
        w
    }

    /// The largest phase error per node, | |K| - |k|h |, over every direction (sampled densely)
    /// and every |k|·h up to `band`, computed exactly rather than linearised.
    pub fn phase_error(&self, band: f64) -> f64 {
        let directions = directions(12);
        let mut worst = 0.0f64;
        for s in 1..=200 {
            let xi = band * s as f64 / 200.0;
            for d in &directions {
                let k: f64 = d
                    .iter()
                    .map(|&di| self.response(xi * di).powi(2))
                    .sum::<f64>()
                    .sqrt();
                worst = worst.max((k - xi).abs());
            }
        }
        worst
    }

    /// The largest relative amplitude error of the interpolation onto a node, up to `band`.
    pub fn interpolation_error(&self, band: f64) -> f64 {
        (1..=400)
            .map(|s| {
                let x = band * s as f64 / 400.0;
                let r: f64 = self
                    .interpolation
                    .iter()
                    .enumerate()
                    .map(|(n, w)| 2.0 * w * ((n as f64 + 0.5) * x).cos())
                    .sum();
                (r - 1.0).abs()
            })
            .fold(0.0, f64::max)
    }

    /// K(x) = 2 Σ cₙ sin((2n-1) x/2): what the operator makes of a wave with k·h = x along
    /// its axis, in units of 1/h.
    pub fn response(&self, x: f64) -> f64 {
        self.coefficients
            .iter()
            .enumerate()
            .map(|(n, c)| 2.0 * c * ((2 * n + 1) as f64 * x / 2.0).sin())
            .sum()
    }

    /// Its order of accuracy.
    pub fn order(&self) -> usize {
        2 * self.coefficients.len()
    }

    /// Nodes it reaches on each side; the halo must be at least this wide.
    pub fn half_length(&self) -> usize {
        self.coefficients.len()
    }

    pub fn coefficients(&self) -> &[f64] {
        &self.coefficients
    }

    /// The coefficients in single precision, as [`Stencil::forward`] and
    /// [`Stencil::backward`] apply them.
    pub fn single(&self) -> &[f32] {
        &self.single
    }

    /// Σ |cₙ|, which sets the stability limit: dt ≤ 1 / (vp Σ|cₙ| √(Σ 1/h²)).
    pub fn sum_abs(&self) -> f64 {
        self.coefficients.iter().map(|c| c.abs()).sum()
    }

    /// The derivative half a node forward of each node, at p + s/2. Always inlined, so the
    /// solver's AVX2 sweep (crate::simd) compiles it with AVX2.
    #[inline(always)]
    pub fn forward<A: Arithmetic>(
        &self,
        out: &mut [f32],
        f: &[f32],
        start: usize,
        stride: usize,
        inv_h: f32,
    ) {
        let len = out.len();
        for (k, &c) in self.single.iter().enumerate() {
            let a = c * inv_h;
            let hi = &f[start + (k + 1) * stride..start + (k + 1) * stride + len];
            let lo = &f[start - k * stride..start - k * stride + len];
            accumulate::<A>(out, a, hi, lo, k == 0);
        }
    }

    /// The derivative half a node backward of each node, at p - s/2.
    #[inline(always)]
    pub fn backward<A: Arithmetic>(
        &self,
        out: &mut [f32],
        f: &[f32],
        start: usize,
        stride: usize,
        inv_h: f32,
    ) {
        let len = out.len();
        for (k, &c) in self.single.iter().enumerate() {
            let a = c * inv_h;
            let hi = &f[start + k * stride..start + k * stride + len];
            let lo = &f[start - (k + 1) * stride..start - (k + 1) * stride + len];
            accumulate::<A>(out, a, hi, lo, k == 0);
        }
    }
}

/// Lagrange's weights at 0 for values at `points`: exact for polynomials of degree below
/// their number.
pub fn lagrange_at_zero(points: &[f64]) -> Vec<f64> {
    points
        .iter()
        .enumerate()
        .map(|(i, &x)| {
            points
                .iter()
                .enumerate()
                .filter(|&(j, _)| j != i)
                .map(|(_, &y)| -y / (x - y))
                .product::<f64>()
        })
        .collect()
}

/// Weights that differentiate at 0 the polynomial through values at `points`: the
/// derivatives at 0 of Lagrange's basis polynomials.
pub fn lagrange_slope_at_zero(points: &[f64]) -> Vec<f64> {
    (0..points.len())
        .map(|i| {
            (0..points.len())
                .filter(|&k| k != i)
                .map(|k| {
                    let others: f64 = (0..points.len())
                        .filter(|&j| j != i && j != k)
                        .map(|j| -points[j] / (points[i] - points[j]))
                        .product();
                    others / (points[i] - points[k])
                })
                .sum()
        })
        .collect()
}

/// Interpolation weights at 0 for values at `points` (2L of them, the two nearest 0 at
/// indices L - 1 and L) minimising the weighted largest complex error over the band, exact for
/// constants and linear functions: Lagrange's weights plus a combination of vectors that
/// change neither Σ w nor Σ w p, each moving weight from the two nearest points to one other.
fn interpolation_design(points: &[f64], band: f64) -> Vec<f64> {
    let base = lagrange_at_zero(points);
    let (a, b) = (points.len() / 2 - 1, points.len() / 2);
    let basis: Vec<Vec<f64>> = (0..points.len())
        .filter(|&k| k != a && k != b)
        .map(|k| {
            let mut v = vec![0.0; points.len()];
            v[k] = 1.0;
            v[b] = (points[a] - points[k]) / (points[b] - points[a]);
            v[a] = -1.0 - v[b];
            v
        })
        .collect();
    let mut rows = Vec::with_capacity(400);
    for s in 1..=200 {
        let x = band * s as f64 / 200.0;
        let sum = |w: &[f64], f: fn(f64) -> f64| {
            w.iter().zip(points).map(|(w, p)| w * f(x * p)).sum::<f64>()
        };
        // Real part: Σ wⱼ cos(x pⱼ) - 1; imaginary: Σ wⱼ sin(x pⱼ).
        for (f, target) in [(f64::cos as fn(f64) -> f64, 1.0), (f64::sin, 0.0)] {
            let row = basis.iter().map(|v| sum(v, f)).collect();
            rows.push(weighted(row, target - sum(&base, f), x / band));
        }
    }
    let (z, _) = minimax_complex(&rows);
    let mut w = base;
    for (v, z) in basis.iter().zip(&z) {
        for (w, v) in w.iter_mut().zip(v) {
            *w += z * v;
        }
    }
    w
}

/// Lagrange weights at 0 for values at ±½, ±3/2, ..., ±(m - ½): the value at k - ½ (and its
/// mirror) gets weight k - 1.
pub fn midpoint_weights(m: usize) -> Vec<f64> {
    let points: Vec<f64> = (1..=m).map(|k| k as f64 - 0.5).collect();
    points
        .iter()
        .map(|&x| {
            points
                .iter()
                .flat_map(|&y| [y, -y])
                .filter(|&y| y != x)
                .map(|y| -y / (x - y))
                .product::<f64>()
        })
        .collect()
}

/// Unit vectors with d₁ ≥ d₂ ≥ d₃ ≥ 0, which stand for every direction: the errors are the
/// same under permutations and sign changes. From the axis to the diagonal, (1, i/m, j/m)
/// normalised for 0 ≤ j ≤ i ≤ m.
fn directions(m: usize) -> Vec<[f64; 3]> {
    let mut out = Vec::new();
    for i in 0..=m {
        for j in 0..=i {
            let v = [1.0, i as f64 / m as f64, j as f64 / m as f64];
            let norm = v.iter().map(|x| x * x).sum::<f64>().sqrt();
            out.push(v.map(|x| x / norm));
        }
    }
    out
}

/// A linear function of the unknowns and its target: the error is a·x - b.
type Row = (Vec<f64>, f64);

/// The phase error per node, linearised: |K| - |k|h ≈ Σᵢ dᵢ δ(ξ dᵢ) with δ(x) = K(x) - x,
/// for c₂..c_L with c₁ eliminated.
fn derivative_rows(l: usize, band: f64) -> Vec<Row> {
    // K(x) = 2 sin(x/2) + Σ_{n≥2} cₙ gₙ(x), gₙ(x) = 2 sin((2n-1)x/2) - 2(2n-1) sin(x/2).
    let g = |n: usize, x: f64| {
        let m = (2 * n - 1) as f64;
        2.0 * (m * x / 2.0).sin() - 2.0 * m * (x / 2.0).sin()
    };
    let mut rows = Vec::new();
    for d in directions(6) {
        for s in 1..=48 {
            let xi = band * s as f64 / 48.0;
            let a = (2..=l)
                .map(|n| d.iter().map(|&di| di * g(n, xi * di)).sum())
                .collect();
            let b = -d
                .iter()
                .map(|&di| di * (2.0 * (xi * di / 2.0).sin() - xi * di))
                .sum::<f64>();
            rows.push(weighted(a, b, xi / band));
        }
    }
    rows
}

/// The interpolation's amplitude error, R(x) - 1 with R(x) = Σ 2wₙ cos((n - ½)x), for
/// w₂..w_L with w₁ eliminated.
fn interpolation_rows(l: usize, band: f64) -> Vec<Row> {
    (1..=200)
        .map(|s| {
            let x = band * s as f64 / 200.0;
            let a = (2..=l)
                .map(|n| 2.0 * (((n as f64 - 0.5) * x).cos() - (x / 2.0).cos()))
                .collect();
            weighted(a, 1.0 - (x / 2.0).cos(), x / band)
        })
        .collect()
}

/// How fast the allowed error falls towards long waves: as (k/k_band)^p. With p = 0 (a plain
/// minimax) every wavenumber in the band may carry the whole tolerance, and P waves and the
/// dominant frequency, well inside the band, came out worse than with Taylor's coefficients;
/// p = 2 keeps them as good while giving up little at the band's edge (research log).
const FLATNESS: i32 = 2;

/// A row scaled so that its allowed error is (k/k_band)^p of the tolerance at fraction
/// `f` = k/k_band of the band.
fn weighted(a: Vec<f64>, b: f64, f: f64) -> Row {
    let w = f.powi(-FLATNESS);
    (a.into_iter().map(|v| v * w).collect(), b * w)
}

/// The x minimising max |a·x - b| over the rows, by Lawson's algorithm: weighted least
/// squares, each row's weight multiplied by its error every round, which concentrates the
/// weight on the extremal rows (Rice & Usow 1968). Returns the best iterate and its largest
/// error.
fn minimax(rows: &[Row]) -> (Vec<f64>, f64) {
    lawson(rows, 1)
}

/// [`minimax`] for a complex error: the rows come in pairs, its real and imaginary parts at
/// one point, and the largest modulus √(re² + im²) is minimised; each pair shares a weight.
/// Taken part by part, the minimax isn't unique (any imaginary part below the largest real
/// one would do), and the design could pick weights with a needless phase error.
fn minimax_complex(rows: &[Row]) -> (Vec<f64>, f64) {
    lawson(rows, 2)
}

/// Lawson's algorithm on the errors of the rows in groups of `group`: one row each (their
/// absolute errors), or two (the moduli of complex errors, real and imaginary part).
fn lawson(rows: &[Row], group: usize) -> (Vec<f64>, f64) {
    let groups = rows.len() / group;
    let mut weights = vec![1.0 / groups as f64; groups];
    let mut best = (f64::INFINITY, Vec::new());
    for _ in 0..400 {
        let row_weights: Vec<f64> = weights
            .iter()
            .flat_map(|&w| std::iter::repeat_n(w, group))
            .collect();
        let x = weighted_least_squares(rows, &row_weights);
        let error = |(a, b): &Row| a.iter().zip(&x).map(|(p, q)| p * q).sum::<f64>() - b;
        let errors: Vec<f64> = rows
            .chunks(group)
            .map(|g| match g {
                [row] => error(row).abs(),
                [re, im] => error(re).hypot(error(im)),
                _ => unreachable!("groups of one or two rows"),
            })
            .collect();
        let worst = errors.iter().copied().fold(0.0, f64::max);
        if worst < best.0 {
            best = (worst, x);
        }
        let total: f64 = weights.iter().zip(&errors).map(|(w, e)| w * e).sum();
        if total == 0.0 {
            break;
        }
        for (w, e) in weights.iter_mut().zip(&errors) {
            *w *= e / total;
        }
    }
    (best.1, best.0)
}

/// The x minimising Σ wᵢ (aᵢ·x - bᵢ)², by Householder QR on columns scaled to unit norm.
fn weighted_least_squares(rows: &[Row], weights: &[f64]) -> Vec<f64> {
    let n = rows[0].0.len();
    let mut a: Vec<Vec<f64>> = rows
        .iter()
        .zip(weights)
        .map(|((r, _), w)| r.iter().map(|v| v * w.sqrt()).collect())
        .collect();
    let mut b: Vec<f64> = rows
        .iter()
        .zip(weights)
        .map(|((_, v), w)| v * w.sqrt())
        .collect();
    let scale: Vec<f64> = (0..n)
        .map(|j| {
            let s = a.iter().map(|r| r[j] * r[j]).sum::<f64>().sqrt();
            if s > 0.0 { 1.0 / s } else { 1.0 }
        })
        .collect();
    for r in &mut a {
        for (v, s) in r.iter_mut().zip(&scale) {
            *v *= s;
        }
    }
    let m = a.len();
    for k in 0..n {
        let norm = (k..m).map(|i| a[i][k] * a[i][k]).sum::<f64>().sqrt();
        if norm == 0.0 {
            continue;
        }
        let alpha = if a[k][k] > 0.0 { -norm } else { norm };
        let mut v: Vec<f64> = (k..m).map(|i| a[i][k]).collect();
        v[0] -= alpha;
        let vv: f64 = v.iter().map(|x| x * x).sum();
        if vv == 0.0 {
            continue;
        }
        for j in k..n {
            let f = 2.0 * (k..m).map(|i| v[i - k] * a[i][j]).sum::<f64>() / vv;
            for i in k..m {
                a[i][j] -= f * v[i - k];
            }
        }
        let f = 2.0 * (k..m).map(|i| v[i - k] * b[i]).sum::<f64>() / vv;
        for i in k..m {
            b[i] -= f * v[i - k];
        }
    }
    let mut x = vec![0.0; n];
    for k in (0..n).rev() {
        let s: f64 = (k + 1..n).map(|j| a[k][j] * x[j]).sum();
        x[k] = if a[k][k] != 0.0 {
            (b[k] - s) / a[k][k]
        } else {
            0.0
        };
    }
    x.iter().zip(&scale).map(|(x, s)| x * s).collect()
}

/// out = a (hi - lo), or out += a (hi - lo) after the first term, the product and the sum
/// rounded as `A` says.
#[inline(always)]
fn accumulate<A: Arithmetic>(out: &mut [f32], a: f32, hi: &[f32], lo: &[f32], first: bool) {
    let n = out.len();
    let (hi, lo) = (&hi[..n], &lo[..n]);
    if first {
        for i in 0..n {
            out[i] = a * (hi[i] - lo[i]);
        }
    } else {
        for i in 0..n {
            out[i] = A::product_plus(a, hi[i] - lo[i], out[i]);
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::simd::Separate;

    /// A remembered design is the design, bit for bit, and comes back far faster.
    #[test]
    fn remembers_designs() {
        // A band no other test designs for, so the first call designs it; wide enough that
        // the design searches for the band it can have, as runs at few points per wavelength do.
        let band = 0.85 * PI + 1e-9;
        let t = std::time::Instant::now();
        let first = Stencil::optimised(12, band, 1e-4).unwrap();
        let designing = t.elapsed();
        let t = std::time::Instant::now();
        let again = Stencil::optimised(12, band, 1e-4).unwrap();
        let remembering = t.elapsed();
        assert_eq!(first, again);
        assert_eq!(first, Stencil::design(12, band, 1e-4).unwrap());
        assert_ne!(first, Stencil::optimised(12, band - 0.1, 1e-4).unwrap());
        println!("designing {designing:?}, remembering {remembering:?}");
        assert!(
            remembering * 10 < designing,
            "{remembering:?} vs {designing:?}"
        );
    }

    #[test]
    fn taylor_coefficients_are_the_known_ones() {
        assert_eq!(Stencil::taylor(2).unwrap().coefficients(), &[1.0]);
        assert_eq!(
            Stencil::taylor(4).unwrap().coefficients(),
            &[9.0 / 8.0, -1.0 / 24.0]
        );
        let eighth = [
            1225.0 / 1024.0,
            -245.0 / 3072.0,
            49.0 / 5120.0,
            -5.0 / 7168.0,
        ];
        for (a, b) in Stencil::taylor(8)
            .unwrap()
            .coefficients()
            .iter()
            .zip(eighth)
        {
            assert!((a - b).abs() < 1e-15, "{a} vs {b}");
        }
        assert!((Stencil::taylor(4).unwrap().sum_abs() - (9.0 / 8.0 + 1.0 / 24.0)).abs() < 1e-15);
        assert!(Stencil::taylor(5).is_err());
        assert!(Stencil::taylor(18).is_err());
        assert!(Stencil::taylor(0).is_err());
    }

    /// The defining conditions: Σ cₙ (2n-1)^(2k-1) = [k = 1] for k = 1..L.
    #[test]
    fn coefficients_cancel_every_error_term_up_to_their_order() {
        for order in (2..=MAX_ORDER).step_by(2) {
            let s = Stencil::taylor(order).unwrap();
            for k in 1..=s.half_length() {
                let terms: Vec<f64> = s
                    .coefficients()
                    .iter()
                    .enumerate()
                    .map(|(i, c)| c * ((2 * i + 1) as f64).powi(2 * k as i32 - 1))
                    .collect();
                let moment: f64 = terms.iter().sum();
                // The terms grow to ~1e15 at order 16: compare with their size.
                let scale: f64 = terms.iter().map(|t| t.abs()).sum();
                let want = if k == 1 { 1.0 } else { 0.0 };
                assert!(
                    (moment - want).abs() < 1e-13 * scale,
                    "order {order}, k {k}: {moment}"
                );
            }
        }
    }

    /// Every order puts the derivative half a node forward or backward of the right node: a
    /// quadratic, which every order differentiates exactly, shows the position.
    #[test]
    fn places_the_derivative_half_a_node_away() {
        let h = 0.5f32;
        let f: Vec<f32> = (0..40)
            .map(|i| {
                let x = i as f32 * h;
                x * x - 2.0 * x
            })
            .collect();
        let df = |x: f32| 2.0 * x - 2.0;
        for order in (2..=MAX_ORDER).step_by(2) {
            let s = Stencil::taylor(order).unwrap();
            let mut out = [0.0f32; 4];
            s.forward::<Separate>(&mut out, &f, 10, 1, 1.0 / h);
            for (i, v) in out.iter().enumerate() {
                let x = (10 + i) as f32 * h + h / 2.0;
                assert!(
                    (v - df(x)).abs() < 1e-3,
                    "order {order} forward at {x}: {v}"
                );
            }
            s.backward::<Separate>(&mut out, &f, 10, 1, 1.0 / h);
            for (i, v) in out.iter().enumerate() {
                let x = (10 + i) as f32 * h - h / 2.0;
                assert!(
                    (v - df(x)).abs() < 1e-3,
                    "order {order} backward at {x}: {v}"
                );
            }
        }
    }

    /// Optimised coefficients beat Taylor's over their band, in every direction, and so does
    /// their interpolation; they stay exact for linear functions and constants. (Weighted
    /// towards long waves, they don't minimise the plain largest error, so the margin varies.)
    #[test]
    fn optimised_coefficients_beat_taylor_over_their_band() {
        for (order, band) in [(4, 1.5), (8, 1.8), (8, 2.2), (12, 2.0), (16, 2.4)] {
            let t = Stencil::taylor(order).unwrap();
            let o = Stencil::optimised(order, band, 1.0).unwrap();
            let (et, eo) = (t.phase_error(band), o.phase_error(band));
            let (it, io) = (t.interpolation_error(band), o.interpolation_error(band));
            eprintln!(
                "order {order}, band {band}: phase error per node {et:.2e} -> {eo:.2e}, \
                 interpolation {it:.2e} -> {io:.2e}, sum|c| {:.4} -> {:.4}",
                t.sum_abs(),
                o.sum_abs()
            );
            assert!(eo < 0.5 * et, "order {order}: {eo} vs {et}");
            assert!(io < it, "order {order}: {io} vs {it}");
            let consistency: f64 = o
                .coefficients()
                .iter()
                .enumerate()
                .map(|(n, c)| (2 * n + 1) as f64 * c)
                .sum();
            assert!((consistency - 1.0).abs() < 1e-12);
            assert!((2.0 * o.interpolation().iter().sum::<f64>() - 1.0).abs() < 1e-12);
            assert_eq!(o.band(), Some(band));
        }
        assert!(Stencil::optimised(8, 3.5, 1.0).is_err());
        assert_eq!(
            Stencil::optimised(2, 1.0, 1.0).unwrap(),
            Stencil::taylor(2).unwrap()
        );
    }

    /// Interpolating from any offset: Lagrange's for Taylor coefficients; for optimised ones
    /// exact for constants and linear functions and better than Lagrange's over the band, and
    /// at δ = ½ what the design for a node gives.
    #[test]
    fn interpolates_from_any_offset() {
        let error = |w: &[f64], delta: f64, band: f64| {
            let l = w.len() / 2;
            (1..=400)
                .map(|s| {
                    let x = band * s as f64 / 400.0;
                    let (re, im) = w.iter().enumerate().fold((0.0, 0.0), |(re, im), (j, w)| {
                        let p = x * (j as f64 - (l - 1) as f64 - delta);
                        (re + w * p.cos(), im + w * p.sin())
                    });
                    ((re - 1.0f64).powi(2) + im * im).sqrt()
                })
                .fold(0.0, f64::max)
        };
        for order in [4, 8, 16] {
            let t = Stencil::taylor(order).unwrap();
            let l = order / 2;
            for delta in [0.1, 0.37, 0.5, 0.81] {
                let points: Vec<f64> = (0..2 * l)
                    .map(|k| k as f64 - (l - 1) as f64 - delta)
                    .collect();
                // At ½ the stencil's own weights, which agree to round-off.
                let gap = t
                    .interpolation_at(delta)
                    .iter()
                    .zip(lagrange_at_zero(&points))
                    .map(|(a, b)| (a - b).abs())
                    .fold(0.0, f64::max);
                assert!(gap < 1e-15, "order {order}, δ {delta}: {gap}");
            }
        }
        for (order, band) in [(4, 1.5), (8, 1.8), (16, 2.4)] {
            let o = Stencil::optimised(order, band, 1.0).unwrap();
            let t = Stencil::taylor(order).unwrap();
            let l = order / 2;
            for delta in [0.1, 0.37, 0.81] {
                let w = o.interpolation_at(delta);
                let p: Vec<f64> = (0..2 * l)
                    .map(|k| k as f64 - (l - 1) as f64 - delta)
                    .collect();
                assert!((w.iter().sum::<f64>() - 1.0).abs() < 1e-12);
                assert!(w.iter().zip(&p).map(|(w, p)| w * p).sum::<f64>().abs() < 1e-12);
                let (eo, et) = (
                    error(&w, delta, band),
                    error(&t.interpolation_at(delta), delta, band),
                );
                eprintln!("order {order}, δ {delta}: largest error {et:.2e} -> {eo:.2e}");
                assert!(eo < et, "order {order}, δ {delta}: {eo} vs {et}");
            }
            // The general design at δ = ½ comes to the node's.
            let p: Vec<f64> = (0..2 * l)
                .map(|k| k as f64 - (l - 1) as f64 - 0.5)
                .collect();
            let general = interpolation_design(&p, band);
            let node = o.interpolation_at(0.5);
            let gap = general
                .iter()
                .zip(&node)
                .map(|(a, b)| (a - b).abs())
                .fold(0.0, f64::max);
            eprintln!("order {order}: the general design at ½ is off the node's by {gap:.1e}");
            assert!(gap < 1e-4, "{gap}");
        }
    }

    /// A band the stencil can't have within the tolerance shrinks to one it can.
    #[test]
    fn the_band_shrinks_to_keep_the_error_within_the_tolerance() {
        let o = Stencil::optimised(8, 2.8, 1e-3).unwrap();
        let band = o.band().unwrap();
        eprintln!("order 8 within 1e-3: band {band:.3}");
        assert!(band < 2.8 && band > 1.8, "{band}");
        // Linearised in the design, exact here: allow a little.
        assert!(o.phase_error(band) < 1.1e-3, "{}", o.phase_error(band));
        // Within the tolerance already: the band stays.
        let o = Stencil::optimised(12, 1.8, 1e-3).unwrap();
        assert_eq!(o.band(), Some(1.8));
    }

    /// For a narrow band the best operator is nearly Taylor's, and the design stays sane.
    #[test]
    fn a_narrow_band_gives_nearly_taylor() {
        let t = Stencil::taylor(8).unwrap();
        let o = Stencil::optimised(8, 0.4, 1.0).unwrap();
        assert!(o.phase_error(0.4) <= t.phase_error(0.4) * 1.01);
        assert!(
            (o.sum_abs() / t.sum_abs() - 1.0).abs() < 0.01,
            "{}",
            o.sum_abs()
        );
    }

    #[test]
    fn follows_the_stride() {
        // Two interleaved lines; differentiating with stride 2 sees only one of them.
        let f: Vec<f32> = (0..40)
            .map(|i| if i % 2 == 0 { i as f32 } else { 1e6 })
            .collect();
        let mut out = [0.0f32; 1];
        for order in [4, 8] {
            let s = Stencil::taylor(order).unwrap();
            s.forward::<Separate>(&mut out, &f, 20, 2, 1.0);
            assert!((out[0] - 2.0).abs() < 1e-3, "order {order}");
            s.backward::<Separate>(&mut out, &f, 20, 2, 1.0);
            assert!((out[0] - 2.0).abs() < 1e-3, "order {order}");
        }
    }
}
