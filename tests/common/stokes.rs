//! Stokes's solution for a point force in a homogeneous full space, its spatial derivatives,
//! and the exact sensitivity kernels of the L2 misfit of one velocity trace: the reference
//! the gradient accuracy study measures the solver's gradients against.
//!
//! The displacement from a force along x_j with time function q (Aki & Richards 2002,
//! eq. 4.23), with r the distance, γ = x/r, a = r/α and b = r/β:
//!
//! u_i = 1/(4πρ) [ (3γᵢγⱼ - δᵢⱼ) N/r³ + γᵢγⱼ q(t - a)/(α² r) - (γᵢγⱼ - δᵢⱼ) q(t - b)/(β² r) ],
//!
//! N = ∫ₐᵇ τ q(t - τ) dτ = a Q₁(t - a) - b Q₁(t - b) + Q₂(t - a) - Q₂(t - b), where Q₁ and Q₂
//! are q's first and second integrals; [`stokes`] differentiates it in space exactly.
//!
//! For χ = ½ ∫₀ᵀ r(t)² dt, r the residual of component n at the receiver, a perturbation of
//! ρ, λ, μ changes χ by ∫ (K_ρ δρ + K_λ δλ + K_μ δμ) dV with (Tromp, Tape & Liu 2005;
//! Fichtner 2011, ch. 9; Born's approximation and reciprocity)
//!
//! K_ρ = -∫ ü·a dτ,  K_λ = -∫ (∇·a)(∇·u) dτ,  K_μ = -∫ ∂ⱼaᵢ (∂ⱼuᵢ + ∂ᵢuⱼ) dτ,
//!
//! u the forward displacement and a the adjoint field: at time τ, the velocity at time T - τ
//! of a force along x_n at the receiver whose time function is the residual reversed,
//! h(s) = r(T - s).

use std::f64::consts::PI;

use rayon::prelude::*;

use super::sub;

#[derive(Clone, Copy, Debug)]
pub struct Medium {
    pub vp: f64,
    pub vs: f64,
    pub rho: f64,
}

/// A time function q, as Stokes's solution needs it at time t: [q, q', Q₁, Q₂], the
/// integrals from -∞.
pub trait Pulse: Sync {
    fn eval(&self, t: f64) -> [f64; 4];
}

/// The `level`-th time derivative of `amplitude` times a Ricker wavelet (1 - 2a τ²) e^(-a τ²),
/// a = (π f0)², τ = t - t0, as the solver's `Wavelet::Ricker` defines it: level 0 gives the
/// displacement of a force with that wavelet, 1 the velocity, 2 the acceleration.
#[derive(Clone, Copy, Debug)]
pub struct Ricker {
    pub f0: f64,
    pub t0: f64,
    pub amplitude: f64,
    pub level: i32,
}

impl Ricker {
    /// The m-th derivative of the wavelet, m ≥ -2 meaning integrals: the wavelet is
    /// -E''/(2a) with E = e^(-a τ²), whose derivatives are Hermite polynomials,
    /// E⁽ⁿ⁾ = (-√a)ⁿ Hₙ(√a τ) E.
    pub fn derivative(&self, m: i32, t: f64) -> f64 {
        let a = (PI * self.f0).powi(2);
        let tau = t - self.t0;
        let n = (m + 2) as usize;
        let x = a.sqrt() * tau;
        let (mut h0, mut h1) = (1.0, 2.0 * x);
        let hermite = if n == 0 {
            1.0
        } else {
            for k in 1..n {
                let h2 = 2.0 * x * h1 - 2.0 * k as f64 * h0;
                h0 = h1;
                h1 = h2;
            }
            h1
        };
        -(-a.sqrt()).powi(n as i32) * hermite * (-a * tau * tau).exp() / (2.0 * a)
    }
}

impl Pulse for Ricker {
    fn eval(&self, t: f64) -> [f64; 4] {
        let k = self.level;
        [k, k + 1, k - 1, k - 2].map(|m| self.amplitude * self.derivative(m, t))
    }
}

/// The displacement u_i and its gradient ∂_l u_i (`[i][l]`) at offset `d` from a force along
/// x_j whose time function is `q`, at time t.
pub fn stokes(
    m: &Medium,
    j: usize,
    d: [f64; 3],
    t: f64,
    q: &impl Pulse,
) -> ([f64; 3], [[f64; 3]; 3]) {
    let r = (d[0] * d[0] + d[1] * d[1] + d[2] * d[2]).sqrt();
    let g = d.map(|x| x / r);
    let (al, be) = (m.vp, m.vs);
    let (a, b) = (r / al, r / be);
    let [p, dp, p1, p2] = q.eval(t - a);
    let [s, ds, s1, s2] = q.eval(t - b);
    let n = a * p1 - b * s1 + p2 - s2;
    let dn = r / (be * be) * s - r / (al * al) * p;
    let c = 1.0 / (4.0 * PI * m.rho);
    let delta = |i: usize, k: usize| if i == k { 1.0 } else { 0.0 };
    let (r2, r3, r4) = (r * r, r * r * r, r * r * r * r);
    let mut u = [0.0; 3];
    let mut grad = [[0.0; 3]; 3];
    for i in 0..3 {
        let gg = g[i] * g[j];
        let near = 3.0 * gg - delta(i, j);
        let far_s = gg - delta(i, j);
        u[i] = c * (near * n / r3 + gg * p / (al * al * r) - far_s * s / (be * be * r));
        for l in 0..3 {
            // ∂_l(γᵢγⱼ)
            let dgg = (delta(i, l) * g[j] + delta(j, l) * g[i] - 2.0 * gg * g[l]) / r;
            grad[i][l] = c
                * (3.0 * dgg * n / r3
                    + near * g[l] * (dn / r3 - 3.0 * n / r4)
                    + dgg * p / (al * al * r)
                    + gg / (al * al) * g[l] * (-dp / (al * r) - p / r2)
                    - dgg * s / (be * be * r)
                    - far_s / (be * be) * g[l] * (-ds / (be * r) - s / r2));
        }
    }
    (u, grad)
}

/// What the kernels are of: a force `amplitude`·w(t) along x_`force` at `source`, w a Ricker
/// wavelet, recorded as velocity component `component` at `receiver` for `duration`, against
/// the same recorded in `observed`.
#[derive(Clone, Copy, Debug)]
pub struct Problem {
    pub medium: Medium,
    pub observed: Medium,
    pub source: [f64; 3],
    pub force: usize,
    pub amplitude: f64,
    pub f0: f64,
    pub t0: f64,
    pub receiver: [f64; 3],
    pub component: usize,
    pub duration: f64,
}

impl Problem {
    fn ricker(&self, level: i32) -> Ricker {
        Ricker {
            f0: self.f0,
            t0: self.t0,
            amplitude: self.amplitude,
            level,
        }
    }

    /// The `k`-th time derivative of the recorded velocity in `medium` at time t.
    pub fn trace_derivative(&self, medium: &Medium, k: i32, t: f64) -> f64 {
        let d = sub(self.receiver, self.source);
        stokes(medium, self.force, d, t, &self.ricker(1 + k)).0[self.component]
    }

    /// The synthetic velocity at time t.
    pub fn synthetic(&self, t: f64) -> f64 {
        self.trace_derivative(&self.medium, 0, t)
    }

    /// The observed velocity at time t.
    pub fn observed(&self, t: f64) -> f64 {
        self.trace_derivative(&self.observed, 0, t)
    }

    /// The `k`-th derivative of the residual, synthetic less observed, at time t.
    pub fn residual(&self, k: i32, t: f64) -> f64 {
        self.trace_derivative(&self.medium, k, t) - self.trace_derivative(&self.observed, k, t)
    }

    /// The adjoint source's time function, the residual reversed, tabulated every `step`.
    pub fn adjoint_pulse(&self, step: f64) -> Table {
        let t = self.duration;
        let n = (t / step).ceil() as usize;
        let step = t / n as f64;
        let mut v: Vec<[f64; 5]> = (0..=n)
            .map(|i| {
                let s = i as f64 * step;
                // h⁽ᵏ⁾(s) = (-1)ᵏ r⁽ᵏ⁾(T - s)
                let r = |k: i32| self.residual(k, t - s) * if k % 2 == 0 { 1.0 } else { -1.0 };
                [r(0), r(1), r(2), r(3), 0.0]
            })
            .collect();
        // ∫₀ˢ h, interval by interval: the trapezium rule with its end correction, exact for
        // the cubic through the values and slopes.
        for i in 1..=n {
            let (a, b) = (v[i - 1], v[i]);
            v[i][4] = a[4] + step / 2.0 * (a[0] + b[0]) + step * step / 12.0 * (a[1] - b[1]);
        }
        Table { step, v }
    }

    /// The kernels' terms at each point, as the solver's staggered grid splits them:
    /// K_ρ by component (-∫ üᵢ aᵢ dτ for i = x, y, z, at the velocities' positions there),
    /// K_λ, K_μ's normal part (-∫ 2 Σᵢ ∂ᵢaᵢ ∂ᵢuᵢ dτ, at the nodes) and its shear parts for
    /// xy, xz, yz (-∫ (∂ⱼaᵢ + ∂ᵢaⱼ)(∂ⱼuᵢ + ∂ᵢuⱼ) dτ, at the shear stresses' positions). The time
    /// integrals by the trapezium rule every `dtau`. The points must keep clear of the source
    /// and the receiver, where the kernels are singular.
    pub fn terms(&self, points: &[[f64; 3]], dtau: f64) -> Vec<[f64; 8]> {
        let pulse = self.adjoint_pulse(1e-5);
        let t = self.duration;
        let n = (t / dtau).ceil() as usize;
        let dtau = t / n as f64;
        let (displacement, acceleration) = (self.ricker(0), self.ricker(2));
        points
            .par_iter()
            .map(|&x| {
                let (ds, dr) = (sub(x, self.source), sub(x, self.receiver));
                let mut k = [0.0; 8];
                for m in 0..=n {
                    let tau = m as f64 * dtau;
                    let w = if m == 0 || m == n { 0.5 } else { 1.0 } * dtau;
                    let (_, gu) = stokes(&self.medium, self.force, ds, tau, &displacement);
                    let (acc, _) = stokes(&self.medium, self.force, ds, tau, &acceleration);
                    let (a, ga) = stokes(&self.medium, self.component, dr, t - tau, &pulse);
                    for i in 0..3 {
                        k[i] -= w * acc[i] * a[i];
                    }
                    let div = |g: &[[f64; 3]; 3]| g[0][0] + g[1][1] + g[2][2];
                    k[3] -= w * div(&ga) * div(&gu);
                    k[4] -= w * 2.0 * (0..3).map(|i| ga[i][i] * gu[i][i]).sum::<f64>();
                    for (slot, (i, j)) in [(5, (0, 1)), (6, (0, 2)), (7, (1, 2))] {
                        k[slot] -= w * (ga[i][j] + ga[j][i]) * (gu[i][j] + gu[j][i]);
                    }
                }
                k
            })
            .collect()
    }

    /// K_ρ, K_λ, K_μ at each point (see [`Problem::terms`]).
    pub fn kernels(&self, points: &[[f64; 3]], dtau: f64) -> Vec<[f64; 3]> {
        self.terms(points, dtau)
            .into_iter()
            .map(|t| [t[0] + t[1] + t[2], t[3], t[4] + t[5] + t[6] + t[7]])
            .collect()
    }

    /// The exact kernels as the solver's grid sees a node's ρ, λ and μ: ρ through the
    /// arithmetic means at the six velocity positions around it (half to each), λ and μ's
    /// normal part at the node, μ's shear parts through the harmonic means at the twelve shear
    /// stress positions around it (a quarter to each in a homogeneous medium). The same
    /// correlations of exact fields, so what this differs from [`Problem::kernels`] by is the
    /// grid's sampling of the model, and what the solver's gradient differs from this by is
    /// the error of its wavefields.
    pub fn grid_kernels(&self, nodes: &[[f64; 3]], h: f64, dtau: f64) -> Vec<[f64; 3]> {
        let e = 0.5 * h;
        let mut offsets: Vec<[f64; 3]> = vec![[0.0; 3]];
        for a in 0..3 {
            for s in [-e, e] {
                let mut o = [0.0; 3];
                o[a] = s;
                offsets.push(o);
            }
        }
        for (i, j) in [(0, 1), (0, 2), (1, 2)] {
            for si in [-e, e] {
                for sj in [-e, e] {
                    let mut o = [0.0; 3];
                    o[i] = si;
                    o[j] = sj;
                    offsets.push(o);
                }
            }
        }
        let points: Vec<[f64; 3]> = nodes
            .iter()
            .flat_map(|x| {
                offsets
                    .iter()
                    .map(move |o| [x[0] + o[0], x[1] + o[1], x[2] + o[2]])
            })
            .collect();
        let terms = self.terms(&points, dtau);
        terms
            .chunks(offsets.len())
            .map(|t| {
                // t[0] at the node; t[1..7] along x-, x+, y-, y+, z-, z+; then four each for
                // the xy, xz and yz edges.
                let rho = (0..3)
                    .map(|a| 0.5 * (t[1 + 2 * a][a] + t[2 + 2 * a][a]))
                    .sum();
                let shear: f64 = (0..3)
                    .map(|p| (0..4).map(|q| 0.25 * t[7 + 4 * p + q][5 + p]).sum::<f64>())
                    .sum();
                [rho, t[0][3], t[0][4] + shear]
            })
            .collect()
    }
}

impl Problem {
    /// Σₙ δ(xₙ) kₙ h³ over `nodes` for the grid kernels kₙ of [`Problem::grid_kernels`],
    /// gathered by position rather than by node: each term once, at the node, the three
    /// velocity positions and the three shear positions after it, weighted by the means of
    /// `delta` over the nodes it averages. `delta` must vanish at the edges of `nodes`.
    pub fn grid_integral(
        &self,
        nodes: &[[f64; 3]],
        h: f64,
        delta: impl Fn([f64; 3]) -> f64 + Sync,
        dtau: f64,
    ) -> [f64; 3] {
        let e = 0.5 * h;
        let shift = |x: [f64; 3], o: [f64; 3]| [x[0] + o[0], x[1] + o[1], x[2] + o[2]];
        let unit = |a: usize, s: f64| {
            let mut o = [0.0; 3];
            o[a] = s;
            o
        };
        // Per node: the node, vx, vy, vz positions, then the xy, xz, yz edges.
        let pairs = [(0, 1), (0, 2), (1, 2)];
        let mut points = Vec::with_capacity(7 * nodes.len());
        let mut weights = Vec::with_capacity(7 * nodes.len());
        for &x in nodes {
            points.push(x);
            weights.push(delta(x));
            for a in 0..3 {
                points.push(shift(x, unit(a, e)));
                weights.push(0.5 * (delta(x) + delta(shift(x, unit(a, h)))));
            }
            for (i, j) in pairs {
                let mut w = 0.0;
                for (si, sj) in [(0.0, 0.0), (h, 0.0), (0.0, h), (h, h)] {
                    w += 0.25 * delta(shift(shift(x, unit(i, si)), unit(j, sj)));
                }
                points.push(shift(shift(x, unit(i, e)), unit(j, e)));
                weights.push(w);
            }
        }
        let terms = self.terms(&points, dtau);
        let mut k = [0.0; 3];
        for (t, w) in terms.chunks(7).zip(weights.chunks(7)) {
            k[1] += w[0] * t[0][3];
            k[2] += w[0] * t[0][4];
            for a in 0..3 {
                k[0] += w[1 + a] * t[1 + a][a];
                k[2] += w[4 + a] * t[4 + a][5 + a];
            }
        }
        k.map(|x| x * h * h * h)
    }
}

/// The kernels for vp, vs (ρ fixed) and ρ (vp and vs fixed), from those for ρ, λ, μ, as
/// `Gradient::velocities` converts the solver's.
pub fn velocity_kernels(m: &Medium, k: [f64; 3]) -> [f64; 3] {
    let [kr, kl, km] = k;
    let (a, b, r) = (m.vp, m.vs, m.rho);
    [
        2.0 * r * a * kl,
        2.0 * r * b * (km - 2.0 * kl),
        kr + (a * a - 2.0 * b * b) * kl + b * b * km,
    ]
}

/// A time function tabulated from 0 every `step`, with its first three derivatives and its
/// integral (`v[i]` = [h, h', h'', h''', ∫₀ h]), zero outside; as a pulse it gives the
/// velocity of a force with time function h: [h', h'', h, ∫h], each interpolated by the
/// cubic through its values and slopes.
pub struct Table {
    pub step: f64,
    pub v: Vec<[f64; 5]>,
}

impl Pulse for Table {
    fn eval(&self, t: f64) -> [f64; 4] {
        let last = self.v.len() - 1;
        if t <= 0.0 {
            return [0.0; 4];
        }
        let x = t / self.step;
        if x >= last as f64 {
            return [0.0, 0.0, 0.0, self.v[last][4]];
        }
        let i = x.floor() as usize;
        let s = x - i as f64;
        let (a, b) = (self.v[i], self.v[i + 1]);
        let h = self.step;
        // Cubic Hermite: values y and slopes d at both ends.
        let hermite = |y0: f64, d0: f64, y1: f64, d1: f64| {
            let s2 = s * s;
            let s3 = s2 * s;
            (2.0 * s3 - 3.0 * s2 + 1.0) * y0
                + (s3 - 2.0 * s2 + s) * h * d0
                + (-2.0 * s3 + 3.0 * s2) * y1
                + (s3 - s2) * h * d1
        };
        [
            hermite(a[1], a[2], b[1], b[2]),
            hermite(a[2], a[3], b[2], b[3]),
            hermite(a[0], a[1], b[0], b[1]),
            hermite(a[4], a[0], b[4], b[0]),
        ]
    }
}
