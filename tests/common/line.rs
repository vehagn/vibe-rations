//! Closed-form solutions for line sources in a homogeneous full space: the 2D (plane-strain)
//! counterparts of the explosion and Stokes's point force in `mod.rs`, for 2D runs.
//!
//! In 2D the scalar Green's function of (∇² - ∂²ₜ/c²) g = -δ(x)δ(t) is
//! g_c = H(t - r/c) / (2π √(t² - r²/c²)), whose tail never ends. With t = (r/c) cosh u it
//! becomes du/2π, so every solution below is an integral over u of the source's time function
//! at t - (r/c) cosh u, with smooth weights, taken by Simpson's rule. The elastic Green's tensor
//! has the same form in any dimension (Aki & Richards 2002, §4.2, written with g_c):
//! G_ij = (1/ρ) [δ_ij g_β/β² + ∂_i ∂_j (I_α - I_β)], I_c the double time integral of g_c, whose
//! second derivatives in 2D are (1/(2πr²)) [γ_i γ_j c t²/Q - (δ_ij - γ_i γ_j) Q/c],
//! Q = √(c²t² - r²). `cross_checks` in `tests/physics.rs` tests them against the 3D solutions
//! integrated along the line.

use std::f64::consts::PI;

use vibe_rations::source::Wavelet;

use super::{RHO, VP, VS};

/// ∫₀^U f(u) du for U = acosh(c t / r), the u that reaches back to the source's start, by
/// Simpson's rule; zero before the wave arrives.
fn tail(c: f64, r: f64, t: f64, f: impl Fn(f64, f64) -> f64) -> f64 {
    let x = c * t / r;
    if x <= 1.0 {
        return 0.0;
    }
    let upper = x.acosh();
    // Fine enough that a step in u is a small part of the wavelet's period even late in the
    // tail, where t - (r/c) cosh u moves by about t per unit u.
    let n = 2 * ((upper * 2000.0).ceil() as usize).max(200);
    let du = upper / n as f64;
    let at = |m: usize| {
        let u = m as f64 * du;
        f(u, t - r / c * u.cosh())
    };
    let inner: f64 = (1..n)
        .map(|m| if m % 2 == 1 { 4.0 } else { 2.0 } * at(m))
        .sum();
    (at(0) + inner + at(n)) * du / 3.0
}

/// Pressure and radial particle velocity at distance r from a line explosion whose moment rate
/// per metre is m0·w(t), in the medium of `mod.rs`: with φ = -(M/ρα²) * g_α,
/// p = -K ∇·u = (K/ρα⁴) m0 w' * g_α and v_r = -(m0/ρα²) ∂_r (w * g_α).
pub fn line_explosion(w: &Wavelet, m0: f64, r: f64, t: f64) -> (f64, f64) {
    let k = RHO * (VP * VP - 4.0 / 3.0 * VS * VS);
    let p = k * m0 / (2.0 * PI * RHO * VP.powi(4)) * tail(VP, r, t, |_, s| w.derivative(s));
    let v = m0 / (2.0 * PI * RHO * VP.powi(3)) * tail(VP, r, t, |u, s| u.cosh() * w.derivative(s));
    (p, v)
}

/// Particle velocity (vx, vz) at offset `d` = (dx, dz) from a line force of f·w(t) per metre
/// along `along` (0 for x, 1 for z):
/// v_i = (f/2πρ) [(1/α²) ∫ (γ_i γ_j cosh 2u - δ_ij sinh²u) w'(t - (r/α) cosh u) du
///              + (1/β²) ∫ (δ_ij cosh²u - γ_i γ_j cosh 2u) w'(t - (r/β) cosh u) du].
pub fn line_force_velocity(w: &Wavelet, f: f64, along: usize, d: [f64; 2], t: f64) -> [f64; 2] {
    let r = d[0].hypot(d[1]);
    let g = d.map(|x| x / r);
    std::array::from_fn(|i| {
        let gg = g[i] * g[along];
        let delta = if i == along { 1.0 } else { 0.0 };
        let p = tail(VP, r, t, |u, s| {
            (gg * (2.0 * u).cosh() - delta * u.sinh().powi(2)) * w.derivative(s)
        });
        let s = tail(VS, r, t, |u, s| {
            (delta * u.cosh().powi(2) - gg * (2.0 * u).cosh()) * w.derivative(s)
        });
        f / (2.0 * PI * RHO) * (p / (VP * VP) + s / (VS * VS))
    })
}

/// ∫ f(y) dy along the line through the point at distance r from the line's foot, over the
/// part of it whose waves (at speed `c`, the fastest) reach the point by time t: with
/// y = r sinh s, dy = √(r² + y²) ds, by Simpson's rule. What a 3D point source's solution,
/// summed over a line of them, gives in 2D.
pub fn along_line(c: f64, r: f64, t: f64, f: impl Fn(f64) -> f64) -> f64 {
    let x = c * t / r;
    if x <= 1.0 {
        return 0.0;
    }
    let upper = x.acosh();
    let n = 2 * ((upper * 600.0).ceil() as usize).max(200);
    let ds = upper / n as f64;
    let at = |m: usize| {
        let s = m as f64 * ds;
        let y = r * s.sinh();
        // Both halves of the line.
        r * s.cosh() * (f(y) + f(-y))
    };
    let inner: f64 = (1..n)
        .map(|m| if m % 2 == 1 { 4.0 } else { 2.0 } * at(m))
        .sum();
    (at(0) + inner + at(n)) * ds / 3.0
}
