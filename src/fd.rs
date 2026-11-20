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

use anyhow::{Result, ensure};

/// The longest stencil offered: 16 points.
pub const MAX_ORDER: usize = 16;

/// A staggered first-derivative operator.
#[derive(Clone, Debug, PartialEq)]
pub struct Stencil {
    coefficients: Vec<f64>,
    /// The same, in single precision for the kernels.
    single: Vec<f32>,
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
        let single = coefficients.iter().map(|&c| c as f32).collect();
        Ok(Self {
            coefficients,
            single,
        })
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

    /// Σ |cₙ|, which sets the stability limit: dt ≤ 1 / (vp Σ|cₙ| √(Σ 1/h²)).
    pub fn sum_abs(&self) -> f64 {
        self.coefficients.iter().map(|c| c.abs()).sum()
    }

    /// The derivative half a node forward of each node, at p + s/2.
    #[inline]
    pub fn forward(&self, out: &mut [f32], f: &[f32], start: usize, stride: usize, inv_h: f32) {
        let len = out.len();
        for (k, &c) in self.single.iter().enumerate() {
            let a = c * inv_h;
            let hi = &f[start + (k + 1) * stride..start + (k + 1) * stride + len];
            let lo = &f[start - k * stride..start - k * stride + len];
            accumulate(out, a, hi, lo, k == 0);
        }
    }

    /// The derivative half a node backward of each node, at p - s/2.
    #[inline]
    pub fn backward(&self, out: &mut [f32], f: &[f32], start: usize, stride: usize, inv_h: f32) {
        let len = out.len();
        for (k, &c) in self.single.iter().enumerate() {
            let a = c * inv_h;
            let hi = &f[start + k * stride..start + k * stride + len];
            let lo = &f[start - (k + 1) * stride..start - (k + 1) * stride + len];
            accumulate(out, a, hi, lo, k == 0);
        }
    }
}

/// out = a (hi - lo), or out += a (hi - lo) after the first term.
#[inline]
fn accumulate(out: &mut [f32], a: f32, hi: &[f32], lo: &[f32], first: bool) {
    let n = out.len();
    let (hi, lo) = (&hi[..n], &lo[..n]);
    if first {
        for i in 0..n {
            out[i] = a * (hi[i] - lo[i]);
        }
    } else {
        for i in 0..n {
            out[i] += a * (hi[i] - lo[i]);
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

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
            s.forward(&mut out, &f, 10, 1, 1.0 / h);
            for (i, v) in out.iter().enumerate() {
                let x = (10 + i) as f32 * h + h / 2.0;
                assert!(
                    (v - df(x)).abs() < 1e-3,
                    "order {order} forward at {x}: {v}"
                );
            }
            s.backward(&mut out, &f, 10, 1, 1.0 / h);
            for (i, v) in out.iter().enumerate() {
                let x = (10 + i) as f32 * h - h / 2.0;
                assert!(
                    (v - df(x)).abs() < 1e-3,
                    "order {order} backward at {x}: {v}"
                );
            }
        }
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
            s.forward(&mut out, &f, 20, 2, 1.0);
            assert!((out[0] - 2.0).abs() < 1e-3, "order {order}");
            s.backward(&mut out, &f, 20, 2, 1.0);
            assert!((out[0] - 2.0).abs() < 1e-3, "order {order}");
        }
    }
}
