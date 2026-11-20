//! Staggered-grid first derivatives along one axis of a padded array, a row of nodes at a time.
//!
//! `start` is the padded index of the row's first node and `stride` the distance in memory
//! between neighbours along the axis being differentiated. Every function fills `out` with
//! the derivative for `out.len()` consecutive nodes (x fastest), each placed half a node
//! forward or backward of the node it is computed for.

#![allow(clippy::needless_range_loop)]

/// Fourth-order staggered coefficients (Levander 1988).
pub const C1: f32 = 9.0 / 8.0;
pub const C2: f32 = -1.0 / 24.0;

/// |C1| + |C2|, which sets the stability limit of the scheme.
pub const STENCIL_SUM: f64 = 9.0 / 8.0 + 1.0 / 24.0;

/// Fourth-order derivative half a node forward: at p + s/2 from f at p - s, p, p + s, p + 2s.
#[inline]
pub fn forward(out: &mut [f32], f: &[f32], start: usize, stride: usize, inv_h: f32) {
    let n = out.len();
    let m1 = &f[start - stride..start - stride + n];
    let p0 = &f[start..start + n];
    let p1 = &f[start + stride..start + stride + n];
    let p2 = &f[start + 2 * stride..start + 2 * stride + n];
    let (a, b) = (C1 * inv_h, C2 * inv_h);
    for i in 0..n {
        out[i] = a * (p1[i] - p0[i]) + b * (p2[i] - m1[i]);
    }
}

/// Fourth-order derivative half a node backward: at p - s/2 from f at p - 2s, p - s, p, p + s.
#[inline]
pub fn backward(out: &mut [f32], f: &[f32], start: usize, stride: usize, inv_h: f32) {
    let n = out.len();
    let m2 = &f[start - 2 * stride..start - 2 * stride + n];
    let m1 = &f[start - stride..start - stride + n];
    let p0 = &f[start..start + n];
    let p1 = &f[start + stride..start + stride + n];
    let (a, b) = (C1 * inv_h, C2 * inv_h);
    for i in 0..n {
        out[i] = a * (p0[i] - m1[i]) + b * (p1[i] - m2[i]);
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// A cubic sampled on a 1D line: the fourth-order stencils differentiate it exactly.
    #[test]
    fn differentiates_a_cubic_exactly() {
        let h = 0.5f32;
        let f: Vec<f32> = (0..12)
            .map(|i| {
                let x = i as f32 * h;
                x * x * x - 2.0 * x
            })
            .collect();
        let df = |x: f32| 3.0 * x * x - 2.0;
        let mut out = [0.0f32; 4];
        forward(&mut out, &f, 3, 1, 1.0 / h);
        for (i, v) in out.iter().enumerate() {
            let x = (3 + i) as f32 * h + h / 2.0;
            assert!((v - df(x)).abs() < 1e-4, "forward at {x}: {v} vs {}", df(x));
        }
        backward(&mut out, &f, 3, 1, 1.0 / h);
        for (i, v) in out.iter().enumerate() {
            let x = (3 + i) as f32 * h - h / 2.0;
            assert!(
                (v - df(x)).abs() < 1e-4,
                "backward at {x}: {v} vs {}",
                df(x)
            );
        }
    }

    #[test]
    fn follows_the_stride() {
        // Two interleaved lines; differentiating with stride 2 sees only one of them.
        let f: Vec<f32> = (0..20)
            .map(|i| if i % 2 == 0 { i as f32 } else { 1e6 })
            .collect();
        let mut out = [0.0f32; 1];
        forward(&mut out, &f, 6, 2, 1.0);
        assert!((out[0] - 2.0).abs() < 1e-4);
        backward(&mut out, &f, 6, 2, 1.0);
        assert!((out[0] - 2.0).abs() < 1e-4);
    }
}
