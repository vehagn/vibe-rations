//! Removing leapfrog's time error: the time-dispersion transforms of Wang & Xu (2015) and
//! Koene et al. (2018).
//!
//! Each step replaces the time derivatives by centred differences over dt, whose symbol at
//! frequency ω is iω′ with ω′ = (2/dt) sin(ω dt/2). Nothing else in a step depends on time
//! (the medium, the stencil, the free surface and where sources enter are the same every step),
//! so the stepped wavefield at ω is the time-continuous ("semi-discrete") one at ω′. Leapfrog's
//! whole time error is that warp of the frequency axis, the same for every wave speed, stencil
//! and position, and two transforms undo it:
//!
//! - before the run, a source's samples get the spectrum S(ω′(ω)) instead of S(ω)
//!   ([`Warp::to_leapfrog`]), so the wavefield at ω is the semi-discrete response to S at ω′;
//! - after it, a trace's spectrum at ω′ is read from the recorded one at
//!   ω = (2/dt) arcsin(ω′ dt/2) ([`Warp::from_leapfrog`]).
//!
//! What is left is the spatial error alone. That holds exactly in the interior and at the free
//! surface; the CPML's memory variables and the sponge are recursive filters of their own,
//! which the warp maps only approximately, and that matters only to waves being absorbed. The
//! wavefield itself (snapshots) stays warped. Time stamps matter: the warp is about t = 0, when
//! the run starts at rest, so each signal is given with the time of its first sample.
//!
//! The pair is Koene et al.'s (2018, eqs. 13–17), with no amplitude factors; they pad to
//! twice the length too. The time stamps of a staggered scheme are this module's concern:
//! the paper doesn't discuss the time origin.
//!
//! Each transform is a real linear map from one sampled signal to another of the same length,
//! band-limited to the frequencies kept and computed as if zero-padded to twice its length, so
//! nothing wraps round. [`Warp::apply_transpose`] is what an adjoint source needs when the
//! misfit is measured on transformed traces.

use std::f64::consts::PI;

/// A linear map from samples x(t₀ + nΔ) to samples y(t₀ + mΔ), n, m < len: y's spectrum at
/// each kept frequency ω is x's at some other frequency, and zero above the band.
#[derive(Clone, Debug)]
pub struct Warp {
    t0: f64,
    interval: f64,
    len: usize,
    /// For each frequency kept, in rad/s: (ω of the output, ω of the input it comes from,
    /// weight of the bin).
    bins: Vec<(f64, f64, f64)>,
}

/// How often the rotating phasors are recomputed from scratch, so round-off can't accumulate.
const RESEED: usize = 1024;

impl Warp {
    /// For a source stepped with time step `dt`: samples whose spectrum at ω is the input's at
    /// ω′ = (2/dt) sin(ω dt/2), up to `f_max` Hz.
    pub fn to_leapfrog(dt: f64, t0: f64, interval: f64, len: usize, f_max: f64) -> Self {
        Self::new(t0, interval, len, f_max, |w| {
            (w * dt <= PI).then(|| 2.0 / dt * (0.5 * w * dt).sin())
        })
    }

    /// For a trace recorded from a run with time step `dt`: samples whose spectrum at ω′ is
    /// the input's at ω = (2/dt) arcsin(ω′ dt/2), up to `f_max` Hz.
    pub fn from_leapfrog(dt: f64, t0: f64, interval: f64, len: usize, f_max: f64) -> Self {
        Self::new(t0, interval, len, f_max, |w| {
            (0.5 * w * dt < 1.0).then(|| 2.0 / dt * (0.5 * w * dt).asin())
        })
    }

    /// The band limit alone, without a warp: to compare like with like.
    pub fn band_limit(t0: f64, interval: f64, len: usize, f_max: f64) -> Self {
        Self::new(t0, interval, len, f_max, Some)
    }

    /// Keeps the frequencies 2πj/P up to `f_max` (and below the samples' Nyquist frequency)
    /// whose source frequency `map` gives, P being twice the signal's length.
    fn new(
        t0: f64,
        interval: f64,
        len: usize,
        f_max: f64,
        map: impl Fn(f64) -> Option<f64>,
    ) -> Self {
        assert!(interval > 0.0 && len > 0, "a warp needs samples");
        let period = 2.0 * len as f64 * interval;
        let nyquist = PI / interval;
        let top = (2.0 * PI * f_max).min(nyquist);
        let mut bins = Vec::new();
        for j in 0.. {
            let out = 2.0 * PI * j as f64 / period;
            if out > top {
                break;
            }
            match map(out) {
                Some(input) if input <= nyquist => {
                    // The inverse DFT over ±ω: the bins above zero count twice.
                    let weight = if j == 0 { 1.0 } else { 2.0 } * interval / period;
                    bins.push((out, input, weight));
                }
                _ => break,
            }
        }
        Self {
            t0,
            interval,
            len,
            bins,
        }
    }

    pub fn len(&self) -> usize {
        self.len
    }

    pub fn is_empty(&self) -> bool {
        self.len == 0
    }

    /// The highest frequency kept, in Hz.
    pub fn f_max(&self) -> f64 {
        self.bins.last().map_or(0.0, |b| b.0 / (2.0 * PI))
    }

    /// y = W x.
    pub fn apply(&self, x: &[f64]) -> Vec<f64> {
        assert_eq!(x.len(), self.len, "the signal's length");
        let mut y = vec![0.0; self.len];
        for &(out, input, weight) in &self.bins {
            let (re, im) = self.transform(x, -input);
            self.add_wave(&mut y, out, weight * re, weight * im);
        }
        y
    }

    /// x = Wᵀ y.
    pub fn apply_transpose(&self, y: &[f64]) -> Vec<f64> {
        assert_eq!(y.len(), self.len, "the signal's length");
        let mut x = vec![0.0; self.len];
        for &(out, input, weight) in &self.bins {
            let (re, im) = self.transform(y, out);
            self.add_wave(&mut x, -input, weight * re, weight * im);
        }
        x
    }

    /// Σₙ x[n] e^{iωtₙ}.
    fn transform(&self, x: &[f64], w: f64) -> (f64, f64) {
        let (mut re, mut im) = (0.0, 0.0);
        for (k, chunk) in x.chunks(RESEED).enumerate() {
            let mut p = self.phasor(w, k * RESEED);
            let r = Phasor::of(w * self.interval);
            for &v in chunk {
                re += v * p.0;
                im += v * p.1;
                p = r.turn(p);
            }
        }
        (re, im)
    }

    /// out[m] += Re((re + i·im) e^{iωtₘ}).
    fn add_wave(&self, out: &mut [f64], w: f64, re: f64, im: f64) {
        for (k, chunk) in out.chunks_mut(RESEED).enumerate() {
            let mut p = self.phasor(w, k * RESEED);
            let r = Phasor::of(w * self.interval);
            for o in chunk {
                *o += re * p.0 - im * p.1;
                p = r.turn(p);
            }
        }
    }

    /// e^{iωtₙ}.
    fn phasor(&self, w: f64, n: usize) -> Phasor {
        Phasor::of(w * (self.t0 + n as f64 * self.interval))
    }
}

/// A unit complex number, (cos, sin).
#[derive(Clone, Copy)]
struct Phasor(f64, f64);

impl Phasor {
    /// e^{i·angle}.
    fn of(angle: f64) -> Phasor {
        let (s, c) = angle.sin_cos();
        Phasor(c, s)
    }

    /// self · p.
    fn turn(self, p: Phasor) -> Phasor {
        Phasor(self.0 * p.0 - self.1 * p.1, self.0 * p.1 + self.1 * p.0)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::source::Wavelet;

    fn ricker(f0: f64, t0: f64, times: impl Iterator<Item = f64>) -> Vec<f64> {
        let w = Wavelet::Ricker { f0, t0 };
        times.map(|t| w.value(t)).collect()
    }

    fn max_abs_diff(a: &[f64], b: &[f64]) -> f64 {
        a.iter()
            .zip(b)
            .map(|(x, y)| (x - y).abs())
            .fold(0.0, f64::max)
    }

    #[test]
    fn the_band_limit_keeps_a_signal_inside_the_band() {
        let (dt, n) = (2e-3, 600);
        let x = ricker(10.0, 0.15, (0..n).map(|k| 0.5 * dt + k as f64 * dt));
        let w = Warp::band_limit(0.5 * dt, dt, n, 50.0);
        assert!((w.f_max() - 50.0).abs() < 1.0);
        assert!(max_abs_diff(&w.apply(&x), &x) < 1e-8);
    }

    #[test]
    fn the_transpose_is_the_transpose() {
        let (dt, n) = (1e-3, 3000);
        let w = Warp::from_leapfrog(dt, dt, 2.0 * dt, n, 60.0);
        // Two arbitrary signals: <W x, y> = <x, Wᵀ y>.
        let x: Vec<f64> = (0..n).map(|k| ((k * 7919) % 101) as f64 - 50.0).collect();
        let y: Vec<f64> = (0..n).map(|k| ((k * 104_729) % 89) as f64 - 44.0).collect();
        let dot = |a: &[f64], b: &[f64]| a.iter().zip(b).map(|(p, q)| p * q).sum::<f64>();
        let (lhs, rhs) = (dot(&w.apply(&x), &y), dot(&x, &w.apply_transpose(&y)));
        assert!((lhs - rhs).abs() < 1e-10 * lhs.abs(), "{lhs} vs {rhs}");
    }

    #[test]
    fn un_warping_undoes_warping() {
        // A coarse step: 20 samples per dominant period. The band, five times the dominant
        // frequency, holds the whole wavelet, which the pair must give back.
        let (dt, n) = (5e-3, 400);
        let x = ricker(10.0, 0.15, (0..n).map(|k| k as f64 * dt));
        let to = Warp::to_leapfrog(dt, 0.0, dt, n, 50.0);
        let from = Warp::from_leapfrog(dt, 0.0, dt, n, 50.0);
        let warped = to.apply(&x);
        assert!(
            max_abs_diff(&warped, &x) > 0.01,
            "the warp should matter here"
        );
        // Measured: 8e-9.
        assert!(max_abs_diff(&from.apply(&warped), &x) < 1e-7);
    }

    /// A 1D staggered velocity–stress run, second order in space, leapfrog in time, timed like
    /// the solver: a force enters the velocity at whole steps, and the velocity is read half a
    /// step later. Returns the receiver's trace.
    fn run_1d(courant: f64, steps: usize, source: impl Fn(usize) -> f64) -> Vec<f64> {
        let (n, at, to) = (1600, 200, 500);
        let mut v = vec![0.0f64; n];
        let mut s = vec![0.0f64; n];
        let mut trace = Vec::with_capacity(steps);
        for step in 0..steps {
            for i in 1..n {
                v[i] += courant * (s[i] - s[i - 1]);
            }
            v[at] += courant * source(step);
            for i in 0..n - 1 {
                s[i] += courant * (v[i + 1] - v[i]);
            }
            trace.push(v[to]);
        }
        trace
    }

    /// With the source warped and the trace un-warped, a coarse time step gives the trace of a
    /// fine one: what is left is the spatial error, which both share.
    #[test]
    fn removes_the_time_error_of_a_staggered_run() {
        // Unit spacing and speed; 300 nodes to the receiver, about ten dominant wavelengths.
        let (f0, t0, f_max) = (0.035, 40.0, 0.0875 * 2.0);
        let run = |courant: f64, duration: f64, transform: bool| {
            let dt = courant;
            let steps = (duration / dt).round() as usize;
            let mut src = ricker(f0, t0, (0..steps).map(|k| k as f64 * dt));
            if transform {
                src = Warp::to_leapfrog(dt, 0.0, dt, steps, f_max).apply(&src);
            }
            let trace = run_1d(courant, steps, |k| src[k]);
            if transform {
                Warp::from_leapfrog(dt, 0.5 * dt, dt, steps, f_max).apply(&trace)
            } else {
                trace
            }
        };
        // Coarse: C = 0.9. Fine: C = 0.9/9, so every ninth fine sample, offset by four, falls
        // on a coarse one: (9k + 4 + ½)·dt/9 = (k + ½)·dt.
        let duration = 0.9 * 600.0;
        let fine = run(0.1, duration, true);
        let fine: Vec<f64> = fine.iter().skip(4).step_by(9).copied().collect();
        let plain = run(0.9, duration, false);
        let fixed = run(0.9, duration, true);
        let peak = fine.iter().fold(0.0f64, |m, v| m.max(v.abs()));
        let err = |a: &[f64]| max_abs_diff(a, &fine) / peak;
        // Measured: 24% of the peak without the transform, 1.1e-8 with it.
        assert!(
            err(&plain) > 0.1,
            "leapfrog's error at C = 0.9: {}",
            err(&plain)
        );
        assert!(err(&fixed) < 1e-6, "after the transform: {}", err(&fixed));
    }
}
