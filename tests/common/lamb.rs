//! Lamb's problem for a Poisson solid: the vertical motion of the free surface of a homogeneous
//! half-space with λ = μ (vp = √3 vs) under a vertical point force on the surface, after
//! Pekeris (1955).
//!
//! For a downward force F H(t), with τ = β t / r (β the S speed, r the distance) and
//! γ² = (3 + √3)/4 the Rayleigh pole (his eq. 27), the downward displacement is (his eqs.
//! 45–48, with his Z, the surface integral of the applied stress, equal to -F):
//!
//! w = 0, τ < 1/√3;
//! w = F/(32πμr) [6 - √3/√(τ² - ¼) - √(3√3 + 5)/√(γ² - τ²) + √(3√3 - 5)/√(τ² - (3 - √3)/4)],
//!     1/√3 < τ < 1;
//! w = F/(16πμr) [6 - √(3√3 + 5)/√(γ² - τ²)], 1 < τ < γ;
//! w = 3F/(8πμr), τ > γ: Boussinesq's static displacement for ν = ¼.
//!
//! The surface first rises at the P arrival, rises without bound as the Rayleigh wave
//! arrives (an integrable singularity) and then sits at the static value.
//!
//! For a force F s(t) the displacement is w * s' and the velocity w * s''. Split at the
//! arrivals, v(t) = ∫ w(t') s''(t - t') dt' over (t_P, t_R) plus w_∞ s'(t - t_R): after t_R
//! the response is the static value, whose integral against s'' is exact. On (t_P, t_S) the
//! substitution t' = t_P + u² and on (t_S, t_R) t' = t_R - u² make the integrand smooth (the
//! second removes the singularity), and Simpson's rule takes it.

use std::f64::consts::PI;

use super::stokes::Ricker;

/// A vertical force `force` (N, downward) times a Ricker wavelet on the surface of a Poisson
/// solid with S speed `vs` and density `rho`.
#[derive(Clone, Copy, Debug)]
pub struct Lamb {
    pub vs: f64,
    pub rho: f64,
    pub force: f64,
    pub wavelet: Ricker,
}

/// The Rayleigh pole: τ = γ is the Rayleigh wave's arrival, its speed vs/γ = 0.9194 vs.
pub fn gamma() -> f64 {
    0.5 * (3.0 + 3f64.sqrt()).sqrt()
}

impl Lamb {
    fn mu(&self) -> f64 {
        self.rho * self.vs * self.vs
    }

    /// Pekeris's downward displacement at distance `r` and time `t` for a downward force F H(t).
    pub fn step(&self, r: f64, t: f64) -> f64 {
        let s3 = 3f64.sqrt();
        let g2 = gamma() * gamma();
        let tau = self.vs * t / r;
        let a = self.force / (PI * self.mu() * r);
        if tau <= 1.0 / s3 {
            0.0
        } else if tau < 1.0 {
            let t2 = tau * tau;
            a / 32.0
                * (6.0 - s3 / (t2 - 0.25).sqrt() - (3.0 * s3 + 5.0).sqrt() / (g2 - t2).sqrt()
                    + (3.0 * s3 - 5.0).sqrt() / (t2 - (3.0 - s3) / 4.0).sqrt())
        } else if tau < gamma() {
            a / 16.0 * (6.0 - (3.0 * s3 + 5.0).sqrt() / (g2 - tau * tau).sqrt())
        } else {
            self.static_displacement(r)
        }
    }

    /// Boussinesq's static displacement, 3F/(8πμr) for ν = ¼.
    pub fn static_displacement(&self, r: f64) -> f64 {
        3.0 * self.force / (8.0 * PI * self.mu() * r)
    }

    /// The downward velocity at distance `r` and time `t`, with `n` Simpson intervals for
    /// each of the two stretches between the arrivals.
    pub fn velocity_with(&self, r: f64, t: f64, n: usize) -> f64 {
        let s = |m: i32, t: f64| self.wavelet.derivative(m, t);
        let (tp, ts) = (r / (3f64.sqrt() * self.vs), r / self.vs);
        let tr = gamma() * r / self.vs;
        // Simpson's rule over u in (0, len) of f(u).
        let simpson = |len: f64, f: &dyn Fn(f64) -> f64| {
            let h = len / n as f64;
            let inner: f64 = (1..n)
                .map(|k| f(k as f64 * h) * if k % 2 == 1 { 4.0 } else { 2.0 })
                .sum();
            (f(0.0) + inner + f(len)) * h / 3.0
        };
        // t' = tp + u²: dt' = 2u du.
        let early = simpson((ts - tp).sqrt(), &|u| {
            let tt = tp + u * u;
            self.step(r, tt) * s(2, t - tt) * 2.0 * u
        });
        // t' = tr - u²: the singular term ~ 1/u, times 2u; at u = 0 its limit.
        let late = simpson((tr - ts).sqrt(), &|u| {
            if u == 0.0 {
                return self.singular_limit(r) * s(2, t - tr);
            }
            let tt = tr - u * u;
            self.step(r, tt) * s(2, t - tt) * 2.0 * u
        });
        early + late + self.static_displacement(r) * s(1, t - tr)
    }

    /// The velocity, with enough intervals that doubling them changes it by less than 1e-9
    /// of its peak (checked in the tests).
    pub fn velocity(&self, r: f64, t: f64) -> f64 {
        self.velocity_with(r, t, 2000)
    }

    /// lim u→0 of 2u w(t_R - u²), from the singular term -F/(16πμr) √(3√3 + 5)/√(γ² - τ²):
    /// γ² - τ² = (vs/r)² (t_R² - t'²) ≈ (vs/r)² 2 t_R u² = (2γ vs/r) u², with t_R = γ r/vs.
    fn singular_limit(&self, r: f64) -> f64 {
        let s3 = 3f64.sqrt();
        let a = self.force / (PI * self.mu() * r);
        let k = (2.0 * gamma() * self.vs / r).sqrt();
        -2.0 * a / 16.0 * (3.0 * s3 + 5.0).sqrt() / k
    }
}

/// The survey the tests share: a Poisson solid (vp 3000 m/s, vs 3000/√3, ρ 2500 kg/m³) with
/// a free surface, a downward force of 1e9 N with a 15 Hz Ricker wavelet on the surface, and
/// vz recorded on the surface at 150, 240 and 330 m along x and 254.6 m along the diagonal.
/// Every position is a multiple of 30 m, so source and receivers sit on nodes at spacings of
/// 10, 7.5, 5 and 3.75 m; the absorbing layers are 120 m thick, and the receivers 100 m clear
/// of them. The Rayleigh wave is 106 m long at f0, 10.6 points at h = 10 m, and the run's
/// shortest S wavelength (at 2.5 f0) 4.6 points.
pub mod survey {
    use super::{Lamb, Ricker};
    use crate::common::{misfit, norm, sub, trace};
    use vibe_rations::receiver::Component;
    use vibe_rations::simulation::Simulation;

    pub const VP: f64 = 3000.0;
    pub const RHO: f64 = 2500.0;
    pub const F0: f64 = 15.0;
    pub const SOURCE: [f64; 3] = [210.0, 240.0, 0.0];
    pub const RECEIVERS: [[f64; 3]; 4] = [
        [360.0, 240.0, 0.0],
        [450.0, 240.0, 0.0],
        [540.0, 240.0, 0.0],
        [390.0, 420.0, 0.0],
    ];

    pub fn vs() -> f64 {
        VP / 3f64.sqrt()
    }

    /// The configuration at spacing `h`, with `scheme` and `time` added to those sections.
    pub fn config(h: f64, scheme: &str, time: &str) -> String {
        let n = |metres: f64| (metres / h).round() as usize + 1;
        let receivers: String = RECEIVERS
            .iter()
            .map(|p| format!("[[receivers.point]]\nposition = {p:?}\n"))
            .collect();
        format!(
            r#"
            [grid]
            points = [{}, {}, {}]
            spacing = {h}
            [time]
            duration = 0.42
            {time}
            [boundary]
            pml_width = {}
            free_surface = true
            [scheme]
            {scheme}
            [[model.layer]]
            top = 0.0
            vp = {VP}
            vs = {}
            rho = {RHO}
            [[source]]
            position = {SOURCE:?}
            force = [0.0, 0.0, 1e9]
            wavelet = {{ kind = "ricker", frequency = {F0} }}
            [receivers]
            components = ["vz"]
            {receivers}
            "#,
            n(750.0),
            n(630.0),
            n(360.0),
            (120.0 / h).round() as usize,
            vs(),
        )
    }

    /// The exact answer: Pekeris's for this survey's force and wavelet (t0 = 1.2/f0, as the
    /// solver's Ricker wavelet has it).
    pub fn lamb() -> Lamb {
        Lamb {
            vs: vs(),
            rho: RHO,
            force: 1e9,
            wavelet: Ricker {
                f0: F0,
                t0: 1.2 / F0,
                amplitude: 1.0,
                level: 0,
            },
        }
    }

    /// Each receiver's distance; its vz trace's misfit against Pekeris's; its amplitude
    /// relative to his, ⟨v, w⟩/⟨w, w⟩; how late its largest swing comes, in samples (a
    /// parabola through the peak); and that swing's size relative to his.
    pub fn compare(sim: &Simulation) -> Vec<[f64; 5]> {
        let lamb = lamb();
        let times = sim.sample_times();
        let peak_at = |v: &[f64]| {
            let k = (1..v.len() - 1)
                .max_by(|&a, &b| v[a].abs().total_cmp(&v[b].abs()))
                .unwrap();
            let (a, b, c) = (v[k - 1].abs(), v[k].abs(), v[k + 1].abs());
            let d = 0.5 * (a - c) / (a - 2.0 * b + c);
            (k as f64 + d, b - 0.25 * (a - c) * d)
        };
        RECEIVERS
            .iter()
            .enumerate()
            .map(|(i, p)| {
                let r = norm(sub(*p, SOURCE));
                let exact: Vec<f64> = times.iter().map(|&t| lamb.velocity(r, t)).collect();
                let v = trace(sim, Component::Vz, i);
                let dot = |a: &[f64], b: &[f64]| a.iter().zip(b).map(|(x, y)| x * y).sum::<f64>();
                let ((k, size), (k_exact, size_exact)) = (peak_at(&v), peak_at(&exact));
                [
                    r,
                    misfit(&v, &exact),
                    dot(&v, &exact) / dot(&exact, &exact),
                    k - k_exact,
                    size / size_exact,
                ]
            })
            .collect()
    }

    /// Each receiver's distance and the misfit of its vz trace against Pekeris's.
    pub fn misfits(sim: &Simulation) -> Vec<(f64, f64)> {
        let lamb = lamb();
        let times = sim.sample_times();
        RECEIVERS
            .iter()
            .enumerate()
            .map(|(i, p)| {
                let r = norm(sub(*p, SOURCE));
                let exact: Vec<f64> = times.iter().map(|&t| lamb.velocity(r, t)).collect();
                (r, misfit(&trace(sim, Component::Vz, i), &exact))
            })
            .collect()
    }
}
