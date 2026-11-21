//! Plane waves at normal incidence on a planar interface between two solids, and their exact
//! reflection and transmission.
//!
//! A uniform plane of sources makes the wavefield depend on z alone: a plane of σzz sources
//! sends P waves, a plane of x forces S waves. Each side of the source plane takes half the
//! jump the source makes, so a moment $m_0$ per node ($m_0/h^2$ per unit area) sends
//! $v_z = \pm m_0 w(t) / (2 \rho \alpha^2 h^2)$ and a force $f$ per node sends
//! $v_x = f w(t) / (2 \rho \beta h^2)$. At the interface particle velocity is reflected by
//! $R = (Z_1 - Z_2)/(Z_1 + Z_2)$ and transmitted by $T = 1 + R$, with $Z = \rho\alpha$ for P
//! and $\rho\beta$ for S; nothing converts at normal incidence.
//!
//! The run keeps every z-plane exactly uniform: the values outside the grid's sides copy the
//! nearest node's, so no x or y derivative sees an edge, and the absorbing layers along x and
//! y, which act only on those derivatives, do nothing. The grid needs only to be as wide as
//! its layers allow.

use vibe_rations::config::{Config, Interfaces};
use vibe_rations::grid::Grid;
use vibe_rations::model::Model;
use vibe_rations::receiver::Component;
use vibe_rations::simulation::Simulation;
use vibe_rations::solver::{Field, Fields, PointOp};
use vibe_rations::source::Wavelet;

use super::trace;

/// The wave a plane of sources sends.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Wave {
    P,
    S,
}

/// The media of the interface tests, [vp, vs, ρ]: R = -0.304 for P, -0.238 for S.
pub const UPPER: [f64; 3] = [2000.0, 1155.0, 2000.0];
pub const LOWER: [f64; 3] = [3000.0, 1500.0, 2500.0];

const MOMENT: f64 = 1e13;
const FORCE: f64 = 1e9;
const PML: usize = 10;

/// A plane wave sent down from `source` onto an interface at `interface`, recorded at
/// `above` (between the two: the direct and reflected waves) and `below` (the transmitted
/// wave). Media are [vp, vs, ρ]; depths in m.
#[derive(Clone, Debug)]
pub struct PlaneWave {
    pub wave: Wave,
    pub spacing: f64,
    pub f0: f64,
    pub order: usize,
    /// How the layers meet the nodes when [`PlaneWave::run`] builds them.
    pub interfaces: Interfaces,
    pub upper: [f64; 3],
    pub lower: [f64; 3],
    pub source: f64,
    pub interface: f64,
    pub above: f64,
    pub below: f64,
}

/// Arrivals of the wavelet: (amplitude, delay).
pub type Arrivals = Vec<(f64, f64)>;

/// The traces of a run, and their times.
pub struct Traces {
    pub above: Vec<f64>,
    pub below: Vec<f64>,
    pub times: Vec<f64>,
}

/// How a trace compares with one exact arrival, over that arrival's window.
#[derive(Clone, Copy, Debug)]
pub struct Fit {
    /// ‖simulated - exact‖ / ‖arrival‖, the other arrivals' exact values taken off first.
    pub misfit: f64,
    /// The least-squares amplitude of the simulated arrival relative to the exact one.
    pub ratio: f64,
    /// How much later the simulated arrival comes, from the correlation's peak, in s.
    pub shift: f64,
}

impl PlaneWave {
    /// The interface tests' setup on a 10 m grid at `scale`: 10/`scale` Hz, distances in
    /// nodes multiplied by `scale`, so that `scale` multiplies the points per wavelength
    /// (P 20, S 11.6 at f0 for 1). The direct and reflected waves don't overlap at the
    /// receiver above.
    pub fn interface(wave: Wave, scale: usize, delta: f64) -> Self {
        let s = scale;
        let gaps = match wave {
            Wave::P => [10 * s, 30 * s, 20 * s],
            Wave::S => [10 * s, 17 * s, 20 * s],
        };
        Self::new(wave, 10.0, 10.0 / s as f64, UPPER, LOWER, gaps, delta)
    }

    /// The source a wavelength below the top layer, the upper receiver `gap_above` nodes
    /// below it, the interface `gap` nodes below that plus `delta` of a node, the lower
    /// receiver `gap_below` nodes below the interface's node.
    pub fn new(
        wave: Wave,
        spacing: f64,
        f0: f64,
        upper: [f64; 3],
        lower: [f64; 3],
        [gap_above, gap, gap_below]: [usize; 3],
        delta: f64,
    ) -> Self {
        let k_source = (PML + 5) as f64;
        let k_above = k_source + gap_above as f64;
        let k_interface = k_above + gap as f64;
        Self {
            wave,
            spacing,
            f0,
            order: 8,
            interfaces: Interfaces::Nodes,
            upper,
            lower,
            source: k_source * spacing,
            interface: (k_interface + delta) * spacing,
            above: k_above * spacing,
            below: (k_interface + gap_below as f64) * spacing,
        }
    }

    fn speeds(&self) -> (f64, f64) {
        let i = match self.wave {
            Wave::P => 0,
            Wave::S => 1,
        };
        (self.upper[i], self.lower[i])
    }

    pub fn wavelet(&self) -> Wavelet {
        Wavelet::Ricker {
            f0: self.f0,
            t0: 1.2 / self.f0,
        }
    }

    /// The exact arrivals above and below the interface: (amplitude, delay) of the wavelet.
    /// Above: the direct wave, then the reflected one; below: the transmitted wave.
    pub fn exact(&self) -> (Arrivals, Arrivals) {
        let (c1, c2) = self.speeds();
        let (z1, z2) = (self.upper[2] * c1, self.lower[2] * c2);
        let r = (z1 - z2) / (z1 + z2);
        let h2 = self.spacing * self.spacing;
        let a = match self.wave {
            Wave::P => MOMENT / (2.0 * self.upper[2] * c1 * c1 * h2),
            Wave::S => FORCE / (2.0 * self.upper[2] * c1 * h2),
        };
        let above = vec![
            (a, (self.above - self.source) / c1),
            (
                r * a,
                (2.0 * self.interface - self.source - self.above) / c1,
            ),
        ];
        let below = vec![(
            (1.0 + r) * a,
            (self.interface - self.source) / c1 + (self.below - self.interface) / c2,
        )];
        (above, below)
    }

    /// Long enough for the last arrival to pass.
    pub fn duration(&self) -> f64 {
        let (above, below) = self.exact();
        let last = above.iter().chain(&below).map(|a| a.1).fold(0.0, f64::max);
        last + 2.4 / self.f0 + 0.02
    }

    /// Runs it with the two media as layers, the lower one's top at the interface, meeting the
    /// nodes as `interfaces` says, or with `medium` giving [vp, vs, ρ] at each depth of nodes
    /// instead.
    pub fn run(&self, medium: Option<&dyn Fn(f64) -> [f64; 3]>) -> Traces {
        let h = self.spacing;
        let across = 2 * PML + 5;
        let nz = (self.below / h).round() as usize + 10 + PML;
        let c = (across / 2) as f64 * h;
        let (source, component) = match self.wave {
            Wave::P => (format!("explosion = {MOMENT:e}"), "vz"),
            Wave::S => (format!("force = [{FORCE:e}, 0.0, 0.0]"), "vx"),
        };
        let layer = |top: f64, m: [f64; 3]| {
            format!(
                "[[model.layer]]\ntop = {top}\nvp = {}\nvs = {}\nrho = {}\n",
                m[0], m[1], m[2]
            )
        };
        let toml = format!(
            r#"
            [grid]
            points = [{across}, {across}, {nz}]
            spacing = {h}
            [time]
            duration = {duration}
            [scheme]
            order = {order}
            [boundary]
            pml_width = {PML}
            free_surface = false
            [model]
            interfaces = "{interfaces}"
            {upper}
            {lower}
            [[source]]
            position = [{c}, {c}, {zs}]
            {source}
            wavelet = {{ kind = "ricker", frequency = {f0} }}
            [receivers]
            components = ["{component}"]
            [[receivers.point]]
            position = [{c}, {c}, {za}]
            [[receivers.point]]
            position = [{c}, {c}, {zb}]
            "#,
            duration = self.duration(),
            order = self.order,
            interfaces = match self.interfaces {
                Interfaces::Nodes => "nodes",
                Interfaces::BandLimited => "band-limited",
            },
            upper = layer(0.0, self.upper),
            lower = layer(self.interface, self.lower),
            zs = self.source,
            za = self.above,
            zb = self.below,
            f0 = self.f0,
        );
        let cfg = Config::from_toml(&toml).expect("a valid configuration");
        let mut sim = match medium {
            None => Simulation::new(&cfg, std::path::Path::new(".")).expect("a valid setup"),
            Some(medium) => {
                let grid = Simulation::grid(&cfg).expect("a grid");
                let mut model = Model::homogeneous(grid, 0.0, 0.0, 0.0);
                let plane = grid.n[0] * grid.n[1];
                for k in 0..grid.n[2] {
                    let [vp, vs, rho] = medium(k as f64 * h);
                    let range = k * plane..(k + 1) * plane;
                    model.vp[range.clone()].fill(vp as f32);
                    model.vs[range.clone()].fill(vs as f32);
                    model.rho[range].fill(rho as f32);
                }
                Simulation::with_model(&cfg, model).expect("a valid setup")
            }
        };
        let grid = sim.solver.grid;
        // The point source's weight on each plane of the one field kept, spread evenly over
        // that plane: σzz for P (σxx and σyy would only act at the sides), vx for S.
        let keep = match self.wave {
            Wave::P => Field::Szz,
            Wave::S => Field::Vx,
        };
        let mut weights: Vec<(usize, f32)> = Vec::new();
        for &(f, index, w) in &sim.sources[0].op.terms {
            if f != keep {
                continue;
            }
            let k = index / grid.plane_len() - grid.halo;
            match weights.iter_mut().find(|(kk, _)| *kk == k) {
                Some((_, sum)) => *sum += w,
                None => weights.push((k, w)),
            }
        }
        let mut terms = Vec::new();
        for &(k, w) in &weights {
            for j in 0..grid.n[1] {
                for i in 0..grid.n[0] {
                    terms.push((keep, grid.idx(i, j, k), w));
                }
            }
        }
        sim.sources[0].op = PointOp { terms };
        // Stepped here rather than by `Simulation::run`, to keep the planes uniform.
        let source = &sim.sources[0];
        for _ in 0..sim.steps {
            let n = sim.solver.step;
            sim.solver.step_with(
                |f, t| {
                    if source.is_force() {
                        source.inject(f, n, t);
                    }
                    extend_sideways(f, grid);
                },
                |f, t| {
                    if !source.is_force() {
                        source.inject(f, n, t);
                    }
                    extend_sideways(f, grid);
                },
            );
            if n.is_multiple_of(sim.record_every) {
                sim.receivers.record(&sim.solver);
            }
        }
        sim.finish();
        assert!(
            sim.solver.peak_velocity().is_finite(),
            "the wavefield blew up"
        );
        let c = match self.wave {
            Wave::P => Component::Vz,
            Wave::S => Component::Vx,
        };
        Traces {
            above: trace(&sim, c, 0),
            below: trace(&sim, c, 1),
            times: sim.sample_times(),
        }
    }

    /// The direct, reflected and transmitted waves of `traces` against the exact ones.
    pub fn fits(&self, traces: &Traces) -> [Fit; 3] {
        let (above, below) = self.exact();
        [
            self.fit(&traces.above, &traces.times, &above, 0),
            self.fit(&traces.above, &traces.times, &above, 1),
            self.fit(&traces.below, &traces.times, &below, 0),
        ]
    }

    /// Arrival `which` of `arrivals` against `trace`, over the arrival's window.
    fn fit(&self, trace: &[f64], times: &[f64], arrivals: &[(f64, f64)], which: usize) -> Fit {
        let w = self.wavelet();
        let (a, delay) = arrivals[which];
        let window: Vec<usize> = (0..times.len())
            .filter(|&m| times[m] >= delay && times[m] <= delay + 2.4 / self.f0)
            .collect();
        // The trace with the other arrivals taken off: what this one should account for.
        let rest: Vec<f64> = window
            .iter()
            .map(|&m| {
                let others: f64 = arrivals
                    .iter()
                    .enumerate()
                    .filter(|&(i, _)| i != which)
                    .map(|(_, &(b, d))| b * w.value(times[m] - d))
                    .sum();
                trace[m] - others
            })
            .collect();
        let exact = |m: usize, lag: f64| a * w.value(times[m] - delay - lag);
        let (mut err, mut norm, mut cross) = (0.0, 0.0, 0.0);
        for (&m, &r) in window.iter().zip(&rest) {
            let e = exact(m, 0.0);
            err += (r - e).powi(2);
            norm += e * e;
            cross += r * e;
        }
        // The correlation's peak over lags of up to a quarter period, refined by a parabola.
        let dt = times[1] - times[0];
        let lags = (0.25 / self.f0 / dt).ceil() as i64;
        let corr: Vec<f64> = (-lags..=lags)
            .map(|lag| {
                window
                    .iter()
                    .zip(&rest)
                    .map(|(&m, &r)| r * exact(m, lag as f64 * dt))
                    .sum()
            })
            .collect();
        let best = (0..corr.len())
            .max_by(|&p, &q| corr[p].total_cmp(&corr[q]))
            .expect("lags");
        let mut shift = best as f64 - lags as f64;
        if best > 0 && best + 1 < corr.len() {
            let (p, q, r) = (corr[best - 1], corr[best], corr[best + 1]);
            shift += 0.5 * (p - r) / (p - 2.0 * q + r);
        }
        Fit {
            misfit: (err / norm).sqrt(),
            ratio: cross / norm,
            shift: shift * dt,
        }
    }

    /// The time a wave in the upper medium takes to cross one node spacing: the reflected
    /// wave's shift is twice the interface's error in position over this.
    pub fn node_time(&self) -> f64 {
        self.spacing / self.speeds().0
    }
}

/// Sets the values outside the grid along x and y to those of the nearest node.
fn extend_sideways(fields: &mut Fields, grid: Grid) {
    let [px, py, pz] = grid.padded();
    let [nx, ny, _] = grid.n;
    let l = grid.halo;
    for v in [
        &mut fields.vx,
        &mut fields.vy,
        &mut fields.vz,
        &mut fields.sxx,
        &mut fields.syy,
        &mut fields.szz,
        &mut fields.sxy,
        &mut fields.sxz,
        &mut fields.syz,
    ] {
        for k in 0..pz {
            let base = k * px * py;
            for j in l..l + ny {
                let row = base + j * px;
                let (first, last) = (v[row + l], v[row + l + nx - 1]);
                v[row..row + l].fill(first);
                v[row + l + nx..row + px].fill(last);
            }
            for j in (0..l).chain(l + ny..py) {
                let from = base + j.clamp(l, l + ny - 1) * px;
                v.copy_within(from..from + px, base + j * px);
            }
        }
    }
}
