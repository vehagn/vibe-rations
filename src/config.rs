//! The TOML run configuration. `examples/` has annotated configurations.

use std::fs;
use std::path::{Path, PathBuf};

use anyhow::{Context, Result, bail, ensure};
use serde::{Deserialize, Serialize};

use crate::grid::Grid;
use crate::receiver::Component;
use crate::source::{Mechanism, Wavelet, double_couple};

#[derive(Debug, Clone, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Config {
    pub grid: GridConfig,
    pub time: TimeConfig,
    #[serde(default)]
    pub boundary: BoundaryConfig,
    #[serde(default)]
    pub scheme: SchemeConfig,
    pub model: ModelConfig,
    #[serde(default, rename = "source")]
    pub sources: Vec<SourceConfig>,
    #[serde(default)]
    pub receivers: ReceiversConfig,
    #[serde(default, rename = "snapshot")]
    pub snapshots: Vec<SnapshotConfig>,
    #[serde(default)]
    pub output: OutputConfig,
    /// What `invert` does; nothing else reads it.
    pub inversion: Option<InversionConfig>,
}

impl Config {
    pub fn from_toml(text: &str) -> Result<Self> {
        Ok(toml::from_str(text)?)
    }

    pub fn load(path: &Path) -> Result<Self> {
        let text =
            fs::read_to_string(path).with_context(|| format!("reading {}", path.display()))?;
        Self::from_toml(&text).with_context(|| format!("in {}", path.display()))
    }

    /// The same run with only source `s`: one shot of a survey whose sources fire one at a
    /// time, as an inversion takes them.
    pub fn shot(&self, s: usize) -> Result<Self> {
        ensure!(
            s < self.sources.len(),
            "there is no source {s}; the configuration has {}",
            self.sources.len()
        );
        Ok(Self {
            sources: vec![self.sources[s].clone()],
            ..self.clone()
        })
    }
}

/// The directory, under a run's output directory, of shot `s`'s results.
pub fn shot_dir(s: usize) -> String {
    format!("shot-{s:03}")
}

#[derive(Debug, Clone, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct GridConfig {
    /// Nodes along x, y and z.
    pub points: [usize; 3],
    /// Node spacing in metres: one number, or one per axis.
    pub spacing: Spacing,
}

#[derive(Debug, Clone, Copy, Deserialize)]
#[serde(untagged)]
pub enum Spacing {
    Uniform(f64),
    PerAxis([f64; 3]),
}

impl GridConfig {
    pub fn grid(&self) -> Result<Grid> {
        let h = match self.spacing {
            Spacing::Uniform(h) => [h; 3],
            Spacing::PerAxis(h) => h,
        };
        ensure!(h.iter().all(|&h| h > 0.0), "grid spacing must be positive");
        ensure!(
            self.points.iter().all(|&n| n >= 8),
            "the grid needs at least 8 nodes along each axis"
        );
        Ok(Grid::new(self.points, h))
    }
}

#[derive(Debug, Clone, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct TimeConfig {
    /// Simulated time, in seconds.
    pub duration: f64,
    /// Time step, in seconds. By default `courant` times the stability limit.
    pub dt: Option<f64>,
    /// Fraction of the stability limit to step at when `dt` is not given; by default
    /// [`default_courant`] for the stencil order.
    pub courant: Option<f64>,
}

/// The time step that suits the scheme, as a fraction of the stability limit. With the time
/// error removed from the traces (`remove_time_dispersion`, the default) a large step costs no
/// accuracy. Without it, up to fourth order a large step is still best, because its time error
/// partly cancels the spatial one, but longer stencils need a small one (README, "Accuracy").
pub fn default_courant(order: usize, remove_time_dispersion: bool) -> f64 {
    if order <= 4 || remove_time_dispersion {
        0.85
    } else {
        0.3
    }
}

#[derive(Debug, Clone, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct BoundaryConfig {
    /// Thickness of the absorbing layers, in nodes.
    #[serde(default = "default_pml_width")]
    pub pml_width: usize,
    /// A free surface at z = 0 instead of an absorbing layer.
    #[serde(default = "default_true")]
    pub free_surface: bool,
    /// Theoretical reflection coefficient of the layers; by default from their width.
    pub pml_reflection: Option<f64>,
    #[serde(default = "default_kappa_max")]
    pub pml_kappa_max: f64,
    #[serde(default = "default_pml_power")]
    pub pml_power: f64,
    /// Frequency, in Hz, that sets α_max = π f; by default the sources' dominant frequency.
    pub pml_alpha_frequency: Option<f64>,
    /// Plain damping in the side layers, with a free surface, as a fraction of theirs: it keeps
    /// them stable for surface waves in soft, high-Poisson's-ratio media. 0 turns it off.
    #[serde(default = "default_pml_sponge")]
    pub pml_sponge: f64,
    /// The P-wave speed that sets the layers' damping, in m/s; by default the model's
    /// fastest. An inversion holds it at the starting model's, so the layers don't change
    /// as the model does.
    pub pml_vp_max: Option<f64>,
}

impl Default for BoundaryConfig {
    fn default() -> Self {
        Self {
            pml_width: default_pml_width(),
            free_surface: true,
            pml_reflection: None,
            pml_kappa_max: default_kappa_max(),
            pml_power: default_pml_power(),
            pml_alpha_frequency: None,
            pml_sponge: default_pml_sponge(),
            pml_vp_max: None,
        }
    }
}

fn default_pml_width() -> usize {
    20
}
fn default_true() -> bool {
    true
}
fn default_kappa_max() -> f64 {
    1.0
}
fn default_pml_power() -> f64 {
    2.0
}
fn default_pml_sponge() -> f64 {
    0.2
}

/// The numerical scheme.
#[derive(Debug, Clone, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct SchemeConfig {
    /// Order of the staggered spatial derivatives: 2, 4, ..., 16. Eighth order with the time
    /// step at 0.3 of the limit needs about half the grid points per wavelength of fourth
    /// order for the same accuracy (README, "Accuracy").
    #[serde(default = "default_order")]
    pub order: usize,
    /// Remove leapfrog's time error with the time-dispersion transform (Wang & Xu 2015;
    /// Koene et al. 2018): pre-warp the source wavelets and un-warp the traces after the run
    /// (`tdt`). The traces are then as if time were continuous, so a large time step costs no
    /// accuracy. On by default; snapshots stay warped, and traces are final only after the
    /// last step. `false` records what leapfrog computes, at a small time step.
    #[serde(default = "default_true")]
    pub remove_time_dispersion: bool,
    /// The derivative coefficients: optimised for every wavenumber the run's grid and wavelets
    /// need (`Stencil::optimised`, the default), or Taylor's, exact for the longest waves.
    #[serde(default)]
    pub coefficients: Coefficients,
    /// Optimise for this many grid points per shortest wavelength instead of the run's own,
    /// which depends on the model's slowest wave. An inversion should fix it, so that the
    /// operator stays the same while the model changes.
    pub optimised_for: Option<f64>,
}

impl Default for SchemeConfig {
    fn default() -> Self {
        Self {
            order: default_order(),
            remove_time_dispersion: true,
            coefficients: Coefficients::default(),
            optimised_for: None,
        }
    }
}

/// How the stencil's coefficients are chosen.
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq, Deserialize, Serialize)]
#[serde(rename_all = "lowercase")]
pub enum Coefficients {
    /// Exact for polynomials up to the order: best for the longest waves.
    Taylor,
    /// Least phase error per node over every direction and every wavenumber up to the run's
    /// shortest wavelength, in the manner of Holberg (1987), with leapfrog's time error
    /// assumed removed.
    #[default]
    Optimised,
}

fn default_order() -> usize {
    8
}

/// The elastic model: horizontal layers or `.npy` files, then optional spheres on top.
#[derive(Debug, Clone, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct ModelConfig {
    #[serde(default, rename = "layer")]
    pub layers: Vec<LayerConfig>,
    pub files: Option<ModelFiles>,
    #[serde(default, rename = "sphere")]
    pub spheres: Vec<SphereConfig>,
}

/// A layer from depth `top` down to the next layer's top. Speeds in m/s, density in kg/m³.
#[derive(Debug, Clone, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct LayerConfig {
    pub top: f64,
    pub vp: f64,
    pub vs: f64,
    pub rho: f64,
}

/// `.npy` arrays of shape (nz, ny, nx), float32 or float64, relative to the configuration.
#[derive(Debug, Clone, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct ModelFiles {
    pub vp: PathBuf,
    pub vs: PathBuf,
    pub rho: PathBuf,
}

#[derive(Debug, Clone, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct SphereConfig {
    pub center: [f64; 3],
    pub radius: f64,
    pub vp: f64,
    pub vs: f64,
    pub rho: f64,
}

/// A point source: exactly one of `explosion`, `moment`, `double_couple` and `force`.
///
/// The wavelet is the moment rate of a moment-tensor source and the force of a force
/// source. A Ricker wavelet peaks at 1, so its amplitude is the peak moment rate (N·m/s) or
/// force (N); a Gaussian has unit area, so its amplitude is the final moment (N·m) or impulse.
#[derive(Debug, Clone, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct SourceConfig {
    /// Position in metres; snapped to the nearest node.
    pub position: [f64; 3],
    pub wavelet: WaveletConfig,
    /// An isotropic moment tensor with this amplitude on the diagonal.
    pub explosion: Option<f64>,
    pub moment: Option<MomentConfig>,
    pub double_couple: Option<DoubleCoupleConfig>,
    /// Force along x, y and z.
    pub force: Option<[f64; 3]>,
}

impl SourceConfig {
    pub fn mechanism(&self) -> Result<Mechanism> {
        let given = [
            self.explosion.is_some(),
            self.moment.is_some(),
            self.double_couple.is_some(),
            self.force.is_some(),
        ];
        match given.iter().filter(|&&g| g).count() {
            1 => {}
            0 => bail!("a source needs one of explosion, moment, double_couple and force"),
            _ => bail!("a source takes only one of explosion, moment, double_couple and force"),
        }
        Ok(if let Some(m0) = self.explosion {
            Mechanism::Moment([m0, m0, m0, 0.0, 0.0, 0.0])
        } else if let Some(m) = &self.moment {
            Mechanism::Moment([m.xx, m.yy, m.zz, m.xy, m.xz, m.yz])
        } else if let Some(dc) = &self.double_couple {
            Mechanism::Moment(double_couple(dc.strike, dc.dip, dc.rake, dc.m0))
        } else {
            Mechanism::Force(self.force.expect("checked above"))
        })
    }
}

#[derive(Debug, Clone, Deserialize)]
#[serde(tag = "kind", rename_all = "snake_case", deny_unknown_fields)]
pub enum WaveletConfig {
    /// Peak frequency `frequency`; centred at `delay` seconds, by default 1.2 / frequency.
    Ricker { frequency: f64, delay: Option<f64> },
    /// exp(-(π f (t - delay))²), scaled to unit area.
    Gaussian { frequency: f64, delay: Option<f64> },
}

impl WaveletConfig {
    pub fn wavelet(&self) -> Result<Wavelet> {
        let (WaveletConfig::Ricker { frequency, delay }
        | WaveletConfig::Gaussian { frequency, delay }) = *self;
        ensure!(frequency > 0.0, "wavelet frequency must be positive");
        let t0 = delay.unwrap_or(1.2 / frequency);
        Ok(match self {
            WaveletConfig::Ricker { .. } => Wavelet::Ricker { f0: frequency, t0 },
            WaveletConfig::Gaussian { .. } => Wavelet::Gaussian { f0: frequency, t0 },
        })
    }
}

/// Moment-tensor components in N·m (or N·m/s with a Ricker wavelet).
#[derive(Debug, Clone, Default, Deserialize)]
#[serde(deny_unknown_fields, default)]
pub struct MomentConfig {
    pub xx: f64,
    pub yy: f64,
    pub zz: f64,
    pub xy: f64,
    pub xz: f64,
    pub yz: f64,
}

/// A shear dislocation, angles in degrees (Aki & Richards; x north, y east, z down).
#[derive(Debug, Clone, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct DoubleCoupleConfig {
    pub strike: f64,
    pub dip: f64,
    pub rake: f64,
    pub m0: f64,
}

#[derive(Debug, Clone, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct ReceiversConfig {
    #[serde(default = "default_components")]
    pub components: Vec<Component>,
    #[serde(default, rename = "line")]
    pub lines: Vec<LineConfig>,
    #[serde(default, rename = "point")]
    pub points: Vec<PointConfig>,
}

impl Default for ReceiversConfig {
    fn default() -> Self {
        Self {
            components: default_components(),
            lines: Vec::new(),
            points: Vec::new(),
        }
    }
}

fn default_components() -> Vec<Component> {
    vec![Component::Vx, Component::Vy, Component::Vz]
}

/// `count` receivers evenly spaced from `from` to `to`, both included.
#[derive(Debug, Clone, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct LineConfig {
    pub from: [f64; 3],
    pub to: [f64; 3],
    pub count: usize,
}

#[derive(Debug, Clone, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct PointConfig {
    pub position: [f64; 3],
}

impl ReceiversConfig {
    /// Every receiver position, lines first.
    pub fn positions(&self) -> Result<Vec<[f64; 3]>> {
        let mut out = Vec::new();
        for line in &self.lines {
            ensure!(
                line.count >= 1,
                "a receiver line needs a count of at least 1"
            );
            for r in 0..line.count {
                let f = if line.count == 1 {
                    0.0
                } else {
                    r as f64 / (line.count - 1) as f64
                };
                out.push(std::array::from_fn(|a| {
                    line.from[a] + f * (line.to[a] - line.from[a])
                }));
            }
        }
        out.extend(self.points.iter().map(|p| p.position));
        Ok(out)
    }
}

#[derive(Debug, Clone, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct SnapshotConfig {
    pub field: Component,
    pub plane: PlaneKind,
    /// Coordinate of the plane along its normal, in metres; unused for a volume.
    #[serde(default)]
    pub at: f64,
    /// Time between snapshots, in seconds.
    pub interval: f64,
    #[serde(default = "default_formats")]
    pub formats: Vec<Format>,
    /// Colour scale limit for PNGs; by default the 99.5th percentile of the snapshot's
    /// magnitudes.
    pub clip: Option<f64>,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Deserialize)]
#[serde(rename_all = "lowercase")]
pub enum PlaneKind {
    Xy,
    Xz,
    Yz,
    Volume,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Deserialize)]
#[serde(rename_all = "lowercase")]
pub enum Format {
    Png,
    Npy,
}

fn default_formats() -> Vec<Format> {
    vec![Format::Png]
}

/// A model parameter an inversion can update.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Deserialize, Serialize)]
#[serde(rename_all = "lowercase")]
pub enum Parameter {
    Vp,
    Vs,
    Rho,
}

impl Parameter {
    pub fn name(self) -> &'static str {
        match self {
            Parameter::Vp => "vp",
            Parameter::Vs => "vs",
            Parameter::Rho => "rho",
        }
    }
}

/// What an inversion's misfit compares (`crate::misfit`).
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default, Deserialize, Serialize)]
#[serde(rename_all = "lowercase")]
pub enum MisfitKind {
    /// Half the weighted sum of squared differences, sample by sample.
    #[default]
    L2,
    /// Half the sum of each trace's squared cross-correlation time shift (Luo & Schuster
    /// 1991).
    Traveltime,
    /// The time–frequency phase misfit (Fichtner et al. 2008).
    Phase,
}

impl MisfitKind {
    pub fn name(self) -> &'static str {
        match self {
            MisfitKind::L2 => "l2",
            MisfitKind::Traveltime => "traveltime",
            MisfitKind::Phase => "phase",
        }
    }
}

/// One stage of an inversion: its iterations, on traces low-pass filtered up to
/// `max_frequency`, with its own misfit if it names one.
#[derive(Debug, Clone, PartialEq, Deserialize, Serialize)]
#[serde(deny_unknown_fields)]
pub struct StageConfig {
    /// The corner of the low-pass filter applied to observed and synthetic traces alike, in
    /// Hz: where its amplitude halves. None: the traces as recorded.
    pub max_frequency: Option<f64>,
    pub iterations: usize,
    /// The stage's misfit, instead of the inversion's.
    pub misfit: Option<MisfitKind>,
}

/// An L-BFGS inversion of the observed traces for the model (`invert`).
#[derive(Debug, Clone, Deserialize, Serialize)]
#[serde(deny_unknown_fields)]
pub struct InversionConfig {
    /// The observed traces, relative to the configuration: `<component>.npy` as `run` writes
    /// them for one source; with several, each source is a shot, in `shot-NNN/` as
    /// `run --shots` writes them. They must be sampled like the synthetics: same dt.
    pub observed: PathBuf,
    /// What to update; the rest stays as the configuration's model has it.
    #[serde(default = "default_parameters")]
    pub parameters: Vec<Parameter>,
    /// Iterations, when there are no stages: 10 by default.
    pub iterations: Option<usize>,
    /// What the misfit compares: `l2` (the default), `traveltime` or `phase`.
    #[serde(default)]
    pub misfit: MisfitKind,
    /// For the `traveltime` misfit: the largest time shift sought, in seconds.
    pub max_shift: Option<f64>,
    /// The order of the stages' Butterworth low-pass filters (each of its two passes).
    #[serde(default = "default_filter_order")]
    pub filter_order: usize,
    /// Frequency continuation (Bunks et al. 1995): stages run in turn, each with its own
    /// low-pass filter, iterations and misfit, L-BFGS starting afresh in each. Without them,
    /// one stage of `iterations` on the traces as recorded.
    #[serde(default, rename = "stage")]
    pub stages: Vec<StageConfig>,
    /// Model changes and gradient changes L-BFGS remembers.
    #[serde(default = "default_memory")]
    pub memory: usize,
    /// The largest change of any value in the first iteration (and after a restart), as a
    /// fraction of its starting value: the scale L-BFGS starts from.
    #[serde(default = "default_initial_step")]
    pub initial_step: f64,
    /// No value changes by more than this fraction of its starting value in one iteration.
    #[serde(default = "default_max_step")]
    pub max_step: f64,
    /// Weigh each component of each shot by 1/(its observed energy), so pressure and
    /// velocity, and loud and quiet shots, count alike.
    #[serde(default = "default_true")]
    pub balance: bool,
    /// Hold the model fixed above this depth in metres, e.g. a known water layer.
    pub fixed_above: Option<f64>,
    /// Divide the gradient by the forward wavefield's energy at each node (a source-side
    /// pseudo-Hessian), so that poorly lit nodes aren't starved and the sources' surroundings
    /// don't take every update.
    #[serde(default = "default_true")]
    pub precondition: bool,
    /// The fraction of the largest energy added to every node's before dividing by it.
    #[serde(default = "default_water_level")]
    pub water_level: f64,
    /// Also divide by the receivers' side of the illumination: the energy, at each node, of
    /// waves sent out from every receiver at once with random signs, averaged over this many
    /// runs at the start of each stage (a receiver-side pseudo-Hessian). 0, the default, for
    /// the source side alone.
    #[serde(default)]
    pub receiver_runs: usize,
    /// Factors on the preconditioner for each inverted parameter: a change of variables that
    /// scales that parameter's steps, e.g. `rho = 0.3` where density takes up what the
    /// velocities should explain. 1 by default.
    #[serde(default)]
    pub scale: ParameterScale,
    /// Bounds on the inverted values; by default half the smallest (non-zero, for vs) and twice
    /// the largest of the starting model's, and for vp no faster than the time step allows.
    /// λ ≥ 0 too: vs stays below vp/√2, or vp above √2 vs where vs isn't inverted. The
    /// starting model must lie within them.
    pub vp_range: Option<[f64; 2]>,
    pub vs_range: Option<[f64; 2]>,
    pub rho_range: Option<[f64; 2]>,
}

/// A factor for each parameter an inversion can update.
#[derive(Debug, Clone, Copy, PartialEq, Deserialize, Serialize)]
#[serde(deny_unknown_fields)]
pub struct ParameterScale {
    #[serde(default = "one")]
    pub vp: f64,
    #[serde(default = "one")]
    pub vs: f64,
    #[serde(default = "one")]
    pub rho: f64,
}

impl Default for ParameterScale {
    fn default() -> Self {
        Self {
            vp: 1.0,
            vs: 1.0,
            rho: 1.0,
        }
    }
}

impl ParameterScale {
    pub fn of(&self, p: Parameter) -> f64 {
        match p {
            Parameter::Vp => self.vp,
            Parameter::Vs => self.vs,
            Parameter::Rho => self.rho,
        }
    }

    /// Whether every factor is 1.
    pub fn is_unit(&self) -> bool {
        *self == Self::default()
    }
}

fn one() -> f64 {
    1.0
}

impl InversionConfig {
    /// The stages to run, each with its misfit: the configured ones, or one of `iterations`
    /// (10 by default) on the traces as recorded.
    pub fn stages(&self) -> Result<Vec<StageConfig>> {
        if self.stages.is_empty() {
            return Ok(vec![StageConfig {
                max_frequency: None,
                iterations: self.iterations.unwrap_or(10),
                misfit: Some(self.misfit),
            }]);
        }
        ensure!(
            self.iterations.is_none(),
            "give the iterations of each [[inversion.stage]], not of [inversion]"
        );
        Ok(self
            .stages
            .iter()
            .map(|s| StageConfig {
                misfit: Some(s.misfit.unwrap_or(self.misfit)),
                ..s.clone()
            })
            .collect())
    }
}

fn default_parameters() -> Vec<Parameter> {
    vec![Parameter::Vp, Parameter::Vs]
}
fn default_filter_order() -> usize {
    4
}
fn default_memory() -> usize {
    5
}
fn default_initial_step() -> f64 {
    0.02
}
fn default_max_step() -> f64 {
    0.1
}
fn default_water_level() -> f64 {
    1e-3
}

#[derive(Debug, Clone, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct OutputConfig {
    /// Where results go, relative to the working directory.
    #[serde(default = "default_out")]
    pub dir: PathBuf,
    /// Record the receivers every this many steps.
    #[serde(default = "default_record_every")]
    pub record_every: usize,
}

impl Default for OutputConfig {
    fn default() -> Self {
        Self {
            dir: default_out(),
            record_every: default_record_every(),
        }
    }
}

fn default_out() -> PathBuf {
    PathBuf::from("out")
}
fn default_record_every() -> usize {
    1
}

#[cfg(test)]
mod tests {
    use super::*;

    const MINIMAL: &str = r#"
        [grid]
        points = [20, 20, 20]
        spacing = 10.0
        [time]
        duration = 0.1
        [[model.layer]]
        top = 0.0
        vp = 2000.0
        vs = 1000.0
        rho = 2000.0
        [[source]]
        position = [100.0, 100.0, 100.0]
        explosion = 1.0
        wavelet = { kind = "ricker", frequency = 10.0 }
    "#;

    #[test]
    fn a_shot_keeps_one_source() {
        let two = MINIMAL.to_string()
            + r#"
        [[source]]
        position = [50.0, 100.0, 100.0]
        force = [0.0, 0.0, 1.0]
        wavelet = { kind = "ricker", frequency = 12.0 }
        "#;
        let c = Config::from_toml(&two).unwrap();
        let shot = c.shot(1).unwrap();
        assert_eq!(shot.sources.len(), 1);
        assert_eq!(shot.sources[0].position, [50.0, 100.0, 100.0]);
        assert_eq!(shot.grid.points, c.grid.points);
        assert!(c.shot(2).is_err());
        assert_eq!(shot_dir(7), "shot-007");
    }

    #[test]
    fn fills_in_defaults() {
        let c = Config::from_toml(MINIMAL).unwrap();
        assert_eq!(c.boundary.pml_width, 20);
        assert!(c.boundary.free_surface);
        assert_eq!(c.scheme.order, 8);
        assert_eq!(c.time.courant, None);
        assert!(c.scheme.remove_time_dispersion);
        assert_eq!(c.scheme.coefficients, Coefficients::Optimised);
        assert_eq!(c.scheme.optimised_for, None);
        assert_eq!(default_courant(8, false), 0.3);
        assert_eq!(default_courant(4, false), 0.85);
        assert_eq!(default_courant(8, true), 0.85);
        assert_eq!(c.receivers.components.len(), 3);
        assert_eq!(c.output.dir, PathBuf::from("out"));
        let w = c.sources[0].wavelet.wavelet().unwrap();
        assert!(matches!(w, Wavelet::Ricker { f0, t0 } if f0 == 10.0 && (t0 - 0.12).abs() < 1e-12));
        assert!(matches!(c.sources[0].mechanism().unwrap(), Mechanism::Moment(m) if m[0] == 1.0));
    }

    #[test]
    fn rejects_typos() {
        let text = MINIMAL.replace("duration", "durration");
        assert!(Config::from_toml(&text).is_err());
    }

    #[test]
    fn a_source_needs_exactly_one_mechanism() {
        let both = MINIMAL.replace(
            "explosion = 1.0",
            "explosion = 1.0\nforce = [0.0, 0.0, 1.0]",
        );
        let c = Config::from_toml(&both).unwrap();
        assert!(c.sources[0].mechanism().is_err());
        let none = MINIMAL.replace("explosion = 1.0", "");
        let c = Config::from_toml(&none).unwrap();
        assert!(c.sources[0].mechanism().is_err());
    }

    #[test]
    fn spacing_can_differ_per_axis() {
        let text = MINIMAL.replace("spacing = 10.0", "spacing = [10.0, 20.0, 5.0]");
        let g = Config::from_toml(&text).unwrap().grid.grid().unwrap();
        assert_eq!(g.h, [10.0, 20.0, 5.0]);
    }

    #[test]
    fn spreads_receivers_along_lines() {
        let r = ReceiversConfig {
            lines: vec![LineConfig {
                from: [0.0, 0.0, 0.0],
                to: [100.0, 50.0, 0.0],
                count: 3,
            }],
            points: vec![PointConfig {
                position: [1.0, 2.0, 3.0],
            }],
            ..Default::default()
        };
        let p = r.positions().unwrap();
        assert_eq!(
            p,
            vec![
                [0.0, 0.0, 0.0],
                [50.0, 25.0, 0.0],
                [100.0, 50.0, 0.0],
                [1.0, 2.0, 3.0]
            ]
        );
    }
}
