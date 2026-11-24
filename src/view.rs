//! Live pictures of a run, for the web app: a section of the wavefield and the gather so far,
//! as RGBA pixels (the browser draws them) or PNG (a server sends them). The browser's solver
//! and `vibe-rations serve` share it, so both draw the same pictures. An inversion on
//! synthetic data ([`InversionView`]) is stepped a unit of work at a time and pictured the same
//! way: the model, its update, the true update and the gradient.

use std::path::{Path, PathBuf};
use std::sync::Arc;

use anyhow::{Context, Result, bail, ensure};
use serde_json::json;

use crate::config::{Config, InversionConfig, Parameter, PlaneKind, StageConfig};
use crate::grid::Grid;
use crate::image;
use crate::inversion::{Inversion, Observed};
use crate::model::Model;
use crate::npy;
use crate::receiver::Component;
use crate::simulation::Simulation;
use crate::snapshot::Snapshot;

/// Pixels to draw.
#[derive(Clone, Debug)]
pub struct Frame {
    pub width: usize,
    pub height: usize,
    /// RGBA, row by row from the top.
    pub rgba: Vec<u8>,
    /// The value drawn at full colour: the top of the colour scale.
    pub clip: f32,
    /// The bottom of the colour scale: -clip for signed pictures.
    pub low: f32,
}

impl Frame {
    /// The frame as a PNG file.
    pub fn png(&self) -> Result<Vec<u8>> {
        image::png_rgba(self.width, self.height, &self.rgba)
    }
}

/// A simulation and the views of it.
pub struct View {
    pub sim: Simulation,
    /// The largest colour scale each field has had this run.
    largest_clip: [f32; 4],
}

/// Sections are drawn at about this many pixels along their longer side.
const SECTION_PIXELS: usize = 600;

/// A colour scale never drops below this share of the largest one the field has had, so the
/// picture fades as the waves leave rather than blowing up what little is left.
const CLIP_MEMORY: f32 = 0.05;

impl View {
    pub fn new(sim: Simulation) -> Self {
        Self {
            sim,
            largest_clip: [0.0; 4],
        }
    }

    /// The run's settings: `run.json`'s fields, plus the receivers' positions.
    pub fn info(&self) -> serde_json::Value {
        let mut info = serde_json::to_value(&self.sim.info).expect("RunInfo serialises");
        let g = self.sim.solver.grid;
        let receivers: Vec<[f64; 3]> = self
            .sim
            .receivers
            .points
            .iter()
            .map(|p| std::array::from_fn(|a| p[a] * g.h[a]))
            .collect();
        info["receiver_positions"] = serde_json::json!(receivers);
        info
    }

    /// Takes up to `steps` steps, recording the receivers and, given a directory, writing the
    /// snapshots that fall due into it; returns how many it took. Fails if the wavefield has
    /// blown up. After the last step the traces are final (`Simulation::finish`).
    pub fn advance(&mut self, steps: usize, snapshots: Option<&Path>) -> Result<usize> {
        let n = steps.min(self.sim.remaining());
        for _ in 0..n {
            match snapshots {
                Some(dir) => self.sim.advance_writing(dir)?,
                None => self.sim.advance(),
            }
        }
        self.sim.finish();
        ensure!(
            self.sim.solver.peak_velocity().is_finite(),
            "the wavefield blew up by step {}; check the model, or lower dt or courant",
            self.sim.solver.step
        );
        Ok(n)
    }

    /// `field` (vx, vy, vz or p) on `plane` (xy, xz or yz) through node `index` along its
    /// normal, scaled so `gain` 1 shows the 99.5th percentile at full colour.
    pub fn section(&mut self, field: &str, plane: &str, index: usize, gain: f32) -> Result<Frame> {
        let plane = PlaneKind::plane(plane)?;
        let n = self.sim.solver.grid.n[plane.normal().expect("a plane")];
        ensure!(index < n, "plane index {index} is off the grid (0..{n})");
        ensure!(gain > 0.0, "gain must be positive");
        let component = Component::from_name(field)?;
        let snapshot = Snapshot {
            component,
            plane,
            index,
            every: 1,
            png: false,
            npy: false,
            clip: None,
        };
        let capture = snapshot.capture(&self.sim.solver);
        let largest = &mut self.largest_clip[component as usize];
        let clip = snapshot.clip_for(&capture).max(CLIP_MEMORY * *largest);
        *largest = largest.max(clip);
        let clip = clip / gain;
        let (rows, cols) = (capture.shape[0], capture.shape[1]);
        let factor = (SECTION_PIXELS / rows.max(cols)).clamp(1, 8);
        let values = upsample(&capture.values, cols, rows, factor);
        let shade = upsample(&capture.shade, cols, rows, factor);
        let section = image::Section {
            width: (cols - 1) * factor + 1,
            height: (rows - 1) * factor + 1,
            values: &values,
            clip,
            shade: &shade,
        };
        Ok(Frame {
            width: section.width,
            height: section.height,
            rgba: image::section_rgba(&section),
            clip,
            low: -clip,
        })
    }

    /// The gather of `name` so far, one column per receiver and `rows` rows of time.
    pub fn gather(&self, name: &str, rows: usize) -> Result<Frame> {
        let r = &self.sim.receivers;
        let c = self.recorded(name)?;
        let rows = rows.clamp(1, r.samples.max(1));
        Ok(Frame {
            width: r.points.len(),
            height: rows,
            rgba: image::gather_rgba(r.points.len(), r.samples, r.recorded(), &r.traces[c], rows),
            clip: 1.0,
            low: -1.0,
        })
    }

    /// The traces of `name` as an `.npy` file, (receivers, samples).
    pub fn traces_npy(&self, name: &str) -> Result<Vec<u8>> {
        let r = &self.sim.receivers;
        let i = self.recorded(name)?;
        npy::to_bytes(&[r.points.len(), r.samples], &r.traces[i])
    }

    /// Which of the receivers' trace arrays holds `name`.
    fn recorded(&self, name: &str) -> Result<usize> {
        let c = Component::from_name(name)?;
        self.sim
            .receivers
            .components
            .iter()
            .position(|&x| x == c)
            .with_context(|| format!("the receivers don't record {name}"))
    }
}

/// A section through values on the interior nodes of `grid`, on `plane` (xy, xz or yz) through
/// node `index` along its normal, row-major from the top, and its width and height.
fn cut(grid: &Grid, values: &[f32], plane: &str, index: usize) -> Result<(Vec<f32>, usize, usize)> {
    let normal = PlaneKind::plane(plane)?.normal().expect("a plane");
    // The other two axes, the first across the section and the second down it.
    let [across, down] = match normal {
        0 => [1, 2],
        1 => [0, 2],
        _ => [0, 1],
    };
    let n = grid.n;
    ensure!(
        index < n[normal],
        "plane index {index} is off the grid (0..{})",
        n[normal]
    );
    let (w, h) = (n[across], n[down]);
    let mut out = Vec::with_capacity(w * h);
    for r in 0..h {
        for c in 0..w {
            let mut ijk = [0; 3];
            ijk[normal] = index;
            ijk[across] = c;
            ijk[down] = r;
            out.push(values[grid.interior_idx(ijk[0], ijk[1], ijk[2])]);
        }
    }
    Ok((out, w, h))
}

/// A section scaled up for the screen as wavefield sections are, coloured by `colour`.
fn picture(
    grid: &Grid,
    values: &[f32],
    plane: &str,
    index: usize,
    colour: impl Fn(&[f32]) -> Vec<u8>,
) -> Result<(Vec<u8>, usize, usize)> {
    let (v, w, h) = cut(grid, values, plane, index)?;
    let factor = (SECTION_PIXELS / w.max(h)).clamp(1, 8);
    let up = upsample(&v, w, h, factor);
    Ok((colour(&up), (w - 1) * factor + 1, (h - 1) * factor + 1))
}

/// What a picture of an inversion shows.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Show {
    /// The model as it stands.
    Model,
    /// Its change from the start.
    Update,
    /// The true model less the start: what the update would be if the inversion found it.
    TrueUpdate,
    /// The gradient's opposite, -m₀ ∂χ/∂m: the way the misfit falls fastest.
    Descent,
}

impl Show {
    pub fn parse(name: &str) -> Result<Self> {
        Ok(match name {
            "model" => Show::Model,
            "update" => Show::Update,
            "true-update" => Show::TrueUpdate,
            "descent" => Show::Descent,
            _ => bail!("unknown picture {name}; use model, update, true-update or descent"),
        })
    }
}

/// An inversion's models and gradient as they last stood, for pictures. Cloning shares them,
/// so a server can keep a copy to draw from while the inversion moves on.
#[derive(Clone)]
pub struct Pictures {
    start: Arc<Model>,
    truth: Option<Arc<Model>>,
    /// The model as it stands; the start's own until the first iteration.
    model: Arc<Model>,
    /// m₀ ∂χ/∂m for vp, vs and ρ, where inverted, once computed.
    gradient: [Option<Arc<Vec<f32>>>; 3],
}

impl Pictures {
    /// `parameter` (vp, vs or rho) on `plane` (xy, xz or yz) through node `index`, showing
    /// `show` (model, update, true-update or descent). Models are coloured from the smallest
    /// to the largest value of the start and the truth, updates on one scale, the largest true
    /// change (or else the largest change so far) divided by `gain`, and the gradient on its
    /// largest value divided by `gain`.
    pub fn section(
        &self,
        parameter_name: &str,
        show: &str,
        plane: &str,
        index: usize,
        gain: f32,
    ) -> Result<Frame> {
        let p = Parameter::from_name(parameter_name)?;
        ensure!(gain > 0.0, "gain must be positive");
        let grid = &self.start.grid;
        let start = self.start.values(p);
        let largest = |v: &mut dyn Iterator<Item = f32>| v.fold(0.0f32, |m, x| m.max(x.abs()));
        let difference =
            |m: &Model| -> Vec<f32> { m.values(p).iter().zip(start).map(|(a, b)| a - b).collect() };
        let true_change = self.truth.as_ref().map(|t| difference(t));
        let update_clip = || {
            let true_largest = true_change
                .as_ref()
                .map_or(0.0, |d| largest(&mut d.iter().copied()));
            if true_largest > 0.0 {
                true_largest
            } else {
                largest(&mut difference(&self.model).into_iter())
            }
        };
        let signed = |v: &[f32], clip: f32| -> Result<Frame> {
            let clip = if clip > 0.0 { clip / gain } else { 1.0 };
            let (rgba, width, height) =
                picture(grid, v, plane, index, |u| image::signed_rgba(u, clip))?;
            Ok(Frame {
                width,
                height,
                rgba,
                clip,
                low: -clip,
            })
        };
        match Show::parse(show)? {
            Show::Model => {
                let mut range = [f32::INFINITY, f32::NEG_INFINITY];
                for m in std::iter::once(&self.start).chain(&self.truth) {
                    for &v in m.values(p) {
                        range = [range[0].min(v), range[1].max(v)];
                    }
                }
                let [mut lo, mut hi] = range;
                if hi <= lo {
                    // A homogeneous model: a band around it.
                    (lo, hi) = (0.95 * lo, 1.05 * hi.max(1e-3));
                }
                let (rgba, width, height) =
                    picture(grid, self.model.values(p), plane, index, |u| {
                        image::sequential_rgba(u, lo, hi)
                    })?;
                Ok(Frame {
                    width,
                    height,
                    rgba,
                    clip: hi,
                    low: lo,
                })
            }
            Show::Update => signed(&difference(&self.model), update_clip()),
            Show::TrueUpdate => {
                let d = true_change
                    .as_ref()
                    .context("this inversion has no true model")?;
                signed(d, update_clip())
            }
            Show::Descent => {
                let g = self.gradient[p as usize]
                    .as_ref()
                    .with_context(|| format!("there is no gradient for {parameter_name} (yet)"))?;
                let descent: Vec<f32> = g.iter().map(|v| -v).collect();
                signed(&descent, largest(&mut descent.iter().copied()))
            }
        }
    }

    /// The model as it stands, `parameter` (vp, vs or rho) as an `.npy` file shaped (nz, ny,
    /// nx) like model files.
    pub fn npy(&self, parameter_name: &str) -> Result<Vec<u8>> {
        npy::to_bytes(
            &self.model.grid.shape(),
            self.model.values(Parameter::from_name(parameter_name)?),
        )
    }

    /// Bytes held: each model and gradient once, however many share it.
    pub fn memory_bytes(&self) -> usize {
        let mut bytes =
            self.start.memory_bytes() + self.truth.as_deref().map_or(0, Model::memory_bytes);
        if !Arc::ptr_eq(&self.model, &self.start) {
            bytes += self.model.memory_bytes();
        }
        bytes
            + self
                .gradient
                .iter()
                .flatten()
                .map(|g| 4 * g.len())
                .sum::<usize>()
    }
}

/// Where an inversion on synthetic data is: each stage is a unit of work.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Stage {
    /// The true model's traces of this shot (from 0).
    Observing(usize),
    /// The misfit and gradient at the start.
    Starting,
    /// L-BFGS iteration this (from 1, counted over every stage).
    Iterating(usize),
    /// The start of this stage of frequency continuation (from 1; stage 0 starts with the
    /// inversion): its misfit and gradient on its own filtered traces.
    Staging(usize),
    Done,
}

/// An inversion of synthetic data: a start (a configuration with `[inversion]`) and a true
/// model (the same configuration with another model), whose traces are the observed data. The
/// browser's solver and `vibe-rations serve` step it a unit at a time ([`Self::advance`]) so
/// they can report progress, picture it and stop between units.
pub struct InversionView {
    cfg: Config,
    truth: Config,
    base: PathBuf,
    settings: InversionConfig,
    /// The inversion's stages, each with its iterations and misfit.
    stages: Vec<StageConfig>,
    observed: Vec<Observed>,
    inversion: Option<Inversion>,
    stage: Stage,
    /// Iterations taken in all, and in the stage under way.
    iterations: usize,
    in_stage: usize,
    /// Why it stopped before its last iteration, if it did.
    stopped: Option<String>,
    first_misfit: Option<f64>,
    pub pictures: Pictures,
}

impl InversionView {
    /// Sets up nothing yet but the two models: `cfg` is the start, with its `[inversion]`;
    /// `truth` the same survey with the true model. Both must be sampled alike (the same time
    /// step: `settings::synthetic` in the web crate pins it for both).
    pub fn new(cfg: &Config, truth: &Config, base: &Path) -> Result<Self> {
        // No iterations (the page's Gradient): the gradient at the start of the first stage.
        let settings = cfg
            .inversion
            .clone()
            .context("the configuration has no [inversion] section")?
            .gradient_only();
        let stages = settings.stages()?;
        ensure!(
            truth.sources.len() == cfg.sources.len(),
            "the true model's survey has {} source(s), the start's {}",
            truth.sources.len(),
            cfg.sources.len()
        );
        let grid = Simulation::grid(cfg)?;
        ensure!(
            Simulation::grid(truth)? == grid,
            "the true model must be on the start's grid"
        );
        let start = Arc::new(Model::from_config(grid, &cfg.model, base)?);
        let true_model =
            Model::from_config(grid, &truth.model, base).context("in the true model")?;
        Ok(Self {
            cfg: cfg.clone(),
            truth: truth.clone(),
            base: base.to_path_buf(),
            observed: Vec::with_capacity(cfg.sources.len()),
            inversion: None,
            stage: Stage::Observing(0),
            stages,
            iterations: 0,
            in_stage: 0,
            stopped: None,
            first_misfit: None,
            pictures: Pictures {
                model: start.clone(),
                start,
                truth: Some(Arc::new(true_model)),
                gradient: [None, None, None],
            },
            settings,
        })
    }

    pub fn stage(&self) -> Stage {
        self.stage
    }

    /// Whether the inversion is set up: its starting misfit and gradient computed.
    pub fn started(&self) -> bool {
        self.inversion.is_some()
    }

    /// Units of work in all: a forward run per shot, the start, every stage's iterations
    /// and the start of every stage after the first.
    pub fn units(&self) -> usize {
        self.cfg.sources.len() + 1 + self.planned_iterations() + self.stages.len() - 1
    }

    /// The iterations the stages ask for, in all.
    fn planned_iterations(&self) -> usize {
        self.stages.iter().map(|s| s.iterations).sum()
    }

    /// The stage under way (from 0).
    fn inversion_stage(&self) -> usize {
        self.inversion.as_ref().map_or(0, Inversion::stage)
    }

    /// Units done.
    pub fn done(&self) -> usize {
        let shots = self.cfg.sources.len();
        match self.stage {
            Stage::Observing(s) => s,
            Stage::Starting => shots,
            Stage::Iterating(_) | Stage::Staging(_) => {
                shots + 1 + self.iterations + self.inversion_stage()
            }
            Stage::Done => self.units(),
        }
    }

    /// What follows the end of the stage under way: the next stage's start, or the end.
    fn after_stage(&self) -> Stage {
        let next = self.inversion_stage() + 1;
        if next < self.stages.len() {
            Stage::Staging(next)
        } else {
            Stage::Done
        }
    }

    /// What follows the start of stage `s`: its first iteration, or, if it has none, what
    /// follows it.
    fn entering_stage(&self, s: usize) -> Stage {
        if self.stages[s].iterations > 0 {
            Stage::Iterating(self.iterations + 1)
        } else {
            self.after_stage()
        }
    }

    /// Does the next unit of work; false once there is none left.
    pub fn advance(&mut self) -> Result<bool> {
        let shots = self.cfg.sources.len();
        self.stage = match self.stage {
            Stage::Observing(s) => {
                let mut sim = Simulation::new(&self.truth.shot(s)?, &self.base)
                    .context("setting up the true model")?;
                sim.run(None, false)?;
                let r = &sim.receivers;
                let shape = vec![r.points.len(), r.samples];
                self.observed.push(Observed {
                    shapes: Some(vec![shape; r.components.len()]),
                    times: Some(sim.sample_times().iter().map(|&t| t as f32).collect()),
                    traces: std::mem::take(&mut sim.receivers.traces),
                });
                if s + 1 < shots {
                    Stage::Observing(s + 1)
                } else {
                    Stage::Starting
                }
            }
            Stage::Starting => {
                let observed = std::mem::take(&mut self.observed);
                let inversion = Inversion::new(
                    &self.cfg,
                    (*self.pictures.start).clone(),
                    observed,
                    &self.settings,
                )?;
                self.first_misfit = Some(inversion.misfit());
                self.inversion = Some(inversion);
                self.picture();
                self.entering_stage(0)
            }
            Stage::Iterating(_) => {
                let inversion = self.inversion.as_mut().expect("set up at the start");
                let s = inversion.stage();
                let stepped = inversion.iterate()?;
                if stepped.is_some() {
                    self.iterations += 1;
                    self.in_stage += 1;
                } else {
                    self.stopped = Some(if self.stages.len() > 1 {
                        format!("no step lowers stage {}'s misfit further", s + 1)
                    } else {
                        "no step lowers the misfit further".into()
                    });
                }
                self.picture();
                if stepped.is_some() && self.in_stage < self.stages[s].iterations {
                    Stage::Iterating(self.iterations + 1)
                } else {
                    self.after_stage()
                }
            }
            Stage::Staging(s) => {
                let inversion = self.inversion.as_mut().expect("set up at the start");
                inversion.begin_stage(s)?;
                self.in_stage = 0;
                self.picture();
                self.entering_stage(s)
            }
            Stage::Done => Stage::Done,
        };
        Ok(self.stage != Stage::Done)
    }

    /// Takes the pictures of the inversion as it now stands.
    fn picture(&mut self) {
        let Some(inversion) = &self.inversion else {
            return;
        };
        // Freed first, so the old and the new are never held together.
        self.pictures.gradient = [None, None, None];
        for p in Parameter::ALL {
            self.pictures.gradient[p as usize] = inversion.gradient_of(p).map(Arc::new);
        }
        if inversion.iteration > 0 {
            self.pictures.model = Arc::new(inversion.model.clone());
        }
    }

    /// Where it is, as JSON: the stage (`observing` with its `shot`, `starting`, `iterating`
    /// with its `iteration`, `done`), the units done of all, the misfit at the start and now,
    /// every iteration (`inversion::Iteration`), why it stopped early if it did, and the most
    /// memory it has held.
    pub fn status(&self) -> serde_json::Value {
        let (stage, shot, iteration) = match self.stage {
            Stage::Observing(s) => ("observing", Some(s), None),
            Stage::Starting => ("starting", None, None),
            Stage::Iterating(i) => ("iterating", None, Some(i)),
            Stage::Staging(_) => ("staging", None, None),
            Stage::Done => ("done", None, None),
        };
        let inv = self.inversion.as_ref();
        json!({
            "stage": stage,
            "shot": shot,
            "iteration": iteration,
            "shots": self.cfg.sources.len(),
            "iterations": self.planned_iterations(),
            "stages": self.stages,
            "inversion_stage": match self.stage {
                Stage::Staging(s) => s,
                _ => self.inversion_stage(),
            },
            "stage_misfits": inv.map_or(json!([]), |i| json!(i.stage_misfits)),
            "parameters": self.settings.parameters,
            "units": self.units(),
            "done": self.done(),
            "unknowns": inv.map(Inversion::unknowns),
            "starting_misfit": self.first_misfit,
            "misfit": inv.map(Inversion::misfit),
            "history": inv.map_or(json!([]), |i| json!(i.history)),
            "stopped": self.stopped,
            "memory_bytes": self.memory_bytes(),
        })
    }

    /// The most memory held at once so far: the inversion's ([`Inversion::memory_bytes`])
    /// with the pictures.
    pub fn memory_bytes(&self) -> usize {
        self.inversion.as_ref().map_or(0, |i| i.memory_bytes) + self.pictures.memory_bytes()
    }

    /// What `invert` writes after each iteration: the model as `vp.npy`, `vs.npy` and
    /// `rho.npy`, and `inversion.json` with the settings and every iteration's misfit.
    pub fn write_outputs(&self, out: &Path) -> Result<()> {
        std::fs::create_dir_all(out).with_context(|| format!("creating {}", out.display()))?;
        let m = &self.pictures.model;
        let shape = [m.grid.n[2], m.grid.n[1], m.grid.n[0]];
        for (name, values) in [("vp", &m.vp), ("vs", &m.vs), ("rho", &m.rho)] {
            npy::write(&out.join(format!("{name}.npy")), &shape, values)?;
        }
        let inv = self.inversion.as_ref();
        let summary = json!({
            "settings": self.settings,
            "unknowns": inv.map(Inversion::unknowns),
            "starting_misfit": self.first_misfit,
            "iterations": inv.map_or(json!([]), |i| json!(i.history)),
            "memory_bytes": inv.map(|i| i.memory_bytes),
            "stopped": self.stopped,
        });
        std::fs::write(
            out.join("inversion.json"),
            serde_json::to_string_pretty(&summary)? + "\n",
        )?;
        Ok(())
    }
}

/// Bilinear interpolation of a row-major w × h grid to (w - 1)·f + 1 × (h - 1)·f + 1, which
/// keeps the original values at every f-th point. Interpolating values rather than colours
/// keeps zero crossings smooth.
fn upsample(v: &[f32], w: usize, h: usize, f: usize) -> Vec<f32> {
    if f == 1 {
        return v.to_vec();
    }
    let (fw, fh) = ((w - 1) * f + 1, (h - 1) * f + 1);
    // The cell a fine index falls in, and how far across it; the last index ends the last cell.
    let cell = |i: usize, n: usize| {
        let c = (i / f).min(n - 2);
        (c, (i - c * f) as f32 / f as f32)
    };
    let mut out = Vec::with_capacity(fw * fh);
    for y in 0..fh {
        let (r, ty) = cell(y, h);
        for x in 0..fw {
            let (c, tx) = cell(x, w);
            let at = |r: usize, c: usize| v[r * w + c];
            let top = at(r, c) * (1.0 - tx) + at(r, c + 1) * tx;
            let bottom = at(r + 1, c) * (1.0 - tx) + at(r + 1, c + 1) * tx;
            out.push(top * (1.0 - ty) + bottom * ty);
        }
    }
    out
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::config::Config;

    pub(crate) const CONFIG: &str = r#"
        [grid]
        points = [24, 20, 16]
        spacing = 10.0
        [time]
        duration = 0.03
        [boundary]
        pml_width = 4
        [[model.layer]]
        top = 0.0
        vp = 2000.0
        vs = 1000.0
        rho = 2000.0
        [[source]]
        position = [120.0, 100.0, 50.0]
        explosion = 1e12
        wavelet = { kind = "ricker", frequency = 20.0 }
        [receivers]
        components = ["vx", "vz"]
        [[receivers.line]]
        from = [60.0, 100.0, 0.0]
        to = [180.0, 100.0, 0.0]
        count = 5
    "#;

    fn view() -> View {
        let cfg = Config::from_toml(CONFIG).unwrap();
        View::new(Simulation::new(&cfg, Path::new(".")).unwrap())
    }

    #[test]
    fn steps_and_stops_at_the_end() {
        let mut v = view();
        let total = v.sim.steps;
        assert_eq!(v.advance(5, None).unwrap(), 5);
        assert_eq!(v.advance(10_000, None).unwrap(), total - 5);
        assert_eq!(v.advance(1, None).unwrap(), 0);
        assert_eq!(v.sim.receivers.recorded(), v.sim.receivers.samples);
    }

    #[test]
    fn draws_sections_of_each_plane() {
        let mut v = view();
        v.advance(10, None).unwrap();
        // Scaled up by 600 / 24 = 25 → capped at 8, keeping the corner nodes.
        for (plane, index, w, h) in [("xz", 10, 24, 16), ("xy", 0, 24, 20), ("yz", 12, 20, 16)] {
            let f = v.section("vz", plane, index, 1.0).unwrap();
            let (w, h) = ((w - 1) * 8 + 1, (h - 1) * 8 + 1);
            assert_eq!((f.width, f.height), (w, h), "{plane}");
            assert_eq!(f.rgba.len(), w * h * 4);
            assert!(f.clip > 0.0);
            assert!(f.png().unwrap().starts_with(b"\x89PNG"));
        }
        assert!(v.section("vz", "xz", 20, 1.0).is_err(), "y has 20 nodes");
        assert!(v.section("vq", "xz", 0, 1.0).is_err());
        assert!(v.section("vz", "volume", 0, 1.0).is_err());
        // More gain, a smaller clip.
        let one = v.section("p", "xz", 10, 1.0).unwrap().clip;
        let four = v.section("p", "xz", 10, 4.0).unwrap().clip;
        assert!((one / four - 4.0).abs() < 1e-4);
    }

    #[test]
    fn upsampling_keeps_the_nodes_and_interpolates_between() {
        let v = [0.0, 2.0, 4.0, 6.0]; // 2 × 2
        let up = upsample(&v, 2, 2, 2);
        assert_eq!(up, vec![0.0, 1.0, 2.0, 2.0, 3.0, 4.0, 4.0, 5.0, 6.0]);
        assert_eq!(upsample(&v, 2, 2, 1), v.to_vec());
    }

    #[test]
    fn the_colour_scale_remembers_the_run() {
        let mut v = view();
        v.advance(12, None).unwrap();
        let busy = v.section("vz", "xz", 10, 1.0).unwrap().clip;
        v.sim.solver.fields.vz.iter_mut().for_each(|x| *x *= 1e-6);
        let quiet = v.section("vz", "xz", 10, 1.0).unwrap().clip;
        assert!(
            (quiet - CLIP_MEMORY * busy).abs() < 1e-6 * busy,
            "{quiet} vs {busy}"
        );
    }

    #[test]
    fn draws_gathers_and_exports_traces() {
        let mut v = view();
        v.advance(10, None).unwrap();
        let f = v.gather("vz", 7).unwrap();
        assert_eq!((f.width, f.height, f.rgba.len()), (5, 7, 5 * 7 * 4));
        assert!(v.gather("p", 7).is_err(), "p isn't recorded");
        let bytes = v.traces_npy("vx").unwrap();
        assert!(bytes.starts_with(b"\x93NUMPY"));
        assert_eq!(bytes.len(), 128 + 4 * 5 * v.sim.receivers.samples);
    }

    /// A start and a true model with a faster sphere: two shots, receivers on the surface,
    /// the top held fixed, sampled alike.
    fn inversion_pair(iterations: usize) -> (Config, Config) {
        let config = |sphere: &str, inversion: &str| {
            Config::from_toml(&format!(
                r#"
                [grid]
                points = [24, 20, 18]
                spacing = 10.0
                [time]
                duration = 0.15
                dt = 5e-4
                [boundary]
                pml_width = 4
                pml_vp_max = 2700.0
                [scheme]
                coefficients = "taylor"
                [[model.layer]]
                top = 0.0
                vp = 2400.0
                vs = 1300.0
                rho = 2100.0
                {sphere}
                [[source]]
                position = [80.0, 90.0, 30.0]
                force = [0.0, 0.0, 1e9]
                wavelet = {{ kind = "ricker", frequency = 15.0 }}
                [[source]]
                position = [150.0, 100.0, 30.0]
                force = [1e9, 0.0, 0.0]
                wavelet = {{ kind = "ricker", frequency = 15.0 }}
                [receivers]
                components = ["vx", "vz"]
                [[receivers.line]]
                from = [60.0, 70.0, 0.0]
                to = [170.0, 70.0, 0.0]
                count = 10
                [[receivers.line]]
                from = [60.0, 120.0, 0.0]
                to = [170.0, 120.0, 0.0]
                count = 10
                {inversion}
                "#
            ))
            .unwrap()
        };
        let sphere = "[[model.sphere]]\ncenter = [115.0, 95.0, 100.0]\nradius = 25.0\n\
                      vp = 2600.0\nvs = 1450.0\nrho = 2100.0";
        let inversion = format!(
            "[inversion]\nobserved = \"observed\"\niterations = {iterations}\nfixed_above = 40.0\n\
             memory = 1"
        );
        (config("", &inversion), config(sphere, ""))
    }

    /// With frequency continuation the inversion starts each later stage as a unit of its
    /// own and counts iterations over all of them; asked for none (the page's Gradient), it
    /// stops after the first stage's start.
    #[test]
    fn steps_an_inversion_in_stages() {
        let base = Path::new(".");
        let staged = |extra: &str| {
            let (mut start, truth) = inversion_pair(1);
            let inv = start.inversion.as_mut().unwrap();
            inv.iterations = None;
            inv.stages = toml::from_str::<StagesToml>(extra).unwrap().stage;
            (start, truth)
        };
        let (start, truth) = staged(
            "[[stage]]\nmax_frequency = 10.0\niterations = 1\n\
             [[stage]]\niterations = 1\nmisfit = \"traveltime\"",
        );
        let mut v = InversionView::new(&start, &truth, base).unwrap();
        // Two shots, the start, two iterations and the second stage's start.
        assert_eq!(v.units(), 6);
        let mut stages = vec![v.stage()];
        while v.advance().unwrap() {
            stages.push(v.stage());
            assert!(v.done() <= v.units());
        }
        stages.push(v.stage());
        assert_eq!(
            stages,
            [
                Stage::Observing(0),
                Stage::Observing(1),
                Stage::Starting,
                Stage::Iterating(1),
                Stage::Staging(1),
                Stage::Iterating(2),
                Stage::Done
            ]
        );
        let status = v.status();
        assert_eq!(status["stages"].as_array().unwrap().len(), 2);
        assert_eq!(status["stage_misfits"].as_array().unwrap().len(), 2);
        assert_eq!(status["history"][1]["stage"], 1);
        assert_eq!(status["done"], status["units"]);

        // No iterations: the gradient at the first stage's start.
        let (mut start, truth) =
            staged("[[stage]]\nmax_frequency = 10.0\niterations = 3\n[[stage]]\niterations = 2");
        start.inversion.as_mut().unwrap().iterations = Some(0);
        // The page plans the memory first, which must take it the same way (it refused the
        // iterations beside the stages once).
        crate::memory::estimate_synthetic_inversion(&start, base).unwrap();
        let mut v = InversionView::new(&start, &truth, base).unwrap();
        assert_eq!(v.units(), 3);
        while v.advance().unwrap() {}
        let status = v.status();
        assert_eq!(status["stages"][0]["max_frequency"], 10.0);
        assert!(status["history"].as_array().unwrap().is_empty());
        assert!(v.pictures.gradient.iter().any(Option::is_some));
    }

    /// `[[stage]]` tables, read on their own.
    #[derive(serde::Deserialize)]
    struct StagesToml {
        stage: Vec<StageConfig>,
    }

    /// The inversion goes through its stages a unit at a time, lowers the misfit, draws every
    /// picture, writes what `invert` writes, and holds what the estimate says.
    #[test]
    fn steps_an_inversion_of_synthetic_data() {
        let (start, truth) = inversion_pair(2);
        let base = Path::new(".");
        let mut v = InversionView::new(&start, &truth, base).unwrap();
        assert_eq!((v.units(), v.done()), (5, 0));
        assert!(v.pictures.section("vs", "descent", "xz", 9, 1.0).is_err());
        let mut stages = vec![v.stage()];
        while v.advance().unwrap() {
            stages.push(v.stage());
        }
        stages.push(v.stage());
        assert_eq!(
            stages,
            [
                Stage::Observing(0),
                Stage::Observing(1),
                Stage::Starting,
                Stage::Iterating(1),
                Stage::Iterating(2),
                Stage::Done
            ]
        );
        assert!(!v.advance().unwrap(), "nothing left to do");
        let s = v.status();
        assert_eq!(
            (s["stage"].as_str(), s["done"].as_u64()),
            (Some("done"), Some(5))
        );
        let first = s["starting_misfit"].as_f64().unwrap();
        let last = s["misfit"].as_f64().unwrap();
        assert_eq!(s["history"].as_array().unwrap().len(), 2);
        assert!(last < 0.8 * first, "{last} vs {first}");

        // Pictures: the sizes of a wavefield section, models on the start and truth's range.
        let model = v.pictures.section("vs", "model", "xz", 9, 1.0).unwrap();
        assert_eq!((model.width, model.height), (23 * 8 + 1, 17 * 8 + 1));
        assert_eq!((model.low, model.clip), (1300.0, 1450.0));
        let truth = v
            .pictures
            .section("vs", "true-update", "xz", 9, 2.0)
            .unwrap();
        assert_eq!(truth.clip, 75.0, "the largest true change over the gain");
        let update = v.pictures.section("vs", "update", "xz", 9, 2.0).unwrap();
        assert_eq!(update.clip, truth.clip, "on the same scale");
        let descent = v.pictures.section("vp", "descent", "xy", 10, 1.0).unwrap();
        assert!(descent.clip > 0.0 && descent.low == -descent.clip);
        assert!(
            v.pictures.section("rho", "descent", "xz", 9, 1.0).is_err(),
            "not inverted"
        );
        assert!(
            v.pictures.section("vs", "model", "xz", 20, 1.0).is_err(),
            "y has 20 nodes"
        );
        assert!(v.pictures.section("vq", "model", "xz", 0, 1.0).is_err());
        assert!(v.pictures.npy("vs").unwrap().starts_with(b"\x93NUMPY"));

        // The memory planned for it, exactly, with the pictures it holds by now (and L-BFGS's
        // one pair, as full as the plan counts it).
        let planned = crate::memory::estimate_synthetic_inversion(&start, base).unwrap();
        assert_eq!(v.memory_bytes() as u64, planned.peak);

        let dir = std::env::temp_dir().join(format!("vibe-rations-view-{}", std::process::id()));
        v.write_outputs(&dir).unwrap();
        for f in ["vp.npy", "vs.npy", "rho.npy", "inversion.json"] {
            assert!(dir.join(f).is_file(), "{f}");
        }
        std::fs::remove_dir_all(dir).unwrap();
    }

    /// Observed traces sampled at other times than the synthetics are refused at the start.
    #[test]
    fn refuses_a_truth_sampled_otherwise() {
        let (start, mut truth) = inversion_pair(0);
        truth.time.dt = Some(4e-4);
        let mut v = InversionView::new(&start, &truth, Path::new(".")).unwrap();
        v.advance().unwrap();
        v.advance().unwrap();
        let e = format!("{:#}", v.advance().unwrap_err());
        assert!(e.contains("sampl"), "{e}");
        let mut other = truth.clone();
        other.grid.points = [24, 20, 16];
        assert!(InversionView::new(&start, &other, Path::new(".")).is_err());
    }

    #[test]
    fn describes_the_run_with_receiver_positions() {
        let info = view().info();
        assert_eq!(info["points"], serde_json::json!([24, 20, 16]));
        assert_eq!(
            info["receiver_positions"][4],
            serde_json::json!([180.0, 100.0, 0.0])
        );
    }
}
