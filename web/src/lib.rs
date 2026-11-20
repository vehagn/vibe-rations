//! The solver in the browser. The web app's worker drives a [`Simulator`]: it sets a run up
//! from TOML, advances it a chunk of steps at a time, and asks for pictures of a section and
//! of the gather, as RGBA pixels, between chunks.
//!
//! [`View`] does the work in plain Rust, so it can be tested natively; [`Simulator`] only
//! turns its errors into JavaScript ones.

use std::path::Path;

use anyhow::{Context, Result, bail, ensure};
use vibe_rations::config::{Config, PlaneKind};
use vibe_rations::image;
use vibe_rations::npy;
use vibe_rations::receiver::Component;
use vibe_rations::simulation::Simulation;
use vibe_rations::snapshot::Snapshot;
use wasm_bindgen::prelude::*;

#[wasm_bindgen]
extern "C" {
    #[wasm_bindgen(js_namespace = console)]
    fn error(message: &str);
}

/// Starts rayon's thread pool on Web Workers: `await initThreadPool(navigator.hardwareConcurrency)`
/// once before running anything.
#[cfg(feature = "threads")]
pub use wasm_bindgen_rayon::init_thread_pool;

/// Sends panics to the console, which otherwise shows only "unreachable".
#[wasm_bindgen(start)]
pub fn start() {
    std::panic::set_hook(Box::new(|info| error(&info.to_string())));
}

/// Pixels to draw.
#[wasm_bindgen]
pub struct Frame {
    width: usize,
    height: usize,
    rgba: Vec<u8>,
    clip: f32,
}

#[wasm_bindgen]
impl Frame {
    #[wasm_bindgen(getter)]
    pub fn width(&self) -> usize {
        self.width
    }

    #[wasm_bindgen(getter)]
    pub fn height(&self) -> usize {
        self.height
    }

    /// The value drawn at full colour.
    #[wasm_bindgen(getter)]
    pub fn clip(&self) -> f32 {
        self.clip
    }

    /// The pixels, RGBA, row by row from the top; hands them over without copying.
    pub fn take_rgba(&mut self) -> Vec<u8> {
        std::mem::take(&mut self.rgba)
    }
}

/// A simulation the page can step and look at.
#[wasm_bindgen]
pub struct Simulator(View);

fn js(e: anyhow::Error) -> JsError {
    JsError::new(&format!("{e:#}"))
}

#[wasm_bindgen]
impl Simulator {
    /// Sets up a run from a TOML configuration.
    #[wasm_bindgen(constructor)]
    pub fn new(config: &str) -> Result<Simulator, JsError> {
        View::new(config).map(Simulator).map_err(js)
    }

    /// The run's settings as JSON: `run.json`'s fields, plus the receivers' positions.
    pub fn info(&self) -> String {
        self.0.info()
    }

    /// Takes up to `steps` steps, recording the receivers; returns how many it took. Fails if
    /// the wavefield has blown up.
    pub fn advance(&mut self, steps: usize) -> Result<usize, JsError> {
        self.0.advance(steps).map_err(js)
    }

    /// Steps taken.
    #[wasm_bindgen(getter)]
    pub fn step(&self) -> usize {
        self.0.sim.solver.step
    }

    /// Steps the run takes in all.
    #[wasm_bindgen(getter, js_name = totalSteps)]
    pub fn total_steps(&self) -> usize {
        self.0.sim.steps
    }

    /// Simulated time, in seconds.
    #[wasm_bindgen(getter)]
    pub fn time(&self) -> f64 {
        self.0.sim.solver.step as f64 * self.0.sim.solver.dt
    }

    /// Interior grid nodes, for working out speeds.
    #[wasm_bindgen(getter)]
    pub fn cells(&self) -> usize {
        self.0.sim.solver.grid.interior_len()
    }

    /// `field` (vx, vy, vz or p) on `plane` (xy, xz or yz) through node `index` along its
    /// normal, scaled so `gain` 1 shows the 99.5th percentile at full colour.
    pub fn section(
        &mut self,
        field: &str,
        plane: &str,
        index: usize,
        gain: f32,
    ) -> Result<Frame, JsError> {
        self.0.section(field, plane, index, gain).map_err(js)
    }

    /// The gather of `component` so far, one column per receiver and `rows` rows of time.
    pub fn gather(&self, component: &str, rows: usize) -> Result<Frame, JsError> {
        self.0.gather(component, rows).map_err(js)
    }

    /// The traces of `component` as an `.npy` file, (receivers, samples).
    #[wasm_bindgen(js_name = tracesNpy)]
    pub fn traces_npy(&self, component: &str) -> Result<Vec<u8>, JsError> {
        self.0.traces_npy(component).map_err(js)
    }
}

/// The simulation and the views of it.
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

fn component(name: &str) -> Result<Component> {
    Ok(match name {
        "vx" => Component::Vx,
        "vy" => Component::Vy,
        "vz" => Component::Vz,
        "p" => Component::P,
        _ => bail!("unknown field {name}; use vx, vy, vz or p"),
    })
}

impl View {
    pub fn new(config: &str) -> Result<Self> {
        let cfg = Config::from_toml(config)?;
        ensure!(
            cfg.model.files.is_none(),
            "the browser can't read model files; use [[model.layer]]s and [[model.sphere]]s"
        );
        Ok(Self {
            sim: Simulation::new(&cfg, Path::new("."))?,
            largest_clip: [0.0; 4],
        })
    }

    pub fn info(&self) -> String {
        let mut info = serde_json::to_value(&self.sim.info).expect("RunInfo serialises");
        let g = self.sim.solver.grid;
        let receivers: Vec<[f64; 3]> = self
            .sim
            .receivers
            .nodes
            .iter()
            .map(|&n| g.position(n))
            .collect();
        info["receiver_positions"] = serde_json::json!(receivers);
        info.to_string()
    }

    pub fn advance(&mut self, steps: usize) -> Result<usize> {
        let n = steps.min(self.sim.remaining());
        for _ in 0..n {
            self.sim.advance();
        }
        ensure!(
            self.sim.solver.peak_velocity().is_finite(),
            "the wavefield blew up by step {}; check the model, or lower dt or courant",
            self.sim.solver.step
        );
        Ok(n)
    }

    pub fn section(&mut self, field: &str, plane: &str, index: usize, gain: f32) -> Result<Frame> {
        let (plane, axis) = match plane {
            "xy" => (PlaneKind::Xy, 2),
            "xz" => (PlaneKind::Xz, 1),
            "yz" => (PlaneKind::Yz, 0),
            _ => bail!("unknown plane {plane}; use xy, xz or yz"),
        };
        let n = self.sim.solver.grid.n[axis];
        ensure!(index < n, "plane index {index} is off the grid (0..{n})");
        ensure!(gain > 0.0, "gain must be positive");
        let component = component(field)?;
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
        })
    }

    pub fn gather(&self, name: &str, rows: usize) -> Result<Frame> {
        let r = &self.sim.receivers;
        let c = self.recorded(name)?;
        let rows = rows.clamp(1, r.samples.max(1));
        Ok(Frame {
            width: r.nodes.len(),
            height: rows,
            rgba: image::gather_rgba(r.nodes.len(), r.samples, r.recorded(), &r.traces[c], rows),
            clip: 1.0,
        })
    }

    pub fn traces_npy(&self, name: &str) -> Result<Vec<u8>> {
        let r = &self.sim.receivers;
        let i = self.recorded(name)?;
        npy::to_bytes(&[r.nodes.len(), r.samples], &r.traces[i])
    }

    /// Which of the receivers' trace arrays holds `name`.
    fn recorded(&self, name: &str) -> Result<usize> {
        let c = component(name)?;
        self.sim
            .receivers
            .components
            .iter()
            .position(|&x| x == c)
            .with_context(|| format!("the receivers don't record {name}"))
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

    const CONFIG: &str = r#"
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

    #[test]
    fn steps_and_stops_at_the_end() {
        let mut v = View::new(CONFIG).unwrap();
        let total = v.sim.steps;
        assert_eq!(v.advance(5).unwrap(), 5);
        assert_eq!(v.advance(10_000).unwrap(), total - 5);
        assert_eq!(v.advance(1).unwrap(), 0);
        assert_eq!(v.sim.receivers.recorded(), v.sim.receivers.samples);
    }

    #[test]
    fn draws_sections_of_each_plane() {
        let mut v = View::new(CONFIG).unwrap();
        v.advance(10).unwrap();
        // Scaled up by 600 / 24 = 25 → capped at 8, keeping the corner nodes.
        for (plane, index, w, h) in [("xz", 10, 24, 16), ("xy", 0, 24, 20), ("yz", 12, 20, 16)] {
            let mut f = v.section("vz", plane, index, 1.0).unwrap();
            let (w, h) = ((w - 1) * 8 + 1, (h - 1) * 8 + 1);
            assert_eq!((f.width, f.height), (w, h), "{plane}");
            assert_eq!(f.take_rgba().len(), w * h * 4);
            assert!(f.clip > 0.0);
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
        let mut v = View::new(CONFIG).unwrap();
        v.advance(12).unwrap();
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
        let mut v = View::new(CONFIG).unwrap();
        v.advance(10).unwrap();
        let f = v.gather("vz", 7).unwrap();
        assert_eq!((f.width, f.height, f.rgba.len()), (5, 7, 5 * 7 * 4));
        assert!(v.gather("p", 7).is_err(), "p isn't recorded");
        let bytes = v.traces_npy("vx").unwrap();
        assert!(bytes.starts_with(b"\x93NUMPY"));
        assert_eq!(bytes.len(), 128 + 4 * 5 * v.sim.receivers.samples);
    }

    #[test]
    fn describes_the_run_with_receiver_positions() {
        let v = View::new(CONFIG).unwrap();
        let info: serde_json::Value = serde_json::from_str(&v.info()).unwrap();
        assert_eq!(info["points"], serde_json::json!([24, 20, 16]));
        assert_eq!(
            info["receiver_positions"][4],
            serde_json::json!([180.0, 100.0, 0.0])
        );
    }

    #[test]
    fn refuses_model_files() {
        let text = CONFIG.replace(
            "[[model.layer]]\n        top = 0.0\n        vp = 2000.0\n        vs = 1000.0\n        rho = 2000.0",
            "[model.files]\nvp = \"a.npy\"\nvs = \"b.npy\"\nrho = \"c.npy\"",
        );
        let e = format!("{:#}", View::new(&text).err().unwrap());
        assert!(e.contains("model files"), "{e}");
    }
}
