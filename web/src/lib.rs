//! The solver in the browser. The web app's worker drives a [`Simulator`]: it sets a run up
//! from TOML, advances it a chunk of steps at a time, and asks for pictures of a section and
//! of the gather, as RGBA pixels, between chunks.
//!
//! The pictures come from [`vibe_rations::view::View`], which `vibe-rations serve` shares;
//! this crate only checks what a browser can't do and turns errors into JavaScript ones. A
//! browser has no files, so the page fetches the model files a configuration names and hands
//! them over as [`ModelFiles`].

use std::path::Path;

use anyhow::{Context, Result, ensure};
use vibe_rations::config::Config;
use vibe_rations::memory::{self, Memory};
use vibe_rations::model::Files;
use vibe_rations::simulation::Simulation;
use vibe_rations::view::{self, InversionView, View};
use wasm_bindgen::prelude::*;

pub mod settings;

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
pub struct Frame(view::Frame);

#[wasm_bindgen]
impl Frame {
    #[wasm_bindgen(getter)]
    pub fn width(&self) -> usize {
        self.0.width
    }

    #[wasm_bindgen(getter)]
    pub fn height(&self) -> usize {
        self.0.height
    }

    /// The value drawn at full colour: the top of the colour scale.
    #[wasm_bindgen(getter)]
    pub fn clip(&self) -> f32 {
        self.0.clip
    }

    /// The bottom of the colour scale: -clip for signed pictures.
    #[wasm_bindgen(getter)]
    pub fn low(&self) -> f32 {
        self.0.low
    }

    /// The pixels, RGBA, row by row from the top; hands them over without copying.
    pub fn take_rgba(&mut self) -> Vec<u8> {
        std::mem::take(&mut self.0.rgba)
    }
}

/// The model files a configuration names (`Simulator.modelFiles`), by those names, as the
/// page fetched them.
#[wasm_bindgen]
#[derive(Default)]
pub struct ModelFiles(Vec<(String, Vec<u8>)>);

#[wasm_bindgen]
impl ModelFiles {
    #[wasm_bindgen(constructor)]
    pub fn new() -> ModelFiles {
        ModelFiles::default()
    }

    /// Adds the `.npy` file the configuration calls `name`.
    pub fn add(&mut self, name: String, bytes: Vec<u8>) {
        self.0.push((name, bytes));
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
    /// The model files a configuration names, as it names them (relative to it); none if it
    /// has none, or doesn't parse (setting it up then says why).
    #[wasm_bindgen(js_name = modelFiles)]
    pub fn model_files(config: &str) -> Vec<String> {
        model_files(config)
    }

    /// The memory a configuration would need, as JSON (`memory::Memory`: parts, resident,
    /// peak, plus the browser's limit), without setting it up; `files` holds its model files.
    #[wasm_bindgen(js_name = estimateMemory)]
    pub fn estimate_memory(config: &str, files: &ModelFiles) -> Result<String, JsError> {
        let m = estimate(config, &files.0).map_err(js)?;
        let mut value = serde_json::to_value(&m).expect("Memory serialises");
        value["limit"] = serde_json::json!(BROWSER_LIMIT);
        Ok(value.to_string())
    }

    /// Sets up a run from a TOML configuration and the model files it names.
    #[wasm_bindgen(constructor)]
    pub fn new(config: &str, files: &ModelFiles) -> Result<Simulator, JsError> {
        setup(config, &files.0).map(Simulator).map_err(js)
    }

    /// The run's settings as JSON: `run.json`'s fields, plus the receivers' positions.
    pub fn info(&self) -> String {
        self.0.info().to_string()
    }

    /// Takes up to `steps` steps, recording the receivers; returns how many it took. Fails if
    /// the wavefield has blown up.
    pub fn advance(&mut self, steps: usize) -> Result<usize, JsError> {
        self.0.advance(steps, None).map_err(js)
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
        self.0
            .section(field, plane, index, gain)
            .map(Frame)
            .map_err(js)
    }

    /// The gather of `component` so far, one column per receiver and `rows` rows of time.
    pub fn gather(&self, component: &str, rows: usize) -> Result<Frame, JsError> {
        self.0.gather(component, rows).map(Frame).map_err(js)
    }

    /// The traces of `component` as an `.npy` file, (receivers, samples).
    #[wasm_bindgen(js_name = tracesNpy)]
    pub fn traces_npy(&self, component: &str) -> Result<Vec<u8>, JsError> {
        self.0.traces_npy(component).map_err(js)
    }
}

/// An inversion of synthetic data the page can step a unit at a time and look at.
#[wasm_bindgen]
pub struct Inverter(InversionView);

#[wasm_bindgen]
impl Inverter {
    /// The memory an inversion of the start `config` (with `[inversion]`) would need, as
    /// JSON like [`Simulator::estimate_memory`]'s.
    #[wasm_bindgen(js_name = estimateMemory)]
    pub fn estimate_memory(config: &str) -> Result<String, JsError> {
        let m = estimate_inversion(config).map_err(js)?;
        let mut value = serde_json::to_value(&m).expect("Memory serialises");
        value["limit"] = serde_json::json!(BROWSER_LIMIT);
        Ok(value.to_string())
    }

    /// An inversion of `start`'s model (with its `[inversion]`) against the traces of
    /// `truth` (the same survey with the true model), as [`settings::synthetic`] makes them.
    #[wasm_bindgen(constructor)]
    pub fn new(start: &str, truth: &str) -> Result<Inverter, JsError> {
        setup_inversion(start, truth).map(Inverter).map_err(js)
    }

    /// Does the next unit of work: a shot of observed traces, the start, an iteration. False
    /// once there is none left.
    pub fn advance(&mut self) -> Result<bool, JsError> {
        self.0.advance().map_err(js)
    }

    /// Where it is, as JSON (`InversionView::status`).
    pub fn status(&self) -> String {
        self.0.status().to_string()
    }

    /// A picture of `parameter` (vp, vs, rho) showing `show` (model, update, true-update,
    /// descent) on `plane` through node `index`.
    pub fn picture(
        &self,
        parameter: &str,
        show: &str,
        plane: &str,
        index: usize,
        gain: f32,
    ) -> Result<Frame, JsError> {
        self.0
            .pictures
            .section(parameter, show, plane, index, gain)
            .map(Frame)
            .map_err(js)
    }

    /// The model as it stands, `parameter` as an `.npy` file (nz, ny, nx).
    #[wasm_bindgen(js_name = modelNpy)]
    pub fn model_npy(&self, parameter: &str) -> Result<Vec<u8>, JsError> {
        self.0.pictures.npy(parameter).map_err(js)
    }
}

/// What the configuration `text` sets, defaults filled in, as JSON ([`settings::settings`]).
#[wasm_bindgen(js_name = readSettings)]
pub fn read_settings(text: &str) -> Result<String, JsError> {
    settings::settings(text).map(|v| v.to_string()).map_err(js)
}

/// The configuration `text` with `edits` (JSON: `[{path, value}]`) applied
/// ([`settings::edit`]).
#[wasm_bindgen(js_name = editConfig)]
pub fn edit_config(text: &str, edits: &str) -> Result<String, JsError> {
    let edits: serde_json::Value = serde_json::from_str(edits)
        .context("edits must be JSON")
        .map_err(js)?;
    settings::edit(text, &edits).map_err(js)
}

/// The commands that run a configuration without the page, as JSON ([`settings::commands`]).
#[wasm_bindgen]
pub fn commands(input: &str) -> Result<String, JsError> {
    let input: serde_json::Value = serde_json::from_str(input)
        .context("the input must be JSON")
        .map_err(js)?;
    settings::commands(&input)
        .map(|v| v.to_string())
        .map_err(js)
}

/// A true model for a synthetic inversion of `text`: its model with a sphere added.
#[wasm_bindgen(js_name = defaultTruth)]
pub fn default_truth(text: &str) -> Result<String, JsError> {
    settings::default_truth(text).map_err(js)
}

/// The `[model]` of a configuration alone, as TOML ([`settings::model_of`]).
#[wasm_bindgen(js_name = modelOf)]
pub fn model_of(text: &str) -> Result<String, JsError> {
    settings::model_of(text).map_err(js)
}

/// The start and the truth of a synthetic inversion, as JSON `{start, truth}`
/// ([`settings::synthetic`]).
#[wasm_bindgen]
pub fn synthetic(text: &str, truth: &str, observed: &str) -> Result<String, JsError> {
    let (start, truth) = settings::synthetic(text, truth, observed).map_err(js)?;
    Ok(serde_json::json!({ "start": start, "truth": truth }).to_string())
}

/// What WebAssembly can address: 4 GiB. The page refuses runs that would need more, less a
/// margin for the program and the pictures, before setting them up.
pub const BROWSER_LIMIT: u64 = 4 << 30;
const MARGIN: u64 = 128 << 20;

/// The model files `config` names, in the order vp, vs, rho, each once.
pub fn model_files(config: &str) -> Vec<String> {
    let Ok(cfg) = Config::from_toml(config) else {
        return vec![];
    };
    let mut names = Vec::new();
    for f in cfg.model.file_paths() {
        let name = f.to_string_lossy().into_owned();
        if !names.contains(&name) {
            names.push(name);
        }
    }
    names
}

/// The memory a configuration would need in the browser, part by part, without setting it
/// up; `files` holds its model files by name.
pub fn estimate(config: &str, files: &[(String, Vec<u8>)]) -> Result<Memory> {
    memory::estimate_with(&Config::from_toml(config)?, Files::Memory(files))
}

/// A configuration for an inversion in the browser. Its models come from layers and spheres:
/// `InversionView` reads model files from a directory, which a browser hasn't got.
fn inversion_config(config: &str) -> Result<Config> {
    let cfg = Config::from_toml(config)?;
    ensure!(
        cfg.model.files.is_none(),
        "the browser's inversions can't read model files; use [[model.layer]]s and \
         [[model.sphere]]s, or the run server"
    );
    Ok(cfg)
}

/// The memory an inversion of synthetic data from the start `config` would need in the
/// browser.
pub fn estimate_inversion(config: &str) -> Result<Memory> {
    let cfg = inversion_config(config)?;
    memory::estimate_synthetic_inversion(&cfg, Path::new("."))
}

/// Sets up an inversion of synthetic data for the browser, if it fits.
pub fn setup_inversion(start: &str, truth: &str) -> Result<InversionView> {
    let cfg = inversion_config(start)?;
    let truth = inversion_config(truth).context("in the true model")?;
    let m = memory::estimate_synthetic_inversion(&cfg, Path::new("."))?;
    ensure!(
        m.peak + MARGIN <= BROWSER_LIMIT,
        "this inversion needs {} of memory, and a browser gives WebAssembly at most {}; use a \
         smaller grid, or the run server",
        memory::human(m.peak),
        memory::human(BROWSER_LIMIT)
    );
    InversionView::new(&cfg, &truth, Path::new("."))
}

/// Sets a run up for the browser, its model files from `files`: no more memory than
/// WebAssembly has.
pub fn setup(config: &str, files: &[(String, Vec<u8>)]) -> Result<View> {
    let cfg = Config::from_toml(config)?;
    let m = memory::estimate_with(&cfg, Files::Memory(files))?;
    ensure!(
        m.peak + MARGIN <= BROWSER_LIMIT,
        "this run needs {} of memory, and a browser gives WebAssembly at most {}; use a \
         smaller grid, or the command-line program",
        memory::human(m.peak),
        memory::human(BROWSER_LIMIT)
    );
    Ok(View::new(Simulation::with_files(
        &cfg,
        Files::Memory(files),
    )?))
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
    fn sets_up_and_steps() {
        let mut v = setup(CONFIG, &[]).unwrap();
        assert_eq!(v.advance(5, None).unwrap(), 5);
    }

    #[test]
    fn estimates_memory_and_refuses_runs_too_large_for_the_browser() {
        let m = estimate(CONFIG, &[]).unwrap();
        assert_eq!(m, setup(CONFIG, &[]).unwrap().sim.info.memory);
        let e = format!("{:#}", setup(&grown(CONFIG), &[]).err().unwrap());
        assert!(e.contains("at most 4.00 GiB"), "{e}");
        assert!(estimate(&grown(CONFIG), &[]).unwrap().peak > BROWSER_LIMIT);
    }

    /// The test configuration on a grid far too large for a browser.
    fn grown(text: &str) -> String {
        let start = text.find("points = [").unwrap();
        let end = start + text[start..].find(']').unwrap() + 1;
        format!("{}points = [800, 800, 800]{}", &text[..start], &text[end..])
    }

    /// The page's inversion: the settings module's pair of configurations, stepped to the end;
    /// one too large for the browser is refused before anything is allocated.
    #[test]
    fn inverts_synthetic_data_and_refuses_inversions_too_large() {
        let with = settings::edit(
            CONFIG,
            &serde_json::json!([
                { "path": "inversion.observed", "value": "observed" },
                { "path": "inversion.iterations", "value": 1 },
                { "path": "inversion.parameters", "value": ["vs"] },
                { "path": "scheme.coefficients", "value": "taylor" },
            ]),
        )
        .unwrap();
        let truth = settings::default_truth(&with).unwrap();
        let (start, truth) = settings::synthetic(&with, &truth, "observed").unwrap();
        let m = estimate_inversion(&start).unwrap();
        assert!(m.peak > estimate(&start, &[]).unwrap().peak);
        let mut inv = setup_inversion(&start, &truth).unwrap();
        while inv.advance().unwrap() {}
        let s = inv.status();
        assert_eq!(s["stage"], "done");
        assert!(s["misfit"].as_f64().unwrap() < s["starting_misfit"].as_f64().unwrap());
        let f = inv.pictures.section("vs", "update", "xz", 10, 1.0).unwrap();
        assert_eq!(f.rgba.len(), f.width * f.height * 4);
        let e = format!(
            "{:#}",
            setup_inversion(&grown(&start), &grown(&truth))
                .err()
                .unwrap()
        );
        assert!(e.contains("at most 4.00 GiB"), "{e}");
    }

    #[test]
    fn runs_model_files_it_is_given() {
        let text = CONFIG.replace(
            "[[model.layer]]\n        top = 0.0\n        vp = 2000.0\n        vs = 1000.0\n        rho = 2000.0",
            "[model.files]\nvp = \"m/vp.npy\"\nvs = \"m/vs.npy\"\nrho = \"m/rho.npy\"",
        );
        assert_eq!(model_files(&text), ["m/vp.npy", "m/vs.npy", "m/rho.npy"]);
        assert!(model_files(CONFIG).is_empty() && model_files("not toml").is_empty());
        // Sections, (nz, nx), the same as the layer at every y.
        let section = |v: f32| vibe_rations::npy::to_bytes(&[16, 24], &[v; 16 * 24]).unwrap();
        let files = vec![
            ("m/vp.npy".to_string(), section(2000.0)),
            ("m/vs.npy".to_string(), section(1000.0)),
            ("m/rho.npy".to_string(), section(2000.0)),
        ];
        let (mut from_files, mut layered) =
            (setup(&text, &files).unwrap(), setup(CONFIG, &[]).unwrap());
        assert_eq!(estimate(&text, &files).unwrap(), from_files.sim.info.memory);
        from_files.advance(5, None).unwrap();
        layered.advance(5, None).unwrap();
        assert_eq!(
            from_files.traces_npy("vz").unwrap(),
            layered.traces_npy("vz").unwrap()
        );
        let e = format!("{:#}", setup(&text, &files[..2]).err().unwrap());
        assert!(e.contains("m/rho.npy wasn't given"), "{e}");
    }
}
