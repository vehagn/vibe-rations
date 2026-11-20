//! The solver in the browser. The web app's worker drives a [`Simulator`]: it sets a run up
//! from TOML, advances it a chunk of steps at a time, and asks for pictures of a section and
//! of the gather, as RGBA pixels, between chunks.
//!
//! The pictures come from [`vibe_rations::view::View`], which `vibe-rations serve` shares;
//! this crate only checks what a browser can't do and turns errors into JavaScript ones.

use std::path::Path;

use anyhow::{Result, ensure};
use vibe_rations::config::Config;
use vibe_rations::memory::{self, Memory};
use vibe_rations::simulation::Simulation;
use vibe_rations::view::{self, View};
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

    /// The value drawn at full colour.
    #[wasm_bindgen(getter)]
    pub fn clip(&self) -> f32 {
        self.0.clip
    }

    /// The pixels, RGBA, row by row from the top; hands them over without copying.
    pub fn take_rgba(&mut self) -> Vec<u8> {
        std::mem::take(&mut self.0.rgba)
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
    /// The memory a configuration would need, as JSON (`memory::Memory`: parts, resident,
    /// peak, plus the browser's limit), without setting it up.
    #[wasm_bindgen(js_name = estimateMemory)]
    pub fn estimate_memory(config: &str) -> Result<String, JsError> {
        let m = estimate(config).map_err(js)?;
        let mut value = serde_json::to_value(&m).expect("Memory serialises");
        value["limit"] = serde_json::json!(BROWSER_LIMIT);
        Ok(value.to_string())
    }

    /// Sets up a run from a TOML configuration.
    #[wasm_bindgen(constructor)]
    pub fn new(config: &str) -> Result<Simulator, JsError> {
        setup(config).map(Simulator).map_err(js)
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

/// What WebAssembly can address: 4 GiB. The page refuses runs that would need more, less a
/// margin for the program and the pictures, before setting them up.
pub const BROWSER_LIMIT: u64 = 4 << 30;
const MARGIN: u64 = 128 << 20;

/// The memory a configuration would need in the browser, part by part, without setting it up.
pub fn estimate(config: &str) -> Result<Memory> {
    let cfg = browser_config(config)?;
    memory::estimate(&cfg, Path::new("."))
}

fn browser_config(config: &str) -> Result<Config> {
    let cfg = Config::from_toml(config)?;
    ensure!(
        cfg.model.files.is_none(),
        "the browser can't read model files; use [[model.layer]]s and [[model.sphere]]s"
    );
    Ok(cfg)
}

/// Sets a run up for the browser: no model files, and no more memory than WebAssembly has.
pub fn setup(config: &str) -> Result<View> {
    let cfg = browser_config(config)?;
    let m = memory::estimate(&cfg, Path::new("."))?;
    ensure!(
        m.peak + MARGIN <= BROWSER_LIMIT,
        "this run needs {} of memory, and a browser gives WebAssembly at most {}; use a \
         smaller grid, or the command-line program",
        memory::human(m.peak),
        memory::human(BROWSER_LIMIT)
    );
    Ok(View::new(Simulation::new(&cfg, Path::new("."))?))
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
        let mut v = setup(CONFIG).unwrap();
        assert_eq!(v.advance(5, None).unwrap(), 5);
    }

    #[test]
    fn estimates_memory_and_refuses_runs_too_large_for_the_browser() {
        let m = estimate(CONFIG).unwrap();
        assert_eq!(m, setup(CONFIG).unwrap().sim.info.memory);
        let e = format!("{:#}", setup(&grown(CONFIG)).err().unwrap());
        assert!(e.contains("at most 4.00 GiB"), "{e}");
        assert!(estimate(&grown(CONFIG)).unwrap().peak > BROWSER_LIMIT);
    }

    /// The test configuration on a grid far too large for a browser.
    fn grown(text: &str) -> String {
        let start = text.find("points = [").unwrap();
        let end = start + text[start..].find(']').unwrap() + 1;
        format!("{}points = [800, 800, 800]{}", &text[..start], &text[end..])
    }

    #[test]
    fn refuses_model_files() {
        let text = CONFIG.replace(
            "[[model.layer]]\n        top = 0.0\n        vp = 2000.0\n        vs = 1000.0\n        rho = 2000.0",
            "[model.files]\nvp = \"a.npy\"\nvs = \"b.npy\"\nrho = \"c.npy\"",
        );
        let e = format!("{:#}", setup(&text).err().unwrap());
        assert!(e.contains("model files"), "{e}");
    }
}
