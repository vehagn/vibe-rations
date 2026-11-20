//! Live pictures of a run, for the web app: a section of the wavefield and the gather so far,
//! as RGBA pixels (the browser draws them) or PNG (a server sends them). The browser's solver
//! and `vibe-rations serve` share it, so both draw the same pictures.

use std::path::Path;

use anyhow::{Context, Result, bail, ensure};

use crate::config::PlaneKind;
use crate::image;
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
    /// The value drawn at full colour.
    pub clip: f32,
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
            .nodes
            .iter()
            .map(|&n| g.position(n))
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

    /// The gather of `name` so far, one column per receiver and `rows` rows of time.
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

    /// The traces of `name` as an `.npy` file, (receivers, samples).
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
