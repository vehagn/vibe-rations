//! PNG pictures: wavefield sections over the model, and receiver gathers.

use std::fs::File;
use std::io::BufWriter;
use std::path::Path;

use anyhow::{Context, Result};

/// Blue (negative) through white (zero) to red (positive), for v in [-1, 1].
fn diverging(v: f32) -> [f32; 3] {
    let t = if v.is_nan() { 0.0 } else { v.clamp(-1.0, 1.0) };
    if t >= 0.0 {
        [1.0, 1.0 - t, 1.0 - t]
    } else {
        [1.0 + t, 1.0 + t, 1.0]
    }
}

/// A 2D slice of a wavefield, row-major with the top row first.
pub struct Section<'a> {
    pub width: usize,
    pub height: usize,
    pub values: &'a [f32],
    /// The value drawn at full colour; larger magnitudes saturate.
    pub clip: f32,
    /// Brightness per value, 1 for white; darker shades show the model under the wavefield.
    pub shade: &'a [f32],
}

/// Draws a section, scaled up by a whole number so its longer side is about 512 pixels. The
/// colour follows the signed square root of value / clip, so weak arrivals show next to
/// strong ones.
pub fn write_section(path: &Path, s: &Section) -> Result<()> {
    let scale = (512 / s.width.max(s.height)).max(1);
    let (w, h) = (s.width * scale, s.height * scale);
    let clip = if s.clip > 0.0 { s.clip } else { 1.0 };
    let mut rgb = vec![0u8; w * h * 3];
    for y in 0..h {
        for x in 0..w {
            let p = (y / scale) * s.width + x / scale;
            let v = s.values[p] / clip;
            let c = diverging(v.signum() * v.abs().sqrt());
            for ch in 0..3 {
                rgb[(y * w + x) * 3 + ch] = (c[ch] * s.shade[p] * 255.0).round() as u8;
            }
        }
    }
    write_png(path, w, h, &rgb)
}

/// Draws a gather: one column per receiver, time running down. Each trace is scaled to its
/// own peak (but not by more than 1000 times the gather's peak), so late, weak arrivals show.
pub fn write_gather(path: &Path, receivers: usize, samples: usize, data: &[f32]) -> Result<()> {
    let peak = |t: &[f32]| t.iter().fold(0.0f32, |m, v| m.max(v.abs()));
    let global = peak(data).max(f32::MIN_POSITIVE);
    let scales: Vec<f32> = data
        .chunks(samples)
        .map(|t| 1.0 / peak(t).max(1e-3 * global))
        .collect();
    let col = (600 / receivers).clamp(1, 8);
    let rows = samples.min(1000);
    let (w, h) = (receivers * col, rows);
    let mut rgb = vec![0u8; w * h * 3];
    for y in 0..h {
        let m = y * samples / rows;
        for x in 0..w {
            let r = x / col;
            let c = diverging(data[r * samples + m] * scales[r]);
            for ch in 0..3 {
                rgb[(y * w + x) * 3 + ch] = (c[ch] * 255.0).round() as u8;
            }
        }
    }
    write_png(path, w, h, &rgb)
}

fn write_png(path: &Path, w: usize, h: usize, rgb: &[u8]) -> Result<()> {
    let file = File::create(path).with_context(|| format!("creating {}", path.display()))?;
    let mut encoder = png::Encoder::new(BufWriter::new(file), w as u32, h as u32);
    encoder.set_color(png::ColorType::Rgb);
    encoder.set_depth(png::BitDepth::Eight);
    let mut writer = encoder.write_header()?;
    writer.write_image_data(rgb)?;
    writer.finish()?;
    Ok(())
}
