//! Pictures of wavefield sections over the model, and of receiver gathers: as RGBA pixels
//! (for the web app) and as PNG files.

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

/// Dark violet through blue, teal and green to yellow (viridis's stops), for t in [0, 1]: a
/// sequential scale whose lightness rises steadily, for models.
fn sequential(t: f32) -> [f32; 3] {
    const STOPS: [[f32; 3]; 5] = [
        [0.267, 0.005, 0.329],
        [0.229, 0.322, 0.546],
        [0.128, 0.567, 0.551],
        [0.369, 0.789, 0.383],
        [0.993, 0.906, 0.144],
    ];
    let t = if t.is_nan() { 0.0 } else { t.clamp(0.0, 1.0) } * 4.0;
    let i = (t.floor() as usize).min(3);
    let f = t - i as f32;
    std::array::from_fn(|c| STOPS[i][c] * (1.0 - f) + STOPS[i + 1][c] * f)
}

fn rgba(colours: impl Iterator<Item = [f32; 3]>) -> Vec<u8> {
    let mut out = Vec::new();
    for c in colours {
        out.extend(c.map(|c| (c * 255.0).round() as u8));
        out.push(255);
    }
    out
}

/// Values from `lo` (dark) to `hi` (light), RGBA; beyond them the ends' colours.
pub fn sequential_rgba(values: &[f32], lo: f32, hi: f32) -> Vec<u8> {
    let span = if hi > lo { hi - lo } else { 1.0 };
    rgba(values.iter().map(|&v| sequential((v - lo) / span)))
}

/// Signed values, RGBA: blue (negative) through white to red (positive), linear in value /
/// `clip`, so that two pictures with the same clip compare.
pub fn signed_rgba(values: &[f32], clip: f32) -> Vec<u8> {
    let clip = if clip > 0.0 { clip } else { 1.0 };
    rgba(values.iter().map(|&v| diverging(v / clip)))
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

/// A section's colours at one pixel per value, RGBA. The colour follows the signed square root
/// of value / clip, so weak arrivals show next to strong ones.
pub fn section_rgba(s: &Section) -> Vec<u8> {
    let clip = if s.clip > 0.0 { s.clip } else { 1.0 };
    let mut rgba = Vec::with_capacity(s.values.len() * 4);
    for (&v, &shade) in s.values.iter().zip(s.shade) {
        let v = v / clip;
        let c = diverging(v.signum() * v.abs().sqrt());
        rgba.extend(c.map(|c| (c * shade * 255.0).round() as u8));
        rgba.push(255);
    }
    rgba
}

/// Draws a section, scaled up by a whole number so its longer side is about 512 pixels.
pub fn write_section(path: &Path, s: &Section) -> Result<()> {
    let scale = (512 / s.width.max(s.height)).max(1);
    let (w, h) = (s.width * scale, s.height * scale);
    let pixels = section_rgba(s);
    let mut rgb = vec![0u8; w * h * 3];
    for y in 0..h {
        for x in 0..w {
            let p = (y / scale) * s.width + x / scale;
            rgb[(y * w + x) * 3..][..3].copy_from_slice(&pixels[p * 4..p * 4 + 3]);
        }
    }
    write_png(path, w, h, &rgb)
}

/// A gather's colours, RGBA: one column per receiver and `rows` rows spread over all
/// `samples`, time running down. Only the first `recorded` samples are drawn; the rest stay
/// pale grey. Each trace is scaled to its own peak so far (but not by more than 1000 times
/// the gather's peak), so late, weak arrivals show.
pub fn gather_rgba(
    receivers: usize,
    samples: usize,
    recorded: usize,
    data: &[f32],
    rows: usize,
) -> Vec<u8> {
    let recorded = recorded.min(samples);
    let peak = |t: &[f32]| t.iter().fold(0.0f32, |m, v| m.max(v.abs()));
    let traces: Vec<&[f32]> = data
        .chunks(samples.max(1))
        .map(|t| &t[..recorded])
        .collect();
    let global = traces
        .iter()
        .map(|t| peak(t))
        .fold(f32::MIN_POSITIVE, f32::max);
    let scales: Vec<f32> = traces
        .iter()
        .map(|t| 1.0 / peak(t).max(1e-3 * global))
        .collect();
    let mut rgba = Vec::with_capacity(receivers * rows * 4);
    for y in 0..rows {
        let m = y * samples / rows.max(1);
        for r in 0..receivers {
            if m < recorded {
                let c = diverging(traces[r][m] * scales[r]);
                rgba.extend(c.map(|c| (c * 255.0).round() as u8));
            } else {
                rgba.extend([236, 236, 236]);
            }
            rgba.push(255);
        }
    }
    rgba
}

/// Draws a gather: up to 8 pixels per receiver and up to 1000 rows.
pub fn write_gather(path: &Path, receivers: usize, samples: usize, data: &[f32]) -> Result<()> {
    let col = (600 / receivers).clamp(1, 8);
    let rows = samples.min(1000);
    let pixels = gather_rgba(receivers, samples, samples, data, rows);
    let (w, h) = (receivers * col, rows);
    let mut rgb = vec![0u8; w * h * 3];
    for y in 0..h {
        for x in 0..w {
            let p = y * receivers + x / col;
            rgb[(y * w + x) * 3..][..3].copy_from_slice(&pixels[p * 4..p * 4 + 3]);
        }
    }
    write_png(path, w, h, &rgb)
}

/// RGBA pixels as a PNG file in memory.
pub fn png_rgba(w: usize, h: usize, rgba: &[u8]) -> Result<Vec<u8>> {
    let mut out = Vec::new();
    let mut encoder = png::Encoder::new(&mut out, w as u32, h as u32);
    encoder.set_color(png::ColorType::Rgba);
    encoder.set_depth(png::BitDepth::Eight);
    let mut writer = encoder.write_header()?;
    writer.write_image_data(rgba)?;
    writer.finish()?;
    Ok(out)
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

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn colours_follow_the_sign_and_saturate() {
        assert_eq!(diverging(0.0), [1.0, 1.0, 1.0]);
        assert_eq!(diverging(1.0), [1.0, 0.0, 0.0]);
        assert_eq!(diverging(-1.0), [0.0, 0.0, 1.0]);
        assert_eq!(diverging(7.0), diverging(1.0));
        assert_eq!(diverging(f32::NAN), [1.0, 1.0, 1.0]);
    }

    #[test]
    fn model_and_signed_colours() {
        let rgba = sequential_rgba(&[1000.0, 2000.0, 3000.0, 5000.0], 1000.0, 3000.0);
        assert_eq!(&rgba[..4], &[68, 1, 84, 255]);
        assert_eq!(&rgba[8..12], &[253, 231, 37, 255]);
        assert_eq!(&rgba[12..16], &rgba[8..12], "saturates above hi");
        // Lightness rises along the scale.
        let light = |c: &[u8]| c[0] as u32 + c[1] as u32 + c[2] as u32;
        assert!(light(&rgba[..4]) < light(&rgba[4..8]) && light(&rgba[4..8]) < light(&rgba[8..12]));
        let rgba = signed_rgba(&[-2.0, -1.0, 0.0, 0.5], 2.0);
        assert_eq!(&rgba[..4], &[0, 0, 255, 255]);
        assert_eq!(
            &rgba[4..8],
            &[128, 128, 255, 255],
            "linear, not a square root"
        );
        assert_eq!(&rgba[8..12], &[255, 255, 255, 255]);
    }

    fn size(path: &Path) -> (u32, u32) {
        let file = std::io::BufReader::new(File::open(path).unwrap());
        let reader = png::Decoder::new(file).read_info().unwrap();
        (reader.info().width, reader.info().height)
    }

    #[test]
    fn scales_sections_up_and_fits_gathers() {
        let dir = std::env::temp_dir().join(format!("vibe-rations-image-{}", std::process::id()));
        std::fs::create_dir_all(&dir).unwrap();
        let values = vec![0.5f32; 40 * 20];
        let shade = vec![1.0f32; 40 * 20];
        let section = Section {
            width: 40,
            height: 20,
            values: &values,
            clip: 1.0,
            shade: &shade,
        };
        write_section(&dir.join("s.png"), &section).unwrap();
        // 512 / 40 = 12 times larger.
        assert_eq!(size(&dir.join("s.png")), (480, 240));
        // Three receivers get 8 pixels each; 2000 samples fold into 1000 rows.
        write_gather(&dir.join("g.png"), 3, 2000, &vec![1.0; 6000]).unwrap();
        assert_eq!(size(&dir.join("g.png")), (24, 1000));
        std::fs::remove_dir_all(dir).unwrap();
    }

    #[test]
    fn gathers_leave_unrecorded_samples_pale() {
        // Two receivers, four samples, two recorded; one row per sample.
        let data = [1.0, -1.0, 9.0, 9.0, 0.5, 0.0, 9.0, 9.0];
        let px = gather_rgba(2, 4, 2, &data, 4);
        assert_eq!(px.len(), 2 * 4 * 4);
        let at = |row: usize, col: usize| &px[(row * 2 + col) * 4..][..4];
        assert_eq!(at(0, 0), &[255, 0, 0, 255]); // +peak: red
        assert_eq!(at(1, 0), &[0, 0, 255, 255]); // -peak: blue
        assert_eq!(at(0, 1), &[255, 0, 0, 255]); // scaled to its own peak
        assert_eq!(at(2, 0), &[236, 236, 236, 255]); // not recorded yet
    }
}
