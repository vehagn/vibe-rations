//! Model files from a 2D section (`vibe-rations extrude`). Published models such as Marmousi2
//! (Martin, Wiley & Marfurt 2006) come as SEG-Y files of vp, vs and ρ, a trace per x position
//! and a sample per depth, finely sampled. This averages them onto the grid's nodes and
//! repeats the section along y, for `[model.files]`: a 2.5D model, the same at every y,
//! through which the waves still spread in three dimensions.

use std::ops::Range;
use std::path::Path;

use anyhow::{Context, Result, bail, ensure};

/// One property of a 2D section: `traces` columns of `samples` values each, a trace at a time.
#[derive(Debug, Clone)]
pub struct Section {
    pub traces: usize,
    pub samples: usize,
    pub values: Vec<f32>,
}

impl Section {
    fn at(&self, trace: usize, sample: usize) -> f32 {
        self.values[trace * self.samples + sample]
    }
}

/// Reads a SEG-Y file: revision 0 or 1, big-endian, fixed-length traces of IBM or IEEE floats.
/// The trace headers are skipped: trace `t` is taken to stand `t` sample spacings along.
pub fn read_segy(path: &Path) -> Result<Section> {
    let bytes = std::fs::read(path).with_context(|| format!("reading {}", path.display()))?;
    parse_segy(&bytes).with_context(|| format!("in {}", path.display()))
}

/// The traces of a SEG-Y file held in memory; see [`read_segy`].
pub fn parse_segy(bytes: &[u8]) -> Result<Section> {
    ensure!(bytes.len() >= 3600, "too short for a SEG-Y file");
    // Offsets into the binary header, which follows the 3200-byte text header.
    let be16 = |at: usize| i16::from_be_bytes([bytes[at], bytes[at + 1]]);
    // Unsigned: some writers go past 32767 samples.
    let samples = be16(3220) as u16;
    let format = be16(3224);
    let extended = be16(3504);
    ensure!(
        samples > 0,
        "the binary header gives {samples} samples per trace"
    );
    let float: fn(u32) -> f32 = match format {
        1 => ibm_to_f32,
        5 => f32::from_bits,
        _ => bail!(
            "sample format {format} isn't supported, only 1 (IBM float) and 5 (IEEE float); \
             or the file is little-endian"
        ),
    };
    ensure!(
        extended >= 0,
        "a variable number of extended text headers isn't supported"
    );
    let samples = samples as usize;
    let start = 3600 + 3200 * extended as usize;
    let trace_len = 240 + 4 * samples;
    let body = bytes.len().saturating_sub(start);
    ensure!(
        body > 0 && body.is_multiple_of(trace_len),
        "{body} bytes of traces don't make whole traces of {samples} samples"
    );
    let traces = body / trace_len;
    let mut values = Vec::with_capacity(traces * samples);
    for t in 0..traces {
        let data = &bytes[start + t * trace_len + 240..start + (t + 1) * trace_len];
        values.extend(
            data.as_chunks::<4>()
                .0
                .iter()
                .map(|&w| float(u32::from_be_bytes(w))),
        );
    }
    Ok(Section {
        traces,
        samples,
        values,
    })
}

/// An IBM System/360 single-precision float: a sign, an exponent of 16 biased by 64 and a
/// 24-bit fraction. SEG-Y's sample format 1.
pub fn ibm_to_f32(w: u32) -> f32 {
    let fraction = (w & 0x00ff_ffff) as f64 / (1u32 << 24) as f64;
    let exponent = ((w >> 24) & 0x7f) as i32 - 64;
    let v = fraction * 16f64.powi(exponent);
    (if w >> 31 == 1 { -v } else { v }) as f32
}

/// How a section becomes model files.
#[derive(Debug, Clone)]
pub struct Extrusion {
    /// The section's sample spacing in metres, along x and z alike.
    pub input_spacing: f64,
    /// The nodes' spacing in metres.
    pub spacing: f64,
    /// Nodes along y to repeat the section over; none to keep it a section, (nz, nx), which
    /// `[model.files]` repeats along y itself.
    pub ny: Option<usize>,
    /// The part of the section to keep along x, in metres from its first trace; all of it by
    /// default.
    pub x: Option<[f64; 2]>,
    /// The part to keep along z, in metres from its first sample, which becomes z = 0.
    pub z: Option<[f64; 2]>,
    /// Factors that bring the speeds to m/s and the densities to kg/m³, e.g. 1000 for km/s and
    /// g/cm³.
    pub speed_scale: f64,
    pub rho_scale: f64,
    /// Smooth the solid part afterwards with a Gaussian of this standard deviation, in metres:
    /// a starting model for an inversion.
    pub smooth: Option<f64>,
}

/// A model on the nodes, x fastest, then y, then z: as `[model.files]` reads it.
#[derive(Debug, Clone)]
pub struct Extruded {
    /// Nodes along x, y and z; along y, 1 for a section.
    pub n: [usize; 3],
    /// The arrays' shape: (nz, ny, nx), or (nz, nx) for a section.
    pub shape: Vec<usize>,
    pub vp: Vec<f32>,
    pub vs: Vec<f32>,
    pub rho: Vec<f32>,
    /// How many of the section's nodes (one y plane's) are fluid.
    pub fluid_nodes: usize,
}

/// Averages the section's properties over each node's cell, a square `spacing` wide centred on
/// it (where the section has no sample in a cell, the nearest sample stands for it), and
/// repeats the result along y.
///
/// Density is the arithmetic mean, the P modulus ρvp² the harmonic mean: Backus's (1962)
/// long-wave averages for waves crossing fine layers, which is all a grid this coarse can tell
/// apart. A node is fluid if most of its samples are (vs = 0); otherwise its S modulus ρvs² is
/// the harmonic mean over the solid samples. Backus's average would make every node a fluid
/// sample reaches fluid, moving a seafloor down by up to half a node; by majority it stays
/// within half a node either way.
///
/// With `smooth`, each property is then smoothed over the solid nodes only (see
/// [`smooth_solid`]), so the fluid and the seafloor stay sharp, as an inversion that holds the
/// water fixed wants them.
pub fn extrude(vp: &Section, vs: &Section, rho: &Section, e: &Extrusion) -> Result<Extruded> {
    for s in [vs, rho] {
        ensure!(
            (s.traces, s.samples) == (vp.traces, vp.samples),
            "the sections differ in size: vp has {} traces of {} samples, another {} of {}",
            vp.traces,
            vp.samples,
            s.traces,
            s.samples
        );
    }
    ensure!(
        e.input_spacing > 0.0 && e.spacing > 0.0,
        "spacings must be positive"
    );
    ensure!(e.ny != Some(0), "the model needs at least one node along y");
    let d = e.input_spacing;
    let axis = |range: Option<[f64; 2]>, len: usize, name: &str| -> Result<Vec<Range<usize>>> {
        let end = (len - 1) as f64 * d;
        let [from, to] = range.unwrap_or([0.0, end]);
        let slack = 1e-6 * d;
        ensure!(
            0.0 <= from + slack && from < to && to <= end + slack,
            "{name} from {from} to {to} m isn't inside the section, 0 to {end} m"
        );
        let n = ((to - from) / e.spacing + 1e-9).floor() as usize + 1;
        Ok((0..n)
            .map(|i| cell(from + i as f64 * e.spacing, e.spacing, d, len))
            .collect())
    };
    let xs = axis(e.x, vp.traces, "x")?;
    let zs = axis(e.z, vp.samples, "z")?;
    let (nx, nz) = (xs.len(), zs.len());
    let mut plane = [
        vec![0f32; nx * nz],
        vec![0f32; nx * nz],
        vec![0f32; nx * nz],
    ];
    let mut fluid_nodes = 0;
    for (k, zr) in zs.iter().enumerate() {
        for (i, xr) in xs.iter().enumerate() {
            let (mut n, mut fluid) = (0usize, 0usize);
            let (mut rho_sum, mut p_compliance, mut s_compliance) = (0.0f64, 0.0f64, 0.0f64);
            for t in xr.clone() {
                for s in zr.clone() {
                    let r = rho.at(t, s) as f64 * e.rho_scale;
                    let p = vp.at(t, s) as f64 * e.speed_scale;
                    let sv = vs.at(t, s) as f64 * e.speed_scale;
                    ensure!(
                        r > 0.0 && p > 0.0 && sv >= 0.0,
                        "unphysical sample at trace {t}, sample {s}: vp {p}, vs {sv}, rho {r}"
                    );
                    n += 1;
                    rho_sum += r;
                    p_compliance += 1.0 / (r * p * p);
                    if sv == 0.0 {
                        fluid += 1;
                    } else {
                        s_compliance += 1.0 / (r * sv * sv);
                    }
                }
            }
            let r = rho_sum / n as f64;
            let p = (n as f64 / p_compliance / r).sqrt();
            let sv = if 2 * fluid > n {
                fluid_nodes += 1;
                0.0
            } else {
                ((n - fluid) as f64 / s_compliance / r).sqrt()
            };
            ensure!(
                3.0 * p * p > 4.0 * sv * sv,
                "the average at node ({i}, {k}) is unphysical: vp {p}, vs {sv}"
            );
            let at = k * nx + i;
            plane[0][at] = p as f32;
            plane[1][at] = sv as f32;
            plane[2][at] = r as f32;
        }
    }
    if let Some(sigma) = e.smooth {
        ensure!(sigma > 0.0, "the smoothing length must be positive");
        let solid: Vec<bool> = plane[1].iter().map(|&v| v > 0.0).collect();
        for p in &mut plane {
            smooth_solid(p, &solid, nx, sigma / e.spacing);
        }
    }
    let [vp_max, rho_max] = [&plane[0], &plane[2]].map(|v| v.iter().copied().fold(0.0, f32::max));
    ensure!(
        vp_max >= 100.0,
        "the fastest P speed is {vp_max} m/s: are the speeds in km/s? Scale them by 1000"
    );
    ensure!(
        rho_max >= 100.0,
        "the largest density is {rho_max} kg/m³: is it in g/cm³? Scale it by 1000"
    );
    // The same plane at every y.
    let ny = e.ny.unwrap_or(1);
    let [vp, vs, rho] = plane.map(|p| {
        let mut out = Vec::with_capacity(nx * ny * nz);
        for row in p.chunks_exact(nx) {
            for _ in 0..ny {
                out.extend_from_slice(row);
            }
        }
        out
    });
    Ok(Extruded {
        n: [nx, ny, nz],
        shape: match e.ny {
            Some(ny) => vec![nz, ny, nx],
            None => vec![nz, nx],
        },
        vp,
        vs,
        rho,
        fluid_nodes,
    })
}

/// Smooths the solid nodes of a plane (rows of `nx`, x fastest) with a Gaussian of standard
/// deviation `sigma` nodes, weighting solid nodes only (a normalised convolution): the fluid
/// nodes keep their values, and none of them leaks into the solid. A weighted mean of values
/// that each have vp² > 4/3 vs² has it too, so smoothing keeps the model physical.
pub fn smooth_solid(plane: &mut [f32], solid: &[bool], nx: usize, sigma: f64) {
    let nz = plane.len() / nx;
    let r = (3.0 * sigma).ceil() as usize;
    let w: Vec<f64> = (0..=r)
        .map(|d| (-0.5 * (d as f64 / sigma).powi(2)).exp())
        .collect();
    // One pass of the kernel along x (`step` 1) or z (`step` nx), truncated at the edges.
    let blur = |v: &[f64], along_x: bool| -> Vec<f64> {
        let (n, step) = if along_x { (nx, 1) } else { (nz, nx) };
        (0..v.len())
            .map(|at| {
                let c = if along_x { at % nx } else { at / nx };
                let base = at - c * step;
                (c.saturating_sub(r)..=(c + r).min(n - 1))
                    .map(|m| w[c.abs_diff(m)] * v[base + m * step])
                    .sum()
            })
            .collect()
    };
    let mask: Vec<f64> = solid.iter().map(|&s| if s { 1.0 } else { 0.0 }).collect();
    let masked: Vec<f64> = plane
        .iter()
        .zip(&mask)
        .map(|(&v, m)| v as f64 * m)
        .collect();
    let sum = blur(&blur(&masked, true), false);
    let weight = blur(&blur(&mask, true), false);
    for (at, v) in plane.iter_mut().enumerate() {
        if solid[at] {
            *v = (sum[at] / weight[at]) as f32;
        }
    }
}

/// The samples (spaced `d`, `len` of them) in the cell `[c - h/2, c + h/2)`, or the one nearest
/// `c` if there are none.
fn cell(c: f64, h: f64, d: f64, len: usize) -> Range<usize> {
    let index = |x: f64| ((x / d - 1e-9).ceil().max(0.0) as usize).min(len);
    let (lo, hi) = (index(c - h / 2.0), index(c + h / 2.0));
    if lo < hi {
        lo..hi
    } else {
        let s = ((c / d).round().max(0.0) as usize).min(len - 1);
        s..s + 1
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// An IBM float for a value that has one (tests only: truncates the fraction).
    fn f32_to_ibm(v: f32) -> u32 {
        if v == 0.0 {
            return 0;
        }
        let sign = if v < 0.0 { 1u32 << 31 } else { 0 };
        let mut m = v.abs() as f64;
        let mut e = 64i32;
        while m >= 1.0 {
            m /= 16.0;
            e += 1;
        }
        while m < 1.0 / 16.0 {
            m *= 16.0;
            e -= 1;
        }
        sign | (e as u32) << 24 | (m * (1u32 << 24) as f64) as u32
    }

    /// A SEG-Y file of `traces` (each a column of samples) in `format` 1 or 5.
    fn segy(traces: &[Vec<f32>], format: i16) -> Vec<u8> {
        let samples = traces[0].len();
        let mut b = vec![0u8; 3600];
        b[3220..3222].copy_from_slice(&(samples as i16).to_be_bytes());
        b[3224..3226].copy_from_slice(&format.to_be_bytes());
        for t in traces {
            b.extend_from_slice(&[0u8; 240]);
            for &v in t {
                let w = if format == 1 {
                    f32_to_ibm(v)
                } else {
                    v.to_bits()
                };
                b.extend_from_slice(&w.to_be_bytes());
            }
        }
        b
    }

    #[test]
    fn decodes_ibm_floats() {
        // The textbook example, and values near Marmousi2's.
        assert_eq!(ibm_to_f32(0xc276_a000), -118.625);
        assert_eq!(ibm_to_f32(0x4110_0000), 1.0);
        assert_eq!(ibm_to_f32(0x4264_0000), 100.0);
        assert_eq!(ibm_to_f32(0), 0.0);
        for v in [1.5f32, 2.802, 0.291_08, 1.01, 4700.0] {
            let back = ibm_to_f32(f32_to_ibm(v));
            assert!((back - v).abs() <= 1e-6 * v, "{v} came back as {back}");
        }
    }

    #[test]
    fn reads_ibm_and_ieee_sections_alike() {
        let traces: Vec<Vec<f32>> = (0..3)
            .map(|t| (0..5).map(|s| 1.0 + t as f32 + 0.25 * s as f32).collect())
            .collect();
        for format in [1, 5] {
            let s = parse_segy(&segy(&traces, format)).unwrap();
            assert_eq!((s.traces, s.samples), (3, 5));
            assert_eq!(s.at(2, 3), 3.75);
            assert_eq!(s.values, traces.concat());
        }
        let e = parse_segy(&segy(&traces, 3)).unwrap_err().to_string();
        assert!(e.contains("format 3"), "{e}");
        let mut cut = segy(&traces, 5);
        cut.pop();
        assert!(parse_segy(&cut).is_err());
    }

    fn uniform(traces: usize, samples: usize, f: impl Fn(usize, usize) -> f32) -> Section {
        Section {
            traces,
            samples,
            values: (0..traces)
                .flat_map(|t| (0..samples).map(move |s| (t, s)))
                .map(|(t, s)| f(t, s))
                .collect(),
        }
    }

    fn extrusion(spacing: f64, ny: Option<usize>) -> Extrusion {
        Extrusion {
            input_spacing: 1.0,
            spacing,
            ny,
            x: None,
            z: None,
            speed_scale: 1000.0,
            rho_scale: 1000.0,
            smooth: None,
        }
    }

    #[test]
    fn smooths_the_solid_without_touching_the_fluid() {
        // Two rows of water over a solid whose left half is slower than its right.
        let (nx, nz) = (20, 12);
        let solid: Vec<bool> = (0..nx * nz).map(|at| at / nx >= 2).collect();
        let original: Vec<f32> = (0..nx * nz)
            .map(|at| match (solid[at], at % nx < 10) {
                (false, _) => 1500.0,
                (true, true) => 2000.0,
                (true, false) => 3000.0,
            })
            .collect();
        let mut v = original.clone();
        smooth_solid(&mut v, &solid, nx, 2.0);
        for at in 0..nx * nz {
            let (i, k) = (at % nx, at / nx);
            if !solid[at] {
                assert_eq!(v[at], 1500.0);
            } else if !(3..17).contains(&i) {
                // Far from the step, right under the water too: unchanged, the water left out.
                assert!((v[at] - original[at]).abs() < 1e-3, "({i}, {k}): {}", v[at]);
            } else {
                assert!((2000.0..=3000.0).contains(&v[at]), "({i}, {k}): {}", v[at]);
            }
        }
        // Across the step the values climb steadily.
        let row: Vec<f32> = (0..nx).map(|i| v[6 * nx + i]).collect();
        assert!(row.windows(2).all(|p| p[1] > p[0] - 1e-3), "{row:?}");
        assert!(row[9] > 2100.0 && row[10] < 2900.0, "{row:?}");
    }

    #[test]
    fn averages_cells_and_repeats_them_along_y() {
        // 9 traces × 13 samples 1 m apart, water above z = 4.5 m, two solids below, in km/s
        // and g/cm³.
        let water = |s: usize| s <= 4;
        let vp = uniform(9, 13, |_, s| if water(s) { 1.5 } else { 2.0 });
        let vs = uniform(9, 13, |t, s| match (water(s), t < 4) {
            (true, _) => 0.0,
            (false, true) => 0.5,
            (false, false) => 1.0,
        });
        let rho = uniform(9, 13, |_, s| if water(s) { 1.0 } else { 2.0 });
        let m = extrude(&vp, &vs, &rho, &extrusion(4.0, Some(3))).unwrap();
        // Nodes at x = 0, 4, 8 and z = 0, 4, 8, 12.
        assert_eq!((m.n, m.shape.as_slice()), ([3, 3, 4], &[4, 3, 3][..]));
        assert_eq!(m.vp.len(), 36);
        let at = |v: &[f32], i: usize, j: usize, k: usize| v[(k * 3 + j) * 3 + i];
        for j in 0..3 {
            // z = 4: samples 2 to 5, three of them water: a fluid node, with the harmonic P
            // modulus and the mean density.
            assert_eq!(at(&m.vs, 1, j, 1), 0.0);
            let r = (3.0 * 1000.0 + 2000.0) / 4.0;
            let p = 4.0 / (3.0 / (1000.0 * 1500f64.powi(2)) + 1.0 / (2000.0 * 2000f64.powi(2)));
            assert!((at(&m.rho, 1, j, 1) as f64 - r).abs() < 1e-3);
            assert!((at(&m.vp, 1, j, 1) as f64 - (p / r).sqrt()).abs() < 1e-3);
            // z = 8, x = 4: traces 2 to 5, two at vs = 0.5 km/s and two at 1: the harmonic S
            // modulus (density is the same throughout).
            let s = (2.0 / (1.0 / 500f64.powi(2) + 1.0 / 1000f64.powi(2))).sqrt();
            assert!((at(&m.vs, 1, j, 2) as f64 - s).abs() < 1e-3);
            // The edge cell at x = 0 holds samples 0 and 1 only.
            assert_eq!(at(&m.vs, 0, j, 2), 500.0);
        }
        assert_eq!(m.fluid_nodes, 3 + 3);
    }

    #[test]
    fn crops_and_takes_the_nearest_sample_when_refining() {
        let vp = uniform(10, 10, |t, s| 2.0 + 0.1 * t as f32 + 0.01 * s as f32);
        let vs = uniform(10, 10, |_, _| 1.0);
        let rho = uniform(10, 10, |_, _| 2.0);
        let e = Extrusion {
            x: Some([2.0, 4.0]),
            z: Some([3.0, 4.0]),
            ..extrusion(0.5, None)
        };
        let m = extrude(&vp, &vs, &rho, &e).unwrap();
        // Without ny, a section: (nz, nx).
        assert_eq!((m.n, m.shape.as_slice()), ([5, 1, 3], &[3, 5][..]));
        // The node at x = 2.5 lies between samples and its cell [2.25, 2.75) holds none: the
        // nearest one, rounding up, stands for it.
        assert!((m.vp[1] - 2330.0).abs() < 1e-3, "{}", m.vp[1]);
        assert!((m.vp[5 + 4] - 2440.0).abs() < 1e-3, "{}", m.vp[9]);
    }

    #[test]
    fn refuses_what_it_cannot_make_sense_of() {
        let one = |v: f32| uniform(4, 4, move |_, _| v);
        let (vp, vs, rho) = (one(2.0), one(1.0), one(2.0));
        let err = |e: Extrusion| extrude(&vp, &vs, &rho, &e).unwrap_err().to_string();
        assert!(
            err(Extrusion {
                speed_scale: 1.0,
                ..extrusion(1.0, Some(1))
            })
            .contains("km/s")
        );
        assert!(
            err(Extrusion {
                rho_scale: 1.0,
                ..extrusion(1.0, Some(1))
            })
            .contains("g/cm³")
        );
        assert!(
            err(Extrusion {
                x: Some([1.0, 5.0]),
                ..extrusion(1.0, Some(1))
            })
            .contains("inside the section")
        );
        let short = uniform(4, 3, |_, _| 1.0);
        assert!(
            extrude(&vp, &short, &rho, &extrusion(1.0, Some(1)))
                .unwrap_err()
                .to_string()
                .contains("differ")
        );
        // vs too close to vp: no positive bulk modulus.
        let fast = one(1.9);
        assert!(
            extrude(&vp, &fast, &rho, &extrusion(1.0, Some(1)))
                .unwrap_err()
                .to_string()
                .contains("unphysical")
        );
    }
}
