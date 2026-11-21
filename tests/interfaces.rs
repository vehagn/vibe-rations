//! How well a planar interface between two solids reflects and transmits, depending on
//! where it falls between the nodes and on how a layered model is put onto them:
//! `just interfaces`.
//!
//! Plane P and S waves at normal incidence (`common::plane`), against the exact reflection
//! and transmission. The interface sits a fraction δ of a node spacing below a node. Not
//! part of the regular suite (a few minutes); `tests/physics.rs` checks the interface midway
//! between nodes. The research log's "A planar interface" entry records the results.

mod common;

use std::f64::consts::PI;

use common::plane::{LOWER, PlaneWave, UPPER, Wave};

const H: f64 = 10.0;
const DELTAS: [f64; 4] = [0.0, 0.25, 0.5, 0.75];

/// How a layered model is put onto the nodes.
#[derive(Clone, Copy, Debug)]
enum Sampling {
    /// As `[[model.layer]]`s are now: each node takes the medium it lies in, a node on the
    /// interface the one below.
    Node,
    /// Each node takes its cell's average over [z - h/2, z + h/2]: density arithmetically,
    /// the moduli ρvp² and ρvs² harmonically (as compliances).
    Cell,
    /// The step band-limited to the grid's Nyquist wavenumber, then sampled: a sinc
    /// tapered by a Hann window `w` nodes either side, applied to density and compliances.
    Band(f64),
    /// The same applied to the moduli themselves rather than their inverses.
    BandModuli(f64),
}

impl Sampling {
    fn name(self) -> String {
        match self {
            Sampling::Node => "node".into(),
            Sampling::Cell => "cell".into(),
            Sampling::Band(w) => format!("band ±{w}"),
            Sampling::BandModuli(w) => format!("band moduli ±{w}"),
        }
    }

    /// How much of a node u node spacings below the interface the lower medium takes.
    fn fraction_below(self, u: f64) -> f64 {
        match self {
            Sampling::Node => {
                if u >= 0.0 {
                    1.0
                } else {
                    0.0
                }
            }
            Sampling::Cell => (u + 0.5).clamp(0.0, 1.0),
            Sampling::Band(w) | Sampling::BandModuli(w) => {
                // The tapered sinc's integral from -w to u, normalised, by the midpoint rule.
                let kernel = |x: f64| {
                    let sinc = if x == 0.0 {
                        1.0
                    } else {
                        (PI * x).sin() / (PI * x)
                    };
                    sinc * (0.5 + 0.5 * (PI * x / w).cos())
                };
                let n = 4000;
                let dx = 2.0 * w / n as f64;
                let (mut total, mut part) = (0.0, 0.0);
                for m in 0..n {
                    let x = -w + (m as f64 + 0.5) * dx;
                    total += kernel(x);
                    if x < u {
                        part += kernel(x);
                    }
                }
                part / total
            }
        }
    }

    /// [vp, vs, ρ] at the node at depth z.
    fn medium(self, z: f64, interface: f64) -> [f64; 3] {
        let below = self.fraction_below((z - interface) / H);
        let above = 1.0 - below;
        let rho = above * UPPER[2] + below * LOWER[2];
        let modulus = |a: f64, b: f64| match self {
            Sampling::BandModuli(_) => above * a + below * b,
            _ => 1.0 / (above / a + below / b),
        };
        let m = modulus(UPPER[2] * UPPER[0].powi(2), LOWER[2] * LOWER[0].powi(2));
        let mu = modulus(UPPER[2] * UPPER[1].powi(2), LOWER[2] * LOWER[1].powi(2));
        [(m / rho).sqrt(), (mu / rho).sqrt(), rho]
    }
}

fn header() {
    println!(
        "{:<4} {:<16} {:>5} {:>5} | {:>8} | {:>8} {:>7} {:>7} | {:>8} {:>7} {:>7}",
        "wave", "sampling", "order", "δ", "direct", "R", "ratio", "shift", "T", "ratio", "shift"
    );
}

fn row(wave: Wave, sampling: Sampling, scale: usize, order: usize, delta: f64) {
    let mut p = PlaneWave::interface(wave, scale, delta);
    p.order = order;
    let interface = p.interface;
    let traces = match sampling {
        Sampling::Node => p.run(None),
        _ => p.run(Some(&|z| sampling.medium(z, interface))),
    };
    let [direct, r, t] = p.fits(&traces);
    let node = p.node_time();
    println!(
        "{:<4} {:<16} {:>5} {:>5.2} | {:>7.3}% | {:>7.3}% {:>7.4} {:>6.3}h | {:>7.3}% {:>7.4} {:>6.3}h",
        format!("{wave:?}"),
        sampling.name(),
        order,
        delta,
        100.0 * direct.misfit,
        100.0 * r.misfit,
        r.ratio,
        r.shift / node,
        100.0 * t.misfit,
        t.ratio,
        t.shift / node,
    );
}

/// Every way of sampling, at every δ, at twice the points per wavelength too. Misfits are
/// of each arrival, relative to it; ratios its amplitude relative to the exact; shifts in
/// the time a wave in the upper medium takes to cross a node ("h"), so a reflection's shift
/// is twice the interface's error in position.
#[test]
#[ignore = "a study: `just interfaces`"]
fn interfaces_by_position_and_sampling() {
    let samplings = [
        Sampling::Node,
        Sampling::Cell,
        Sampling::Band(4.0),
        Sampling::BandModuli(4.0),
    ];
    for scale in [1, 2] {
        let p = PlaneWave::interface(Wave::P, scale, 0.0);
        println!(
            "\nf0 = {} Hz, h = {H} m: P {:.1} and S {:.1} points per wavelength at f0",
            p.f0,
            UPPER[0] / p.f0 / H,
            UPPER[1] / p.f0 / H
        );
        header();
        for wave in [Wave::P, Wave::S] {
            for sampling in samplings {
                for delta in DELTAS {
                    row(wave, sampling, scale, 8, delta);
                }
            }
        }
    }
}

/// Whether the stencil's order matters: the node-sampled interface midway between nodes and
/// on one, and the band-limited one.
#[test]
#[ignore = "a study: `just interfaces`"]
fn interfaces_by_order() {
    header();
    for wave in [Wave::P, Wave::S] {
        for order in [4, 8, 12, 16] {
            for sampling in [Sampling::Node, Sampling::Band(4.0)] {
                for delta in [0.0, 0.5] {
                    row(wave, sampling, 1, order, delta);
                }
            }
        }
    }
}

/// How wide the band-limiting's taper needs to be.
#[test]
#[ignore = "a study: `just interfaces`"]
fn interfaces_by_band_limit_width() {
    header();
    for wave in [Wave::P, Wave::S] {
        for w in [2.0, 4.0, 8.0] {
            for delta in [0.0, 0.25] {
                row(wave, Sampling::Band(w), 1, 8, delta);
            }
        }
    }
}
