//! Snapshots: wavefield sections or volumes written during a run.

use std::path::Path;

use anyhow::Result;

use crate::config::PlaneKind;
use crate::image::{Section, write_section};
use crate::npy;
use crate::receiver::Component;
use crate::solver::Solver;

#[derive(Clone, Debug)]
pub struct Snapshot {
    pub component: Component,
    pub plane: PlaneKind,
    /// Node index along the plane's normal (unused for a volume).
    pub index: usize,
    /// Steps between snapshots.
    pub every: usize,
    pub png: bool,
    pub npy: bool,
    pub clip: Option<f32>,
}

impl Snapshot {
    /// Whether to take one after `step` steps.
    pub fn due(&self, step: usize) -> bool {
        step.is_multiple_of(self.every)
    }

    /// A base name such as `vz_xz60` (the section through y node 60) or `p_volume`.
    pub fn name(&self) -> String {
        let plane = match self.plane {
            PlaneKind::Xy => "xy",
            PlaneKind::Xz => "xz",
            PlaneKind::Yz => "yz",
            PlaneKind::Volume => return format!("{}_volume", self.component.name()),
        };
        format!("{}_{plane}{}", self.component.name(), self.index)
    }

    /// The nodes of the section, row by row from the top, as (width, height) and a function
    /// from (column, row) to node; for a volume, the whole grid with x fastest.
    fn nodes(&self, solver: &Solver) -> (Vec<usize>, Vec<[usize; 3]>) {
        let [nx, ny, nz] = solver.grid.n;
        let c = self.index;
        let (shape, nodes): (Vec<usize>, Vec<[usize; 3]>) = match self.plane {
            PlaneKind::Xy => (
                vec![ny, nx],
                grid2(nx, ny).map(|(i, j)| [i, j, c]).collect(),
            ),
            PlaneKind::Xz => (
                vec![nz, nx],
                grid2(nx, nz).map(|(i, k)| [i, c, k]).collect(),
            ),
            PlaneKind::Yz => (
                vec![nz, ny],
                grid2(ny, nz).map(|(j, k)| [c, j, k]).collect(),
            ),
            PlaneKind::Volume => (
                vec![nz, ny, nx],
                (0..nz)
                    .flat_map(|k| grid2(nx, ny).map(move |(i, j)| [i, j, k]))
                    .collect(),
            ),
        };
        (shape, nodes)
    }

    /// Writes the snapshot of the wavefield now into `dir`.
    pub fn write(&self, dir: &Path, solver: &Solver) -> Result<()> {
        let (shape, nodes) = self.nodes(solver);
        let values: Vec<f32> = nodes
            .iter()
            .map(|&n| solver.sample(self.component, n))
            .collect();
        let stem = format!("{}_{:06}", self.name(), solver.step);
        if self.npy {
            npy::write(&dir.join(format!("{stem}.npy")), &shape, &values)?;
        }
        if self.png && self.plane != PlaneKind::Volume {
            let shade = shade(solver, &nodes);
            let clip = self.clip.unwrap_or_else(|| robust_clip(&values));
            write_section(
                &dir.join(format!("{stem}.png")),
                &Section {
                    width: shape[1],
                    height: shape[0],
                    values: &values,
                    clip,
                    shade: &shade,
                },
            )?;
        }
        Ok(())
    }
}

/// A colour scale limit that a few extreme values can't dominate (such as the static stress a
/// source with a non-zero final moment leaves at its node): the 99.5th percentile of the
/// magnitudes, but at least 1e-4 of the largest, for when the wavefield covers almost none of
/// the section yet.
fn robust_clip(values: &[f32]) -> f32 {
    if values.is_empty() {
        return 1.0;
    }
    let mut mags: Vec<f32> = values.iter().map(|v| v.abs()).collect();
    let max = mags.iter().fold(0.0f32, |m, &v| m.max(v));
    let k = ((mags.len() as f64 * 0.995) as usize).min(mags.len() - 1);
    let (_, q, _) = mags.select_nth_unstable_by(k, f32::total_cmp);
    q.max(1e-4 * max)
}

/// (column, row) pairs for a w × h section, row by row.
fn grid2(w: usize, h: usize) -> impl Iterator<Item = (usize, usize)> {
    (0..h).flat_map(move |r| (0..w).map(move |c| (c, r)))
}

/// Shades of grey for the model under a section: faster rock darker, absorbing layers darker
/// still.
fn shade(solver: &Solver, nodes: &[[usize; 3]]) -> Vec<f32> {
    let g = solver.grid;
    let [lo, hi] = solver.model.vp_range();
    let span = (hi - lo).max(1e-9);
    nodes
        .iter()
        .map(|&[i, j, k]| {
            let vp = solver.model.vp[g.interior_idx(i, j, k)] as f64;
            let in_pml = [i, j, k]
                .iter()
                .zip(&solver.cpml)
                .any(|(&n, c)| c.slab(n).is_some());
            let s = 1.0 - 0.3 * (vp - lo) / span;
            (if in_pml { 0.85 * s } else { s }) as f32
        })
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn clips_past_a_few_spikes() {
        let mut v: Vec<f32> = (0..1000).map(|i| (i % 10) as f32 * 0.1).collect();
        // Two spikes a thousand times the wave, as a static source stress leaves.
        v[3] = 900.0;
        v[7] = -900.0;
        assert!((robust_clip(&v) - 0.9).abs() < 1e-6);
        // Almost all zeros: the floor keeps the scale from collapsing to nothing.
        let mut w = vec![0.0f32; 1000];
        w[0] = 2.0;
        assert_eq!(robust_clip(&w), 2e-4);
    }
}
