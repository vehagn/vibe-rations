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
        match self.plane {
            PlaneKind::Volume => format!("{}_volume", self.component.name()),
            plane => format!("{}_{}{}", self.component.name(), plane.name(), self.index),
        }
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

    /// The snapshot of the wavefield now.
    pub fn capture(&self, solver: &Solver) -> Capture {
        let (shape, nodes) = self.nodes(solver);
        let values = nodes
            .iter()
            .map(|&n| solver.sample(self.component, n))
            .collect();
        Capture {
            shape,
            values,
            shade: shade(solver, &nodes),
        }
    }

    /// The colour scale limit for a capture: the configured one, or a robust default.
    pub fn clip_for(&self, capture: &Capture) -> f32 {
        self.clip.unwrap_or_else(|| robust_clip(&capture.values))
    }

    /// Writes the snapshot of the wavefield now into `dir`.
    pub fn write(&self, dir: &Path, solver: &Solver) -> Result<()> {
        let c = self.capture(solver);
        let stem = format!("{}_{:06}", self.name(), solver.step);
        if self.npy {
            npy::write(&dir.join(format!("{stem}.npy")), &c.shape, &c.values)?;
        }
        if self.png && self.plane != PlaneKind::Volume {
            write_section(
                &dir.join(format!("{stem}.png")),
                &c.section(self.clip_for(&c)),
            )?;
        }
        Ok(())
    }
}

/// A snapshot's values, with the model's shading under them.
#[derive(Clone, Debug)]
pub struct Capture {
    /// (rows, columns) for a section, the top row first; (nz, ny, nx) for a volume.
    pub shape: Vec<usize>,
    pub values: Vec<f32>,
    /// Brightness of the model under each value: faster rock and the absorbing layers darker.
    pub shade: Vec<f32>,
}

impl Capture {
    /// The capture as a section to draw, for a plane (not a volume).
    pub fn section(&self, clip: f32) -> Section<'_> {
        Section {
            width: self.shape[1],
            height: self.shape[0],
            values: &self.values,
            clip,
            shade: &self.shade,
        }
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

    use crate::grid::Grid;
    use crate::solver::testing;

    fn snapshot(plane: PlaneKind, index: usize) -> Snapshot {
        Snapshot {
            component: Component::P,
            plane,
            index,
            every: 3,
            png: false,
            npy: true,
            clip: None,
        }
    }

    #[test]
    fn names_say_what_and_where() {
        assert_eq!(snapshot(PlaneKind::Xz, 60).name(), "p_xz60");
        assert_eq!(snapshot(PlaneKind::Xy, 0).name(), "p_xy0");
        assert_eq!(snapshot(PlaneKind::Volume, 0).name(), "p_volume");
        assert!(snapshot(PlaneKind::Yz, 1).due(6));
        assert!(!snapshot(PlaneKind::Yz, 1).due(7));
    }

    /// Sections are row-major with the first row at the top (smallest y or z), so depth runs
    /// down the rows of an xz section.
    #[test]
    fn sections_are_row_major_from_the_top() {
        let g = Grid::new([8, 7, 6], [10.0; 3]);
        let mut s = testing::solver(g, 0, false);
        // p = k at every node.
        for k in 0..6 {
            for j in 0..7 {
                for i in 0..8 {
                    s.fields.sxx[g.idx(i, j, k)] = -3.0 * k as f32;
                }
            }
        }
        let dir = std::env::temp_dir().join(format!("vibe-rations-snap-{}", std::process::id()));
        std::fs::create_dir_all(&dir).unwrap();
        for (plane, index, shape) in [
            (PlaneKind::Xy, 2, vec![7, 8]),
            (PlaneKind::Xz, 3, vec![6, 8]),
            (PlaneKind::Yz, 4, vec![6, 7]),
            (PlaneKind::Volume, 0, vec![6, 7, 8]),
        ] {
            let snap = snapshot(plane, index);
            snap.write(&dir, &s).unwrap();
            let (got, values) =
                npy::read(&dir.join(format!("{}_000000.npy", snap.name()))).unwrap();
            assert_eq!(got, shape, "{plane:?}");
            if plane == PlaneKind::Xz {
                assert!((values[8 * 5] - 5.0).abs() < 1e-5, "row 5 should be k = 5");
            }
        }
        std::fs::remove_dir_all(dir).unwrap();
    }
}
