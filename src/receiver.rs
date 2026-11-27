//! Receivers: traces of the wavefield at points (nodes, or anywhere between them).

use anyhow::{Context, Result};
use serde::{Deserialize, Serialize};

use crate::solver::{Field, PointOp, Solver};

/// A quantity to record.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Deserialize, Serialize)]
#[serde(rename_all = "lowercase")]
pub enum Component {
    /// Particle velocity, m/s.
    Vx,
    Vy,
    Vz,
    /// Pressure, -(σxx + σyy + σzz)/3, Pa. In a solid it carries only P waves.
    P,
}

impl Component {
    /// Every component.
    pub const ALL: [Component; 4] = [Component::Vx, Component::Vy, Component::Vz, Component::P];

    /// The component called `name` (vx, vy, vz or p), as configurations name it.
    pub fn from_name(name: &str) -> Result<Self> {
        Self::ALL
            .into_iter()
            .find(|c| c.name() == name)
            .with_context(|| format!("unknown field {name}; use vx, vy, vz or p"))
    }

    /// The wavefield arrays it is read from: its velocity, or for pressure the three normal
    /// stresses.
    pub fn fields(self) -> &'static [Field] {
        match self {
            Component::Vx => &[Field::Vx],
            Component::Vy => &[Field::Vy],
            Component::Vz => &[Field::Vz],
            Component::P => &[Field::Sxx, Field::Syy, Field::Szz],
        }
    }

    pub fn name(self) -> &'static str {
        match self {
            Component::Vx => "vx",
            Component::Vy => "vy",
            Component::Vz => "vz",
            Component::P => "p",
        }
    }
}

/// Every receiver records every component, `samples` times.
#[derive(Clone, Debug)]
pub struct Receivers {
    /// Where each receiver is, in grid coordinates ([`crate::grid::Grid::coordinates`]).
    pub points: Vec<[f64; 3]>,
    pub components: Vec<Component>,
    pub samples: usize,
    /// One array per component: receiver-major, (receivers, samples).
    pub traces: Vec<Vec<f32>>,
    /// How each component is read at each receiver: component-major, like `traces`.
    pub ops: Vec<Vec<PointOp>>,
    recorded: usize,
}

impl Receivers {
    pub fn new(
        solver: &Solver,
        points: Vec<[f64; 3]>,
        components: Vec<Component>,
        samples: usize,
    ) -> Self {
        let traces = components
            .iter()
            .map(|_| vec![0.0; points.len() * samples])
            .collect();
        let ops = components
            .iter()
            .map(|&c| points.iter().map(|&p| solver.receiver_op(c, p)).collect())
            .collect();
        Self {
            points,
            components,
            samples,
            traces,
            ops,
            recorded: 0,
        }
    }

    /// Records the next sample, if there is room for it.
    pub fn record(&mut self, solver: &Solver) {
        self.record_with(|op| op.sample(&solver.fields));
    }

    /// Records the next sample, if there is room for it, as `read` gives each operator's: it
    /// is asked for them in the order of [`Receivers::ops`], component by component.
    pub fn record_with(&mut self, mut read: impl FnMut(&PointOp) -> f32) {
        if self.recorded == self.samples {
            return;
        }
        let m = self.recorded;
        for (ops, data) in self.ops.iter().zip(&mut self.traces) {
            for (r, op) in ops.iter().enumerate() {
                data[r * self.samples + m] = read(op);
            }
        }
        self.recorded += 1;
    }

    /// Samples recorded so far.
    pub fn recorded(&self) -> usize {
        self.recorded
    }

    /// The trace of `component` at receiver `r`, or None if the component isn't recorded.
    pub fn trace(&self, component: Component, r: usize) -> Option<&[f32]> {
        let c = self.components.iter().position(|&c| c == component)?;
        Some(&self.traces[c][r * self.samples..(r + 1) * self.samples])
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::grid::Grid;
    use crate::solver::testing;

    fn solver() -> Solver {
        testing::solver(Grid::new([10, 10, 10], [10.0; 3]), 3, false)
    }

    fn close(a: &[f32], b: &[f32]) -> bool {
        a.len() == b.len() && a.iter().zip(b).all(|(x, y)| (x - y).abs() < 1e-5)
    }

    #[test]
    fn records_each_sample_into_its_slot_until_full() {
        let mut s = solver();
        let mut r = Receivers::new(
            &s,
            vec![[4.0, 4.0, 4.0], [5.0, 6.0, 7.0]],
            vec![Component::Vx, Component::P],
            3,
        );
        for n in 0..4 {
            // A uniform field reads the same at every node; the fourth sample doesn't fit.
            s.fields.vx.fill(n as f32);
            s.fields.sxx.fill(-3.0 * n as f32);
            r.record(&s);
        }
        assert_eq!(r.recorded(), 3);
        assert!(close(r.trace(Component::Vx, 1).unwrap(), &[0.0, 1.0, 2.0]));
        // p = -(σxx + σyy + σzz)/3, with only σxx set.
        assert!(close(r.trace(Component::P, 0).unwrap(), &[0.0, 1.0, 2.0]));
        assert_eq!(r.trace(Component::Vz, 0), None);
    }

    #[test]
    fn names_match_the_configuration() {
        for c in Component::ALL {
            let parsed: Component = serde_json::from_str(&format!("\"{}\"", c.name())).unwrap();
            assert_eq!(parsed, c);
            assert_eq!(Component::from_name(c.name()).unwrap(), c);
        }
        assert!(Component::from_name("vq").is_err());
    }
}
