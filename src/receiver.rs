//! Receivers: traces of the wavefield at nodes.

use serde::{Deserialize, Serialize};

use crate::solver::{PointOp, Solver};

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
    pub nodes: Vec<[usize; 3]>,
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
        nodes: Vec<[usize; 3]>,
        components: Vec<Component>,
        samples: usize,
    ) -> Self {
        let traces = components
            .iter()
            .map(|_| vec![0.0; nodes.len() * samples])
            .collect();
        let ops = components
            .iter()
            .map(|&c| nodes.iter().map(|&n| solver.receiver_op(c, n)).collect())
            .collect();
        Self {
            nodes,
            components,
            samples,
            traces,
            ops,
            recorded: 0,
        }
    }

    /// Records the next sample, if there is room for it.
    pub fn record(&mut self, solver: &Solver) {
        if self.recorded == self.samples {
            return;
        }
        let m = self.recorded;
        for (ops, data) in self.ops.iter().zip(&mut self.traces) {
            for (r, op) in ops.iter().enumerate() {
                data[r * self.samples + m] = op.sample(&solver.fields);
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
    use crate::cpml::PmlSettings;
    use crate::grid::Grid;
    use crate::model::Model;

    fn solver() -> Solver {
        let g = Grid::new([10, 10, 10], [10.0; 3]);
        let pml = PmlSettings {
            power: 2.0,
            reflection: None,
            kappa_max: 1.0,
            alpha_max: 30.0,
            vp_max: 3000.0,
            sponge: 0.0,
        };
        Solver::new(
            Model::homogeneous(g, 3000.0, 1700.0, 2500.0),
            1e-3,
            3,
            false,
            &pml,
        )
    }

    fn close(a: &[f32], b: &[f32]) -> bool {
        a.len() == b.len() && a.iter().zip(b).all(|(x, y)| (x - y).abs() < 1e-5)
    }

    #[test]
    fn records_each_sample_into_its_slot_until_full() {
        let mut s = solver();
        let mut r = Receivers::new(
            &s,
            vec![[4, 4, 4], [5, 6, 7]],
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
        for c in [Component::Vx, Component::Vy, Component::Vz, Component::P] {
            let parsed: Component = serde_json::from_str(&format!("\"{}\"", c.name())).unwrap();
            assert_eq!(parsed, c);
        }
    }
}
