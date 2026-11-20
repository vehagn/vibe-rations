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
