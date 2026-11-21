//! 3D elastic wave propagation with staggered-grid finite differences.
//!
//! The solver integrates the velocity–stress form of the isotropic elastic wave equation
//! (Virieux 1986, Levander 1988): particle velocities and stresses live on a staggered grid,
//! of order 2–16 in space (8 by default) and second order (leapfrog) in time, whose time error
//! the traces can be cleared of (`tdt`). Convolutional PMLs absorb waves at the sides and
//! bottom, and the top is either absorbing or a free surface (stress imaging).

pub mod adjoint;
pub mod config;
pub mod cpml;
pub mod fd;
pub mod gradient;
pub mod grid;
pub mod image;
pub mod inversion;
pub mod memory;
pub mod misfit;
pub mod model;
pub mod npy;
pub mod receiver;
#[cfg(feature = "cli")]
pub mod remote;
#[cfg(feature = "cli")]
pub mod server;
pub mod simulation;
pub mod snapshot;
pub mod solver;
pub mod source;
pub mod tdt;
pub mod view;
