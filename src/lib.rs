//! 3D elastic wave propagation with a fourth-order staggered-grid finite-difference scheme.
//!
//! The solver integrates the velocity–stress form of the isotropic elastic wave equation
//! (Virieux 1986, Levander 1988): particle velocities and stresses live on a staggered grid,
//! fourth order in space and second order (leapfrog) in time. Convolutional PMLs absorb waves
//! at the sides and bottom, and the top is either absorbing or a free surface (stress imaging).

pub mod adjoint;
pub mod config;
pub mod cpml;
pub mod fd;
pub mod gradient;
pub mod grid;
pub mod image;
pub mod model;
pub mod npy;
pub mod receiver;
pub mod simulation;
pub mod snapshot;
pub mod solver;
pub mod source;
