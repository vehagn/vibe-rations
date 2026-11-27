//! A closed surface one node wide on which the forward wavefield is recorded, so that it can
//! be rebuilt backwards in time inside the surface by injecting the records as sources: the
//! elastodynamic Kirchhoff–Helmholtz reconstruction of Mittet (1994) in the velocity–stress
//! form of Aaker, Raknes, Pedersen & Arntsen (2020), whose sources are the same as Vasmel &
//! Robertsson's (2016) multiple point sources.
//!
//! The idea, discretely: the field truncated to the box Ω, w = Xu with X the box's indicator
//! (½ on its faces), satisfies the solver's own equations plus the commutators of X with the
//! derivatives, B [X, D] σ and C [X, D] v, which live only around the faces. In the continuum
//! they are a force n·σ and a stress source c : (v ⊗ n) on the surface (n its outward normal),
//! which needs only v and n·σ there: six values per face node. On the grid the commutator of
//! an order-2L stencil spreads over 2L nodes on either side, so the one-node surface stands
//! in for it: each recorded value is read onto the face's node plane with the stencil's own
//! interpolation and spread back with its transpose, along the face's normal only (each
//! field keeps its own lattice along the face). At second order this is exact; at higher
//! orders it is right for well-sampled waves and leaks a little, which [`Surface`]'s users
//! measure. The faces sit on node planes, where the indicator of a node-lattice field is ½,
//! so that every lattice's step lies on the face. Along a free surface the box is open: the
//! surface closes it, and the mirrors carry the truncated field as they carry the whole one.
//!
//! The rebuild runs the time-reversed problem forward with the solver's ordinary step, layers
//! and all: with ṽ = v and σ̃ = -σ a forward step of (ṽ, σ̃) undoes a step of (v, σ), if the
//! sources are withdrawn (forces) or added (moments) instead. The rebuilt field leaks outside
//! Ω, where the absorbing layers take it, unlike the band's exact step back
//! ([`crate::solver::Band`]), which can't run its layers backwards and so must keep a band as
//! deep as the stencil reaches.

use anyhow::{Result, ensure};

use crate::grid::Grid;
use crate::solver::{Field, Fields, PointOp, Solver};

/// The box Ω: node planes `lo` and `hi` along each axis, inclusive, or no `lo` along z with
/// a free surface, which closes it there.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct Bounds {
    pub lo: [Option<usize>; 3],
    pub hi: [usize; 3],
}

impl Bounds {
    /// The box `inset` nodes inside `interior` (the nodes outside the absorbing layers) on
    /// every side but a free surface. It must be at least the stencil's half-length `reach`
    /// in from the layers, so that everything the injection reaches is updated without
    /// them.
    pub fn new(
        interior: &[std::ops::Range<usize>; 3],
        free_surface: bool,
        reach: usize,
        inset: usize,
    ) -> Result<Self> {
        ensure!(
            inset >= reach,
            "the surface must be at least {reach} nodes in from the layers, not {inset}"
        );
        let mut lo = [None; 3];
        let mut hi = [0; 3];
        for a in 0..3 {
            let r = &interior[a];
            let low = r.start + inset;
            ensure!(
                r.end > low + inset + 2,
                "the surface {inset} nodes in from the layers leaves no room inside along {}",
                ["x", "y", "z"][a]
            );
            hi[a] = r.end - 1 - inset;
            if !(a == 2 && free_surface) {
                lo[a] = Some(low);
            }
        }
        Ok(Self { lo, hi })
    }

    /// The share of a value at grid coordinate `x` along `axis` inside the box: 1 inside, ½ on
    /// a face, 0 outside. Above a free surface (x < 0) it is inside, as the mirror images are.
    pub fn share(&self, axis: usize, x: f64) -> f32 {
        let hi = self.hi[axis] as f64;
        let inside_lo = match self.lo[axis] {
            Some(lo) => {
                let lo = lo as f64;
                if x < lo {
                    return 0.0;
                }
                x > lo
            }
            None => true,
        };
        if x > hi {
            0.0
        } else if x == hi || !inside_lo {
            0.5
        } else {
            1.0
        }
    }

    /// The share of `field` at padded index `p` inside the box: the product of
    /// [`Bounds::share`] along the three axes.
    pub fn field_share(&self, g: &Grid, field: Field, p: usize) -> f32 {
        let (shift, _) = field.lattice();
        let [_, sy, sz] = g.strides();
        let h = g.halo as isize;
        let c = [
            (p % sy) as isize - h,
            ((p / sy) % g.padded()[1]) as isize - h,
            (p / sz) as isize - h,
        ];
        (0..3)
            .map(|a| self.share(a, c[a] as f64 + shift[a]))
            .product()
    }

    /// The nodes whose gradient reads only values strictly inside the box, for a stencil of
    /// half-length `reach`: `reach` + 1 nodes in from each face, since the shear strain rates
    /// half a node out from a node read `reach` - ½ further.
    pub fn region(&self, reach: usize) -> [std::ops::Range<usize>; 3] {
        std::array::from_fn(|a| {
            let start = self.lo[a].map_or(0, |lo| lo + reach + 1);
            start..self.hi[a].saturating_sub(reach).max(start)
        })
    }
}

/// The values a [`Surface`] records per step and the bytes its operators hold, counted
/// without building it.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct SurfaceSize {
    /// Values in one record: six per face node.
    pub values: usize,
    /// Bytes of its read and injection operators.
    pub operator_bytes: usize,
}

/// One field recorded on one face, and where the record goes back in.
#[derive(Clone, Debug)]
struct Part {
    field: Field,
    /// Padded-index offsets along the normal and weights that read the field onto the face's
    /// node plane from a point there: one term on its own lattice, else the stencil's
    /// interpolation from 2L values.
    read: Vec<(isize, f32)>,
    /// Padded index of each point on the face plane (the field's own lattice along the face).
    points: Vec<usize>,
    targets: Vec<Target>,
}

/// Where a recorded field is injected: a field, offsets along the normal from each point, and
/// a coefficient per point and offset holding the material, the outward normal's sign, 1/h,
/// the spreading weight, the share along the face and the time reversal's sign.
#[derive(Clone, Debug)]
struct Target {
    field: Field,
    offsets: Vec<isize>,
    coefs: Vec<f32>,
}

/// What a part is, before it is built: its field, its face and its points along the face.
struct Layout {
    field: Field,
    /// The face's normal axis, node plane and outward sign.
    axis: usize,
    plane: usize,
    sign: f32,
    /// The two axes along the face, and the range of the field's indices along each.
    along: [(usize, std::ops::Range<usize>); 2],
}

/// The closed surface around [`Bounds`], with its recording and injection operators.
#[derive(Clone, Debug)]
pub struct Surface {
    pub bounds: Bounds,
    parts: Vec<Part>,
    len: usize,
}

/// The normal stress, the shear stress in plane (a, b) and the velocity along an axis.
fn normal_stress(a: usize) -> Field {
    [Field::Sxx, Field::Syy, Field::Szz][a]
}

fn shear_stress(a: usize, b: usize) -> Field {
    match (a.min(b), a.max(b)) {
        (0, 1) => Field::Sxy,
        (0, 2) => Field::Sxz,
        _ => Field::Syz,
    }
}

fn velocity(a: usize) -> Field {
    [Field::Vx, Field::Vy, Field::Vz][a]
}

/// The parts of every face of `bounds`, in the order of a record.
fn layouts(bounds: &Bounds) -> Vec<Layout> {
    let mut out = Vec::new();
    for a in 0..3 {
        let [b, c] = match a {
            0 => [1, 2],
            1 => [0, 2],
            _ => [0, 1],
        };
        let faces = [
            bounds.lo[a].map(|p| (p, -1.0f32)),
            Some((bounds.hi[a], 1.0)),
        ];
        for (plane, sign) in faces.into_iter().flatten() {
            let fields = [
                normal_stress(a),
                shear_stress(a, b),
                shear_stress(a, c),
                velocity(a),
                velocity(b),
                velocity(c),
            ];
            for field in fields {
                let (shift, _) = field.lattice();
                let range = |axis: usize| {
                    // Indices q whose position q + shift lies in the box, faces included.
                    let start = bounds.lo[axis].unwrap_or(0);
                    let end = if shift[axis] == 0.0 {
                        bounds.hi[axis] + 1
                    } else {
                        bounds.hi[axis]
                    };
                    (axis, start..end)
                };
                out.push(Layout {
                    field,
                    axis: a,
                    plane,
                    sign,
                    along: [range(b), range(c)],
                });
            }
        }
    }
    out
}

impl Surface {
    /// The surface `inset` nodes inside the absorbing layers of `solver` (by default, and at
    /// least, its stencil's half-length L), with operators for its model and stencil.
    pub fn new(solver: &Solver, inset: Option<usize>) -> Result<Self> {
        let reach = solver.stencil.half_length();
        let bounds = Bounds::new(
            &solver.interior(),
            solver.free_surface,
            reach,
            inset.unwrap_or(reach),
        )?;
        let g = solver.grid;
        let strides = g.strides();
        let m = &solver.material;
        let mut parts = Vec::new();
        for l in layouts(&bounds) {
            let a = l.axis;
            let (shift, _) = l.field.lattice();
            // Along the normal: on the field's lattice, or read with the stencil's
            // interpolation from the half nodes around the plane.
            let stagger = |parity: f32| -> Vec<(isize, f32)> {
                solver
                    .stagger_weights(a, l.plane, parity)
                    .into_iter()
                    .map(|(m, w)| ((m as isize - l.plane as isize) * strides[a] as isize, w))
                    .collect()
            };
            let parity_along = |f: Field| if a == 2 { f.lattice().1 } else { 1.0 };
            let read = if shift[a] == 0.0 {
                vec![(0, 1.0)]
            } else {
                stagger(parity_along(l.field))
            };
            let [(b, rb), (c, rc)] = &l.along;
            let mut points = Vec::with_capacity(rb.len() * rc.len());
            let mut shares = Vec::with_capacity(points.capacity());
            for qc in rc.clone() {
                for qb in rb.clone() {
                    let mut node = [0; 3];
                    node[a] = l.plane;
                    node[*b] = qb;
                    node[*c] = qc;
                    points.push(g.idx(node[0], node[1], node[2]));
                    shares.push(
                        bounds.share(*b, qb as f64 + shift[*b])
                            * bounds.share(*c, qc as f64 + shift[*c]),
                    );
                }
            }
            let scale = l.sign / g.h[a] as f32;
            let on_surface =
                |p: usize| solver.free_surface && p / g.plane_len() == g.halo && a != 2;
            // A target spread along the normal from each point, or at the point itself.
            let target = |field: Field, spread: bool, coef: &dyn Fn(usize, f32) -> f32| {
                let offsets: Vec<(isize, f32)> = if spread {
                    stagger(parity_along(field))
                } else {
                    vec![(0, 1.0)]
                };
                let mut coefs = Vec::with_capacity(points.len() * offsets.len());
                for (&p, &share) in points.iter().zip(&shares) {
                    for &(o, w) in &offsets {
                        let q = (p as isize + o) as usize;
                        coefs.push(coef(q, scale * share * w));
                    }
                }
                Target {
                    field,
                    offsets: offsets.into_iter().map(|(o, _)| o).collect(),
                    coefs,
                }
            };
            let targets = match l.field {
                // The force n·σ, injected into the velocities with the reversed sign of σ.
                f if f == normal_stress(a) => {
                    let bu = m.buoyancy(a);
                    vec![target(velocity(a), true, &|q, s| -s * bu[q])]
                }
                f if f == shear_stress(a, *b) || f == shear_stress(a, *c) => {
                    let other = if f == shear_stress(a, *b) { *b } else { *c };
                    let bu = m.buoyancy(other);
                    vec![target(velocity(other), false, &|q, s| -s * bu[q])]
                }
                // The stress source c : (v ⊗ n): the normal velocity into the normal stresses,
                // on a free surface with its plane-stress moduli (σzz stays zero there).
                f if f == velocity(a) => (0..3)
                    .map(|d| {
                        target(normal_stress(d), false, &|q, s| {
                            let (lam, l2m) = (m.lam[q], m.l2m[q]);
                            if on_surface(q) {
                                let r = if l2m > 0.0 { lam * lam / l2m } else { 0.0 };
                                s * match d {
                                    2 => 0.0,
                                    _ if d == a => l2m - r,
                                    _ => lam - r,
                                }
                            } else {
                                s * if d == a { l2m } else { lam }
                            }
                        })
                    })
                    .collect(),
                // A tangential velocity into the shear stress in its plane with the normal.
                f => {
                    let other = if f == velocity(*b) { *b } else { *c };
                    let shear = shear_stress(a, other);
                    let mu = match shear {
                        Field::Sxy => &m.mxy,
                        Field::Sxz => &m.mxz,
                        _ => &m.myz,
                    };
                    vec![target(shear, true, &|q, s| s * mu[q])]
                }
            };
            parts.push(Part {
                field: l.field,
                read,
                points,
                targets,
            });
        }
        let len = parts.iter().map(|p| p.points.len()).sum();
        Ok(Self { bounds, parts, len })
    }

    /// [`Surface::new`]'s size for a stencil of half-length `reach`, counted.
    pub fn size(
        interior: &[std::ops::Range<usize>; 3],
        free_surface: bool,
        reach: usize,
        inset: Option<usize>,
    ) -> Result<SurfaceSize> {
        let bounds = Bounds::new(interior, free_surface, reach, inset.unwrap_or(reach))?;
        let (mut values, mut bytes) = (0, 0);
        let spread = 2 * reach;
        for l in layouts(&bounds) {
            let points = l.along[0].1.len() * l.along[1].1.len();
            let a = l.axis;
            let staggered = l.field.lattice().0[a] != 0.0;
            let read = if staggered { spread } else { 1 };
            // Offsets per target: spread along the normal or not.
            let offsets: &[usize] = match l.field {
                f if f == normal_stress(a) => &[spread],
                f if f == velocity(a) => &[1, 1, 1],
                Field::Vx | Field::Vy | Field::Vz => &[spread],
                _ => &[1],
            };
            values += points;
            bytes += read * std::mem::size_of::<(isize, f32)>()
                + points * std::mem::size_of::<usize>()
                + offsets
                    .iter()
                    .map(|&o| o * std::mem::size_of::<isize>() + points * o * 4)
                    .sum::<usize>();
        }
        Ok(SurfaceSize {
            values,
            operator_bytes: bytes,
        })
    }

    /// Values in one record.
    pub fn len(&self) -> usize {
        self.len
    }

    /// Whether a record holds nothing (never, for a surface [`Surface::new`] built).
    pub fn is_empty(&self) -> bool {
        self.len == 0
    }

    /// Bytes its operators hold.
    pub fn memory_bytes(&self) -> usize {
        self.parts
            .iter()
            .map(|p| {
                std::mem::size_of_val(p.read.as_slice())
                    + std::mem::size_of_val(p.points.as_slice())
                    + p.targets
                        .iter()
                        .map(|t| {
                            std::mem::size_of_val(t.offsets.as_slice())
                                + std::mem::size_of_val(t.coefs.as_slice())
                        })
                        .sum::<usize>()
            })
            .sum()
    }

    /// Records the surface's values of `fields` into `out`, one record.
    pub fn save(&self, fields: &Fields, out: &mut [f32]) {
        assert_eq!(out.len(), self.len);
        let mut o = 0;
        for part in &self.parts {
            let v = fields.get(part.field);
            for &p in &part.points {
                out[o] = part
                    .read
                    .iter()
                    .map(|&(off, w)| w * v[(p as isize + off) as usize])
                    .sum();
                o += 1;
            }
        }
    }

    /// Adds the force n·σ of `record` to the velocities, as their update's commutator with
    /// the box would: call it after the velocity update of a reversed step.
    pub fn inject_forces(&self, fields: &mut Fields, record: &[f32]) {
        self.inject(fields, record, |f| !f.is_velocity());
    }

    /// Adds the stress source c : (v ⊗ n) of `record` to the stresses: call it after the
    /// stress update of a reversed step.
    pub fn inject_stresses(&self, fields: &mut Fields, record: &[f32]) {
        self.inject(fields, record, Field::is_velocity);
    }

    fn inject(&self, fields: &mut Fields, record: &[f32], which: impl Fn(Field) -> bool) {
        assert_eq!(record.len(), self.len);
        let mut o = 0;
        for part in &self.parts {
            let n = part.points.len();
            if which(part.field) {
                let values = &record[o..o + n];
                for t in &part.targets {
                    let out = fields.get_mut(t.field);
                    let k = t.offsets.len();
                    for (i, (&p, &value)) in part.points.iter().zip(values).enumerate() {
                        for (&off, &c) in t.offsets.iter().zip(&t.coefs[i * k..(i + 1) * k]) {
                            out[(p as isize + off) as usize] += c * value;
                        }
                    }
                }
            }
            o += n;
        }
    }

    /// Turns the solver's state into the start of the reversed run inside the box: every
    /// field times its share inside ([`Bounds::field_share`]), the stresses negated, and the
    /// absorbing layers' memory cleared.
    pub fn truncate(&self, solver: &mut Solver) {
        let g = solver.grid;
        let [px, py, pz] = g.padded();
        let h = g.halo as isize;
        for field in Field::ALL {
            let (shift, _) = field.lattice();
            let sign = if field.is_velocity() { 1.0 } else { -1.0 };
            let along = |a: usize, n: usize| -> Vec<f32> {
                (0..n)
                    .map(|q| self.bounds.share(a, (q as isize - h) as f64 + shift[a]))
                    .collect()
            };
            let (mx, my, mz) = (along(0, px), along(1, py), along(2, pz));
            let v = solver.fields.get_mut(field);
            for k in 0..pz {
                for j in 0..py {
                    let s = sign * my[j] * mz[k];
                    let row = &mut v[(k * py + j) * px..(k * py + j + 1) * px];
                    for (x, &w) in row.iter_mut().zip(&mx) {
                        *x *= s * w;
                    }
                }
            }
        }
        solver.clear_layers();
    }

    /// `op` (a source's) with each of its terms times its field's share inside the box, as
    /// the truncated field feels it.
    pub fn truncated_op(&self, g: &Grid, op: &PointOp) -> PointOp {
        let mut op = op.clone();
        for (f, p, w) in &mut op.terms {
            *w *= self.bounds.field_share(g, *f, *p);
        }
        op
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// The box sits `inset` in from the layers, open along z under a free surface, and refuses
    /// an inset shorter than the stencil's reach or a box with no room inside.
    #[test]
    fn the_box_sits_inside_the_layers() {
        let interior = [8..40, 8..30, 0..24];
        let b = Bounds::new(&interior, true, 4, 4).unwrap();
        assert_eq!(b.lo, [Some(12), Some(12), None]);
        assert_eq!(b.hi, [35, 25, 19]);
        let closed = Bounds::new(&interior, false, 4, 5).unwrap();
        assert_eq!(closed.lo, [Some(13), Some(13), Some(5)]);
        assert!(Bounds::new(&interior, true, 4, 3).is_err());
        assert!(Bounds::new(&[8..18, 8..30, 0..24], true, 4, 4).is_err());
        assert!(Bounds::new(&[8..19, 8..30, 0..24], true, 4, 4).is_ok());
        // The gradient's region: L + 1 in from each face, from the top under a free surface.
        assert_eq!(b.region(4), [17..31, 17..21, 0..15]);
    }

    /// A value's share is 1 inside, ½ on a face and 0 outside; above a free surface (and on
    /// it) it is inside.
    #[test]
    fn shares_are_halves_on_the_faces() {
        let b = Bounds::new(&[8..40, 8..30, 0..24], true, 4, 4).unwrap();
        let shares: Vec<f32> = [11.5, 12.0, 12.5, 30.0, 35.0, 35.5]
            .iter()
            .map(|&x| b.share(0, x))
            .collect();
        assert_eq!(shares, [0.0, 0.5, 1.0, 1.0, 0.5, 0.0]);
        assert_eq!(b.share(2, -0.5), 1.0);
        assert_eq!(b.share(2, 0.0), 1.0);
        assert_eq!(b.share(2, 19.0), 0.5);
    }
}
