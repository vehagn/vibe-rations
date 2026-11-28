//! An explosion in a homogeneous full space recorded as pressure, and the change of that
//! record for a change of the model at one node, to first order (Born): the closed forms the
//! Hessian's column study measures the Born run and the Gauss–Newton Hessian against.
//!
//! An explosion of moment M(t) radiates a spherical P wave, u = ∇φ with φ = -M(t - r/α)/(c r)
//! and c = 4πρα² (Aki & Richards 2002, §4.3: the radial displacement is
//! [M/r² + Ṁ/(α r)]/c). With γ the unit vector from the source, its strain is
//!
//! ε = Σₖ aₖ M⁽ᵏ⁾(t - r/α),  a₀ = (I - 3γγ)/(c r³),  a₁ = (I - 3γγ)/(c α r²),  a₂ = -γγ/(c α² r),
//!
//! whose trace is a₂'s alone, and the pressure -κ ∇·u it carries is κ M̈(t - r/α)/(c α² r),
//! κ = λ + 2μ/3.
//!
//! A change δρ, δλ, δμ in a volume V at y scatters that wave as the force -δρ ü V and the
//! stress δσ = δλ tr(ε) I + 2 δμ ε over V. The pressure a receiver records from them is, by
//! reciprocity, what the explosion of impulsive moment *at the receiver* does at y: with E
//! its displacement there and ε_E its strain,
//!
//! δp = κ V (δρ E·ü + ε_E : δσ), the products convolutions in time.
//!
//! Both ends being retarded impulses' derivatives, each convolution is a derivative of M
//! delayed by (r + R)/α, R the receiver's distance: a scatterer's share is Σₙ cₙ M⁽ⁿ⁾, n ≤ 4.
//!
//! The staggered grid spreads a node's change as its material averaging does ([`scattering`]):
//! δλ, and δμ for the normal strains, stay on the node; δμ for each shear strain goes a quarter
//! to each of the four edge midpoints around the node in that strain's plane (the harmonic
//! mean of four equal moduli); δρ goes half to each of the two nearest positions of each
//! velocity component, which carries that component of the force.

use std::f64::consts::PI;

use super::stokes::{Medium, Ricker};
use super::{norm, sub};

/// An explosion whose moment *rate* is `rate`·w(t), w a Ricker wavelet of peak frequency `f0`
/// centred at `t0`, as the solver's `explosion` source with a Ricker wavelet is.
#[derive(Clone, Copy, Debug)]
pub struct Explosion {
    pub rate: f64,
    pub f0: f64,
    pub t0: f64,
}

impl Explosion {
    /// The moment's derivatives M⁽ⁿ⁾(t), n = 0 to 4.
    pub fn moment(&self, t: f64) -> [f64; 5] {
        let w = Ricker {
            f0: self.f0,
            t0: self.t0,
            amplitude: 1.0,
            level: 0,
        };
        std::array::from_fn(|n| self.rate * w.derivative(n as i32 - 1, t))
    }

    /// The pressure at distance `r` and time `t`.
    pub fn pressure(&self, m: &Medium, r: f64, t: f64) -> f64 {
        let kappa = m.rho * (m.vp * m.vp - 4.0 / 3.0 * m.vs * m.vs);
        kappa * self.moment(t - r / m.vp)[2] / (4.0 * PI * m.rho * m.vp.powi(4) * r)
    }
}

/// One point scatterer's share of the record per unit moment: Σₙ `c[n]` M⁽ⁿ⁾(t - `delay`).
#[derive(Clone, Copy, Debug)]
pub struct Term {
    pub delay: f64,
    pub c: [f64; 5],
}

/// The record's change at time `t` that `terms` add up to.
pub fn evaluate(terms: &[Term], explosion: &Explosion, t: f64) -> f64 {
    terms
        .iter()
        .map(|term| {
            let m = explosion.moment(t - term.delay);
            term.c.iter().zip(m).map(|(c, m)| c * m).sum::<f64>()
        })
        .sum()
}

/// The strain's coefficients a₀, a₁, a₂.
type Strain = [[[f64; 3]; 3]; 3];
/// The displacement's coefficients of M and Ṁ along each axis.
type Displacement = [[f64; 3]; 2];

/// aₖ of the module's text for a point at offset `d` from the explosion, the displacement's
/// two coefficients along each axis, of M and Ṁ: γ/(c r²) and γ/(c α r), and the distance.
fn field(m: &Medium, d: [f64; 3]) -> (Strain, Displacement, f64) {
    let r = norm(d);
    let g = d.map(|x| x / r);
    let c = 4.0 * PI * m.rho * m.vp * m.vp;
    let mut a = [[[0.0; 3]; 3]; 3];
    for i in 0..3 {
        for j in 0..3 {
            let delta = if i == j { 1.0 } else { 0.0 };
            let gg = g[i] * g[j];
            a[0][i][j] = (delta - 3.0 * gg) / (c * r.powi(3));
            a[1][i][j] = (delta - 3.0 * gg) / (c * m.vp * r * r);
            a[2][i][j] = -gg / (c * m.vp * m.vp * r);
        }
    }
    let u = [g.map(|x| x / (c * r * r)), g.map(|x| x / (c * m.vp * r))];
    (a, u, r)
}

/// The pressure's change at `receiver`, per unit moment of an explosion at `source`, for unit
/// changes of ρ, of λ and of μ (in that order) at the node at `node` of a grid spaced `h`,
/// sampled as the staggered grid samples a node's change (the module's text).
pub fn scattering(
    m: &Medium,
    source: [f64; 3],
    receiver: [f64; 3],
    node: [f64; 3],
    h: f64,
) -> [Vec<Term>; 3] {
    let kappa = m.rho * (m.vp * m.vp - 4.0 / 3.0 * m.vs * m.vs);
    let scale = kappa * h.powi(3);
    let at = |shift: [f64; 3]| {
        let y: [f64; 3] = std::array::from_fn(|a| node[a] + 0.5 * h * shift[a]);
        let (a, u, r) = field(m, sub(y, source));
        let (b, e, big_r) = field(m, sub(y, receiver));
        (a, u, b, e, (r + big_r) / m.vp)
    };
    let mut rho = Vec::with_capacity(6);
    let mut mu = Vec::with_capacity(13);

    // On the node: λ takes the strains' traces, a₂'s and b₂'s alone; μ the normal strains.
    let (a, _, b, _, delay) = at([0.0; 3]);
    let trace = |t: &[[f64; 3]; 3]| t[0][0] + t[1][1] + t[2][2];
    let mut c = [0.0; 5];
    c[4] = scale * trace(&b[2]) * trace(&a[2]);
    let lambda = vec![Term { delay, c }];
    let mut c = [0.0; 5];
    for k in 0..3 {
        for l in 0..3 {
            c[k + l] += scale * 2.0 * (0..3).map(|i| b[k][i][i] * a[l][i][i]).sum::<f64>();
        }
    }
    mu.push(Term { delay, c });

    // The edges: 2 δμ (ε_E,ij ε_ij + ε_E,ji ε_ji) with a quarter of δμ on each of four.
    for (i, j) in [(0, 1), (0, 2), (1, 2)] {
        for (si, sj) in [(-1.0, -1.0), (1.0, -1.0), (-1.0, 1.0), (1.0, 1.0)] {
            let mut shift = [0.0; 3];
            shift[i] = si;
            shift[j] = sj;
            let (a, _, b, _, delay) = at(shift);
            let mut c = [0.0; 5];
            for k in 0..3 {
                for l in 0..3 {
                    c[k + l] += scale * b[k][i][j] * a[l][i][j];
                }
            }
            mu.push(Term { delay, c });
        }
    }

    // The velocity positions: E_i ü_i with half of δρ on each of two per component; E's
    // coefficients are of the impulse and its derivative, ü's of M̈ and M⃛.
    for i in 0..3 {
        for s in [-1.0, 1.0] {
            let mut shift = [0.0; 3];
            shift[i] = s;
            let (_, u, _, e, delay) = at(shift);
            let mut c = [0.0; 5];
            for k in 0..2 {
                for l in 0..2 {
                    c[k + l + 2] += scale * 0.5 * e[k][i] * u[l][i];
                }
            }
            rho.push(Term { delay, c });
        }
    }
    [rho, lambda, mu]
}

/// [`scattering`]'s terms for unit changes of vp, of vs and of ρ (in that order) at the node,
/// the other two held, in a medium `m`: δλ = 2ρα δα - 4ρβ δβ + (α² - 2β²) δρ and
/// δμ = 2ρβ δβ + β² δρ.
pub fn velocity_scattering(
    m: &Medium,
    source: [f64; 3],
    receiver: [f64; 3],
    node: [f64; 3],
    h: f64,
) -> [Vec<Term>; 3] {
    let [rho, lambda, mu] = scattering(m, source, receiver, node, h);
    fn scaled(terms: &[Term], f: f64) -> impl Iterator<Item = Term> + '_ {
        terms.iter().map(move |t| Term {
            delay: t.delay,
            c: t.c.map(|c| f * c),
        })
    }
    let (a, b, r) = (m.vp, m.vs, m.rho);
    let vp = scaled(&lambda, 2.0 * r * a).collect();
    let vs = scaled(&lambda, -4.0 * r * b)
        .chain(scaled(&mu, 2.0 * r * b))
        .collect();
    let density = scaled(&rho, 1.0)
        .chain(scaled(&lambda, a * a - 2.0 * b * b))
        .chain(scaled(&mu, b * b))
        .collect();
    [vp, vs, density]
}
