//! Which compiled version of the update loops a solver runs: the plain loops, compiled for the
//! build's target (SSE2 on x86-64 by default, NEON on Apple silicon, SIMD128 in the browser),
//! or on x86-64 the same loops compiled again with AVX2 enabled, chosen at run time where the
//! CPU has it. The default build therefore runs on any x86-64 CPU.
//!
//! The AVX2 loops are the same source, so they do the same arithmetic in the same order, only
//! eight values at a time instead of four: the compiler vectorises them across the nodes of a
//! row, never reassociates floating-point sums, and Rust never fuses a·b + c into one rounding
//! (FMA is not enabled, and wouldn't be used for it if it were). Both give the same results
//! bit for bit; a test checks it.
//!
//! `VIBE_RATIONS_SIMD=off` runs the plain loops, for comparisons and for machines or virtual
//! machines where AVX2 misbehaves. In a build for the host's CPU (`-C target-cpu=native`) the
//! plain loops are themselves compiled with its features, so running without AVX2 means the
//! default build.
//!
//! A third version fuses: where the loops add a product to something, a·b + c, it rounds the
//! result once ([`Fused`], `f32::mul_add`) instead of the product and then the sum
//! ([`Separate`]). That is one instruction instead of two on CPUs that have it (x86-64 with
//! FMA beside AVX2, any 64-bit ARM), but other bits: the results are close to the others',
//! not the same. So it runs only when asked for, `VIBE_RATIONS_SIMD=fma` or `run --fma`.

use std::sync::OnceLock;

use serde::Serialize;

/// The environment variable that turns the AVX2 loops off (`off`), read once per process.
pub const SIMD_VARIABLE: &str = "VIBE_RATIONS_SIMD";

/// A compiled version of the update loops.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Serialize)]
#[serde(rename_all = "lowercase")]
pub enum Kernels {
    /// Compiled for the build's target.
    Plain,
    /// Compiled with AVX2, on x86-64 CPUs that have it.
    Avx2,
    /// With fused multiply-adds: compiled with AVX2 and FMA on x86-64 CPUs that have both,
    /// for the build's target on 64-bit ARM. Not the others' results bit for bit.
    Fma,
}

impl Kernels {
    /// Whether this CPU can run them.
    pub fn available(self) -> bool {
        match self {
            Kernels::Plain => true,
            Kernels::Avx2 => has_avx2(),
            Kernels::Fma => has_fma(),
        }
    }

    /// What a solver runs unless told otherwise: AVX2 where the CPU has it, unless
    /// [`SIMD_VARIABLE`] says `off`.
    pub fn chosen() -> Kernels {
        choice().0
    }

    /// A warning about [`SIMD_VARIABLE`], if it holds something it can't do.
    pub fn warning() -> Option<&'static str> {
        choice().1.as_deref()
    }

    /// `self` where this CPU can run it, else the plain loops.
    pub fn or_plain(self) -> Kernels {
        if self.available() {
            self
        } else {
            Kernels::Plain
        }
    }

    pub fn name(self) -> &'static str {
        match self {
            Kernels::Plain => "plain",
            Kernels::Avx2 => "avx2",
            Kernels::Fma => "fma",
        }
    }
}

impl std::fmt::Display for Kernels {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.write_str(self.name())
    }
}

/// The kernels to run and a warning about the variable, decided once.
fn choice() -> &'static (Kernels, Option<String>) {
    static CHOICE: OnceLock<(Kernels, Option<String>)> = OnceLock::new();
    CHOICE.get_or_init(|| decide(std::env::var(SIMD_VARIABLE).ok().as_deref()))
}

/// The kernels a value of [`SIMD_VARIABLE`] asks for, given this CPU: `off` (or `plain`) the
/// plain loops; unset, empty, `auto` or `avx2` the AVX2 ones where the CPU has them; `fma`
/// the fused ones. Anything else is ignored with a warning, as are `avx2` and `fma` where the
/// CPU can't run them.
fn decide(value: Option<&str>) -> (Kernels, Option<String>) {
    let best = Kernels::Avx2.or_plain();
    match value.map(str::trim) {
        None | Some("" | "auto") => (best, None),
        Some("off" | "plain") => (Kernels::Plain, None),
        Some("avx2") if best == Kernels::Avx2 => (best, None),
        Some("avx2") => (
            best,
            Some(format!(
                "{SIMD_VARIABLE}=avx2, but this CPU can't run AVX2: running the plain loops"
            )),
        ),
        Some("fma") if Kernels::Fma.available() => (Kernels::Fma, None),
        Some("fma") => (
            best,
            Some(format!(
                "{SIMD_VARIABLE}=fma, but this CPU has no fused multiply-add the loops can use: \
                 running the {best} loops"
            )),
        ),
        Some(other) => (
            best,
            Some(format!(
                "{SIMD_VARIABLE}={other} means nothing (it takes off, avx2 or fma): running the \
                 {best} loops"
            )),
        ),
    }
}

fn has_avx2() -> bool {
    #[cfg(target_arch = "x86_64")]
    {
        std::arch::is_x86_feature_detected!("avx2")
    }
    #[cfg(not(target_arch = "x86_64"))]
    {
        false
    }
}

/// Whether the fused loops have an instruction to run on: FMA beside AVX2 on x86-64 (every CPU
/// with one has had the other, bar a few), always on 64-bit ARM, and nowhere else (without the
/// instruction a fused product is a call into a maths library for every value).
fn has_fma() -> bool {
    #[cfg(target_arch = "x86_64")]
    {
        std::arch::is_x86_feature_detected!("avx2") && std::arch::is_x86_feature_detected!("fma")
    }
    #[cfg(not(target_arch = "x86_64"))]
    {
        cfg!(target_arch = "aarch64")
    }
}

/// How the update loops round a product added to something, a·b + c: the one place where the
/// compiled versions' arithmetic may differ. The loops take it as a type, so each version is
/// compiled with its own.
pub trait Arithmetic {
    /// a·b + c.
    fn product_plus(a: f32, b: f32, c: f32) -> f32;
}

/// The product rounded, then the sum: the reference, the same on every machine.
pub struct Separate;

/// Product and sum rounded once (a fused multiply-add): half a unit in the last place closer
/// to a·b + c each time, and so not [`Separate`]'s bits.
pub struct Fused;

impl Arithmetic for Separate {
    #[inline(always)]
    fn product_plus(a: f32, b: f32, c: f32) -> f32 {
        a * b + c
    }
}

impl Arithmetic for Fused {
    #[inline(always)]
    fn product_plus(a: f32, b: f32, c: f32) -> f32 {
        a.mul_add(b, c)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn the_variable_turns_avx2_off_and_nothing_else_on() {
        let best = Kernels::Avx2.or_plain();
        assert_eq!(decide(None), (best, None));
        assert_eq!(decide(Some("")), (best, None));
        assert_eq!(decide(Some("auto")), (best, None));
        assert_eq!(decide(Some("off")), (Kernels::Plain, None));
        assert_eq!(decide(Some(" plain ")), (Kernels::Plain, None));
        let (k, warning) = decide(Some("avx2"));
        assert_eq!(k, best);
        assert_eq!(warning.is_some(), best == Kernels::Plain);
        let (k, warning) = decide(Some("sse"));
        assert_eq!(k, best);
        assert!(warning.unwrap().contains("sse"));
        // The fused loops only when asked for, and only where they have an instruction.
        let (k, warning) = decide(Some("fma"));
        if Kernels::Fma.available() {
            assert_eq!((k, warning), (Kernels::Fma, None));
        } else {
            assert_eq!(k, best);
            assert!(warning.unwrap().contains("fma"));
        }
    }

    /// A fused product is the exact a·b + c rounded once: where the product alone isn't a
    /// single-precision number, the two differ.
    #[test]
    fn fused_products_round_once() {
        let (a, b, c) = (1.0 + f32::EPSILON, 1.0 + f32::EPSILON, -1.0);
        // (1 + ε)² = 1 + 2ε + ε²: the product drops ε², the fused sum keeps it.
        assert_eq!(Separate::product_plus(a, b, c), 2.0 * f32::EPSILON);
        assert_eq!(
            Fused::product_plus(a, b, c),
            2.0 * f32::EPSILON + f32::EPSILON * f32::EPSILON
        );
        assert_eq!(Separate::product_plus(3.0, 4.0, 5.0), 17.0);
        assert_eq!(Fused::product_plus(3.0, 4.0, 5.0), 17.0);
    }

    #[test]
    fn only_x86_64_runs_avx2() {
        assert!(Kernels::Plain.available());
        if !cfg!(target_arch = "x86_64") {
            assert!(!Kernels::Avx2.available());
            assert_eq!(Kernels::Avx2.or_plain(), Kernels::Plain);
            assert_eq!(Kernels::Fma.available(), cfg!(target_arch = "aarch64"));
        }
    }
}
