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
}

impl Kernels {
    /// Whether this CPU can run them.
    pub fn available(self) -> bool {
        match self {
            Kernels::Plain => true,
            Kernels::Avx2 => has_avx2(),
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
/// plain loops; unset, empty, `auto` or `avx2` the AVX2 ones where the CPU has them. Anything
/// else is ignored with a warning, as is `avx2` where the CPU hasn't got it.
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
        Some(other) => (
            best,
            Some(format!(
                "{SIMD_VARIABLE}={other} means nothing (it takes off or avx2): running the {best} loops"
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
    }

    #[test]
    fn only_x86_64_runs_avx2() {
        assert!(Kernels::Plain.available());
        if !cfg!(target_arch = "x86_64") {
            assert!(!Kernels::Avx2.available());
            assert_eq!(Kernels::Avx2.or_plain(), Kernels::Plain);
        }
    }
}
