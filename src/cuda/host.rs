//! The GPU's kernels on a CPU: the same source compiled by the system's C compiler and run a
//! thread at a time (`host.c`), slowly. It is what lets the kernels, and everything that
//! drives them, be tested against the solver where there is no GPU.

use std::ffi::c_uint;
use std::process::Command;
use std::sync::atomic::{AtomicUsize, Ordering};

use anyhow::{Context, Result, ensure};

use super::driver::Library;
use super::{Args, Device, KERNELS, Kernel};

/// What runs every thread of a launch in turn: `emulate_<kernel>` in `host.c`.
type Emulate = unsafe extern "C" fn(*const Args, *const c_uint, *const c_uint);

const SHIM: &str = include_str!("host.c");

/// A CPU standing in for a GPU.
pub struct Host {
    functions: [Emulate; Kernel::ALL.len()],
    /// What was allocated, in eight-byte words so that anything is aligned in it.
    memory: Vec<Box<[u64]>>,
}

impl Host {
    /// Compiles the kernels with the system's C compiler (`cc`, or what `CC` names) and loads
    /// them. Without fused multiply-adds, as the solver is compiled.
    pub fn build() -> Result<Self> {
        static BUILDS: AtomicUsize = AtomicUsize::new(0);
        let dir = std::env::temp_dir().join(format!(
            "vibe-rations-kernels-{}-{}",
            std::process::id(),
            BUILDS.fetch_add(1, Ordering::Relaxed)
        ));
        std::fs::create_dir_all(&dir).with_context(|| format!("creating {}", dir.display()))?;
        let built = (|| {
            std::fs::write(dir.join("host.c"), SHIM)?;
            std::fs::write(dir.join("kernels.cu"), KERNELS)?;
            let library = dir.join("kernels.so");
            let compiler = std::env::var("CC").unwrap_or_else(|_| "cc".into());
            let output = Command::new(&compiler)
                .args(["-O2", "-ffp-contract=off", "-shared", "-fPIC", "-o"])
                .arg(&library)
                .arg(dir.join("host.c"))
                .output()
                .with_context(|| format!("running the C compiler, {compiler}"))?;
            ensure!(
                output.status.success(),
                "{compiler} did not compile the kernels: {}",
                String::from_utf8_lossy(&output.stderr)
            );
            let library = Library::open(&[library.to_str().context("a path that isn't text")?])?;
            let mut functions = Vec::new();
            for kernel in Kernel::ALL {
                // SAFETY: host.c defines each with this signature.
                functions.push(unsafe { library.function(&format!("emulate_{}", kernel.name()))? });
            }
            Ok(functions)
        })();
        // Once loaded the library needs its file no longer.
        let _ = std::fs::remove_dir_all(&dir);
        let functions: Vec<Emulate> = built?;
        Ok(Self {
            functions: functions.try_into().expect("a function per kernel"),
            memory: Vec::new(),
        })
    }

    /// The allocation `bytes` bytes at `at` lie in, checked.
    fn region(&mut self, at: u64, bytes: usize) -> Result<*mut u8> {
        let found = self.memory.iter_mut().find_map(|m| {
            let start = m.as_ptr() as u64;
            let room = (8 * m.len()) as u64;
            (at >= start && at - start + bytes as u64 <= room).then(|| {
                // SAFETY: inside the allocation, by the test above.
                unsafe { m.as_mut_ptr().cast::<u8>().add((at - start) as usize) }
            })
        });
        found.with_context(|| format!("{bytes} bytes at {at:#x} are not in one allocation"))
    }
}

impl Device for Host {
    fn name(&self) -> String {
        "this CPU, running the GPU's kernels a thread at a time".into()
    }

    fn prepare(&mut self, _half_length: usize) -> Result<()> {
        // Compiled for any stencil: the kernels read its length from their arguments.
        Ok(())
    }

    fn free_bytes(&mut self) -> Result<Option<usize>> {
        Ok(None)
    }

    fn alloc(&mut self, bytes: usize) -> Result<u64> {
        let words = vec![0u64; bytes.div_ceil(8).max(1)].into_boxed_slice();
        let at = words.as_ptr() as u64;
        self.memory.push(words);
        Ok(at)
    }

    fn upload(&mut self, to: u64, data: &[u8]) -> Result<()> {
        let p = self.region(to, data.len())?;
        // SAFETY: `region` checked the room; the allocations are this device's own.
        unsafe { std::ptr::copy_nonoverlapping(data.as_ptr(), p, data.len()) };
        Ok(())
    }

    fn download(&mut self, from: u64, out: &mut [u8]) -> Result<()> {
        let p = self.region(from, out.len())?;
        // SAFETY: as in `upload`.
        unsafe { std::ptr::copy_nonoverlapping(p, out.as_mut_ptr(), out.len()) };
        Ok(())
    }

    fn launch(
        &mut self,
        kernel: Kernel,
        grid: [u32; 3],
        block: [u32; 3],
        args: &Args,
    ) -> Result<()> {
        // SAFETY: the kernels read and write the arrays `args` names, which the solver that
        // launches them allocated here with the sizes they index; nothing else holds them
        // meanwhile.
        unsafe {
            (self.functions[kernel as usize])(args, grid.as_ptr(), block.as_ptr());
        }
        Ok(())
    }

    fn wait(&mut self) -> Result<()> {
        Ok(())
    }
}
