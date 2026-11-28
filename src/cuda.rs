//! The time step on an NVIDIA GPU.
//!
//! A [`GpuSolver`] holds a [`Solver`]'s wavefield, material and absorbing layers on a
//! [`Device`] and steps them there with the kernels of `cuda/kernels.cu`: one thread per node
//! doing what the solver's rows do, in the same order and with the same single-precision
//! operations. A GPU's additions, multiplications and divisions round as a CPU's do
//! (IEEE 754), so as long as no product is fused into a sum the fields and traces of a run on
//! the GPU are the CPU's bit for bit, and the CPU stays the reference: the tests here compare
//! bits, not tolerances.
//!
//! The whole step runs on the device (the updates, the layers, the free surface's mirrors,
//! the sponge, sources and receivers), since copying nine fields across every step would cost
//! more than stepping them where they are. What crosses is small: the sources' values before
//! the run, what the receivers' terms read in chunks, the largest velocity now and then, and
//! the fields when a snapshot or the end of the run wants them.
//!
//! Nothing is linked when this is built: the kernels are CUDA C source held in the binary,
//! compiled when a run starts by NVRTC and launched through the driver's API, both found at
//! run time ([`Cuda`]). A build is the same with or without a GPU, and the same kernels
//! compile for a CPU ([`Host`]), where a test steps them a thread at a time against the
//! solver. Forward runs ([`run`]) and gradients replayed from checkpoints ([`gradient`],
//! `cuda/replay.rs`), and so inversions; rebuilding the forward wavefield, the receivers' runs
//! and the Hessian step on the CPU.

use std::cell::RefCell;
use std::ops::Range;
use std::path::Path;
use std::rc::Rc;
use std::time::{Duration, Instant};

use anyhow::{Context, Result, bail, ensure};

use crate::cpml::Stagger;
use crate::grid::Grid;
use crate::receiver::Receivers;
use crate::simulation::{Simulation, human_time};
use crate::solver::{Field, PointOp, Solver, State};
use crate::source::Source;

mod driver;
mod host;
mod replay;

pub use driver::{Cuda, FMAD_VARIABLE, NVRTC_VARIABLE, compile};
pub use host::Host;
pub use replay::gradient;

/// The kernels' source: CUDA C, in what it shares with C.
pub const KERNELS: &str = include_str!("cuda/kernels.cu");

/// The environment variable that sets the threads along x, y and z of a block of the
/// updates, as `32,4,4`: for tuning.
pub const BLOCK_VARIABLE: &str = "VIBE_RATIONS_CUDA_BLOCK";

/// A kernel of `cuda/kernels.cu`, in the order of [`Kernel::ALL`].
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Kernel {
    Velocity,
    Stress,
    Damp,
    Mirror,
    Inject,
    Gather,
    Peak,
    Layout,
    Strain,
    Correlate,
    Copy,
}

impl Kernel {
    pub const ALL: [Kernel; 11] = [
        Kernel::Velocity,
        Kernel::Stress,
        Kernel::Damp,
        Kernel::Mirror,
        Kernel::Inject,
        Kernel::Gather,
        Kernel::Peak,
        Kernel::Layout,
        Kernel::Strain,
        Kernel::Correlate,
        Kernel::Copy,
    ];

    /// Its function's name in the source.
    pub fn name(self) -> &'static str {
        match self {
            Kernel::Velocity => "velocity",
            Kernel::Stress => "stress",
            Kernel::Damp => "damp",
            Kernel::Mirror => "mirror",
            Kernel::Inject => "inject",
            Kernel::Gather => "gather",
            Kernel::Peak => "peak",
            Kernel::Layout => "layout",
            Kernel::Strain => "strain",
            Kernel::Correlate => "correlate",
            Kernel::Copy => "copy",
        }
    }
}

/// What the kernels read, passed to each by value: `Args` in `cuda/kernels.cu`, member for
/// member (the comments are there). Addresses are the device's.
#[repr(C)]
#[derive(Clone, Copy, Debug, Default)]
pub struct Args {
    f: [u64; 9],
    m: [u64; 8],
    psi_v: [u64; 3],
    psi_s: [u64; 3],
    profile: u64,
    sponge: u64,
    term_field: u64,
    term_index: u64,
    term_weight: u64,
    term_source: u64,
    values: u64,
    gathered: u64,
    peak: u64,
    v: [u64; 3],
    e: [u64; 6],
    ea: [u64; 6],
    before: [u64; 3],
    after: [u64; 3],
    sums: u64,
    lit: u64,
    volume: f64,
    n: [i32; 3],
    halo: [i32; 3],
    padded: [i32; 3],
    lo: [i32; 3],
    hi: [i32; 3],
    profile_at: [i32; 3],
    l: i32,
    free_surface: i32,
    layers: i32,
    which: i32,
    from: i32,
    to: i32,
    row: i32,
    sources: i32,
    check: i32,
    skip: i32,
    pad: i32,
    sign: f32,
    c: [[f32; 8]; 3],
}

/// Somewhere the kernels run: a GPU ([`Cuda`]), or a CPU standing in for one ([`Host`]).
/// Memory is addressed as the device addresses it, and freed with the device.
pub trait Device {
    /// What it is, for people.
    fn name(&self) -> String;

    /// Makes the kernels ready for stencils `half_length` nodes long each way.
    fn prepare(&mut self, half_length: usize) -> Result<()>;

    /// Bytes it can still allocate, if it knows.
    fn free_bytes(&mut self) -> Result<Option<usize>>;

    /// Allocates `bytes` bytes, holding anything, and returns their address.
    fn alloc(&mut self, bytes: usize) -> Result<u64>;

    /// Copies `data` to the device at `to`.
    fn upload(&mut self, to: u64, data: &[u8]) -> Result<()>;

    /// Fills `out` from the device at `from`, once the kernels launched before have run.
    fn download(&mut self, from: u64, out: &mut [u8]) -> Result<()>;

    /// Starts `kernel` on `grid` blocks of `block` threads; it may still be running when this
    /// returns. Kernels run in the order they were launched.
    fn launch(
        &mut self,
        kernel: Kernel,
        grid: [u32; 3],
        block: [u32; 3],
        args: &Args,
    ) -> Result<()>;

    /// Returns once every kernel launched has run.
    fn wait(&mut self) -> Result<()>;
}

impl<D: Device + ?Sized> Device for Box<D> {
    fn name(&self) -> String {
        (**self).name()
    }

    fn prepare(&mut self, half_length: usize) -> Result<()> {
        (**self).prepare(half_length)
    }

    fn free_bytes(&mut self) -> Result<Option<usize>> {
        (**self).free_bytes()
    }

    fn alloc(&mut self, bytes: usize) -> Result<u64> {
        (**self).alloc(bytes)
    }

    fn upload(&mut self, to: u64, data: &[u8]) -> Result<()> {
        (**self).upload(to, data)
    }

    fn download(&mut self, from: u64, out: &mut [u8]) -> Result<()> {
        (**self).download(from, out)
    }

    fn launch(
        &mut self,
        kernel: Kernel,
        grid: [u32; 3],
        block: [u32; 3],
        args: &Args,
    ) -> Result<()> {
        (**self).launch(kernel, grid, block, args)
    }

    fn wait(&mut self) -> Result<()> {
        (**self).wait()
    }
}

/// One device that several hold, each launching on it in turn: what is launched runs in the
/// order launched, whoever launched it. The kernels are compiled once for each stencil
/// length.
pub struct Shared<D>(Rc<RefCell<(D, Option<usize>)>>);

impl<D> Shared<D> {
    pub fn new(device: D) -> Self {
        Self(Rc::new(RefCell::new((device, None))))
    }
}

impl<D> Clone for Shared<D> {
    fn clone(&self) -> Self {
        Self(Rc::clone(&self.0))
    }
}

impl<D: Device> Device for Shared<D> {
    fn name(&self) -> String {
        self.0.borrow().0.name()
    }

    fn prepare(&mut self, half_length: usize) -> Result<()> {
        let mut held = self.0.borrow_mut();
        if held.1 != Some(half_length) {
            held.0.prepare(half_length)?;
            held.1 = Some(half_length);
        }
        Ok(())
    }

    fn free_bytes(&mut self) -> Result<Option<usize>> {
        self.0.borrow_mut().0.free_bytes()
    }

    fn alloc(&mut self, bytes: usize) -> Result<u64> {
        self.0.borrow_mut().0.alloc(bytes)
    }

    fn upload(&mut self, to: u64, data: &[u8]) -> Result<()> {
        self.0.borrow_mut().0.upload(to, data)
    }

    fn download(&mut self, from: u64, out: &mut [u8]) -> Result<()> {
        self.0.borrow_mut().0.download(from, out)
    }

    fn launch(
        &mut self,
        kernel: Kernel,
        grid: [u32; 3],
        block: [u32; 3],
        args: &Args,
    ) -> Result<()> {
        self.0.borrow_mut().0.launch(kernel, grid, block, args)
    }

    fn wait(&mut self) -> Result<()> {
        self.0.borrow_mut().0.wait()
    }
}

/// The values of a slice as bytes, to cross to a device.
fn bytes<T: Plain>(v: &[T]) -> &[u8] {
    // SAFETY: `Plain` types have no padding and any bytes are a value of them.
    unsafe { std::slice::from_raw_parts(v.as_ptr().cast(), size_of_val(v)) }
}

fn bytes_mut<T: Plain>(v: &mut [T]) -> &mut [u8] {
    // SAFETY: as in `bytes`.
    unsafe { std::slice::from_raw_parts_mut(v.as_mut_ptr().cast(), size_of_val(v)) }
}

/// Numbers that are the same bytes on the device.
trait Plain: Copy {}
impl Plain for f32 {}
impl Plain for i32 {}
impl Plain for f64 {}
impl Plain for [f64; 8] {}

/// Threads per block of the short kernels that run along one index.
const LINE: u32 = 256;

/// Bytes of receivers' readings kept on the device before they are fetched.
const GATHERED_BYTES: usize = 64 << 20;

/// A solver's wavefield on a device, stepped there.
pub struct GpuSolver<D: Device> {
    device: D,
    args: Args,
    grid: Grid,
    dt: f64,
    free_surface: bool,
    sponge: bool,
    /// Threads along x, y and z of a block of the sweeps over the nodes.
    block: [u32; 3],
    /// Steps taken.
    pub step: usize,
    /// The terms of the forces, the moment tensors and the receivers, among all the terms on
    /// the device.
    forces: Range<usize>,
    moments: Range<usize>,
    receivers: Range<usize>,
    /// Steps the sources have values for.
    steps: usize,
    /// Samples the device keeps before they must be fetched, and how many it holds now.
    chunk: usize,
    held: usize,
    /// Bytes allocated on the device.
    bytes: usize,
}

impl<D: Device> GpuSolver<D> {
    /// The state of `solver` on `device`, ready to step from there.
    pub fn new(mut device: D, solver: &Solver) -> Result<Self> {
        let g = solver.grid;
        let l = solver.stencil.half_length();
        let sponge = solver.sponge();
        let memory = solver.layers().map(|p| p.each_ref().map(Vec::len));
        let layers: usize = memory.iter().flatten().sum();
        // The kernels index with 32-bit integers.
        ensure!(
            g.padded_len() < i32::MAX as usize && 3 * layers < i32::MAX as usize,
            "{} nodes are more than the GPU's kernels index",
            g.padded_len()
        );
        let need = 4 * (17 * g.padded_len() + layers);
        if let Some(free) = device.free_bytes()? {
            ensure!(
                need <= free,
                "the run needs {:.2} GiB on the GPU, which has {:.2} GiB free",
                need as f64 / (1u64 << 30) as f64,
                free as f64 / (1u64 << 30) as f64
            );
        }
        device.prepare(l)?;

        let mut gpu = Self {
            device,
            args: Args::default(),
            grid: g,
            dt: solver.dt,
            free_surface: solver.free_surface,
            sponge: sponge.is_some(),
            block: block_shape(&g)?,
            step: solver.step,
            forces: 0..0,
            moments: 0..0,
            receivers: 0..0,
            steps: 0,
            chunk: 0,
            held: 0,
            bytes: 0,
        };
        let int = |v: usize| v as i32;
        let a = &mut gpu.args;
        a.n = g.n.map(int);
        a.halo = g.halos().map(int);
        a.padded = g.padded().map(int);
        a.l = int(l);
        a.free_surface = solver.free_surface as i32;
        a.layers = 1;
        a.sign = 1.0;
        for axis in 0..3 {
            let inside = solver.cpml[axis].interior();
            a.lo[axis] = int(inside.start);
            a.hi[axis] = int(g.n[axis] - inside.end);
            // As Stencil::forward and backward scale them, in single precision.
            let ih = (1.0 / g.h[axis]) as f32;
            for (c, &s) in a.c[axis].iter_mut().zip(solver.stencil.single()) {
                *c = s * ih;
            }
        }

        // Each axis's layers: at the nodes a, b and 1/κ for every slab node, then the same
        // at the half nodes.
        let mut profile = Vec::new();
        for axis in 0..3 {
            gpu.args.profile_at[axis] = int(profile.len());
            for at in [Stagger::Node, Stagger::Half] {
                let p = solver.cpml[axis].profile(at);
                for part in [&p.a, &p.b, &p.inv_kappa] {
                    profile.extend_from_slice(part);
                }
            }
        }
        gpu.args.profile = gpu.put(&profile)?;
        if let Some(factors) = sponge {
            gpu.args.sponge = gpu.put(&factors.concat())?;
        }
        let m = &solver.material;
        for (to, values) in [&m.bx, &m.by, &m.bz, &m.lam, &m.l2m, &m.mxy, &m.mxz, &m.myz]
            .into_iter()
            .enumerate()
        {
            gpu.args.m[to] = gpu.put(values)?;
        }
        for q in 0..9 {
            gpu.args.f[q] = gpu.reserve(4 * g.padded_len())?;
        }
        let [of_v, of_s] = memory;
        for (axis, (v, s)) in of_v.into_iter().zip(of_s).enumerate() {
            gpu.args.psi_v[axis] = gpu.reserve(4 * v)?;
            gpu.args.psi_s[axis] = gpu.reserve(4 * s)?;
        }
        gpu.args.peak = gpu.reserve(4 * (g.n[1] * g.n[2]).max(4))?;
        gpu.check_layout()?;
        gpu.upload(solver)?;
        Ok(gpu)
    }

    /// What it runs on.
    pub fn device(&self) -> &D {
        &self.device
    }

    /// Bytes it holds on the device.
    pub fn memory_bytes(&self) -> usize {
        self.bytes
    }

    fn reserve(&mut self, bytes: usize) -> Result<u64> {
        self.bytes += bytes;
        self.device.alloc(bytes)
    }

    /// `values` on the device, in an array of their own.
    fn put<T: Plain>(&mut self, values: &[T]) -> Result<u64> {
        let at = self.reserve(size_of_val(values))?;
        self.device.upload(at, bytes(values))?;
        Ok(at)
    }

    /// Makes sure the kernels read the arguments as they are written here: a kernel sends
    /// back the structure's size and its last members.
    fn check_layout(&mut self) -> Result<()> {
        let mut args = self.args;
        (args.check, args.sign, args.c[2][7]) = (0x5eed, -3.5, 0.625);
        self.device.launch(Kernel::Layout, [1; 3], [1; 3], &args)?;
        let mut back = [0.0f32; 4];
        self.device.download(self.args.peak, bytes_mut(&mut back))?;
        ensure!(
            back == [size_of::<Args>() as f32, 0x5eed as f32, -3.5, 0.625],
            "the kernels read their arguments differently from how they are passed \
             ({back:?} came back)"
        );
        Ok(())
    }

    /// Takes `solver`'s fields, layers' memory and step count.
    pub fn upload(&mut self, solver: &Solver) -> Result<()> {
        ensure!(solver.grid == self.grid, "another grid's solver");
        self.upload_parts(solver.fields.arrays(), solver.layers(), solver.step)
    }

    /// Takes up `state`, one of this solver's ([`Solver::state`]).
    pub fn upload_state(&mut self, state: &State) -> Result<()> {
        self.upload_parts(state.fields().arrays(), state.layers(), state.step())
    }

    fn upload_parts(
        &mut self,
        fields: [&Vec<f32>; 9],
        [psi_v, psi_s]: [&[Vec<f32>; 3]; 2],
        step: usize,
    ) -> Result<()> {
        ensure!(
            fields.iter().all(|f| f.len() == self.grid.padded_len()),
            "another grid's fields"
        );
        for (at, values) in self.args.f.into_iter().zip(fields) {
            self.device.upload(at, bytes(values))?;
        }
        for (at, values) in
            (self.args.psi_v.into_iter().zip(psi_v)).chain(self.args.psi_s.into_iter().zip(psi_s))
        {
            self.device.upload(at, bytes(values))?;
        }
        self.step = step;
        Ok(())
    }

    /// Gives `solver` the fields, the layers' memory and the step count held here.
    pub fn download(&mut self, solver: &mut Solver) -> Result<()> {
        ensure!(solver.grid == self.grid, "another grid's solver");
        // Apart from the solver for a moment, to borrow both.
        let mut fields = std::mem::replace(&mut solver.fields, crate::solver::Fields::zeros(0));
        let done = Self::download_parts(
            &mut self.device,
            &self.args,
            fields.arrays_mut(),
            solver.layers_mut(),
        );
        solver.fields = fields;
        done?;
        solver.step = self.step;
        Ok(())
    }

    /// Gives `state` (one of this solver's) the fields, the layers' memory and the step count
    /// held here.
    pub fn download_state(&mut self, state: &mut State) -> Result<()> {
        let (fields, layers, step) = state.parts_mut();
        ensure!(
            fields
                .arrays()
                .iter()
                .all(|f| f.len() == self.grid.padded_len()),
            "another grid's state"
        );
        Self::download_parts(&mut self.device, &self.args, fields.arrays_mut(), layers)?;
        *step = self.step;
        Ok(())
    }

    fn download_parts(
        device: &mut D,
        args: &Args,
        fields: [&mut Vec<f32>; 9],
        [psi_v, psi_s]: [&mut [Vec<f32>; 3]; 2],
    ) -> Result<()> {
        for (at, values) in args.f.into_iter().zip(fields) {
            device.download(at, bytes_mut(values))?;
        }
        for (at, values) in (args.psi_v.into_iter().zip(psi_v.iter_mut()))
            .chain(args.psi_s.into_iter().zip(psi_s.iter_mut()))
        {
            device.download(at, bytes_mut(values))?;
        }
        Ok(())
    }

    /// Puts the sources of a run of `steps` steps, and its receivers, on the device:
    /// [`GpuSolver::step`] then adds the sources, and [`GpuSolver::record`] reads the
    /// receivers.
    pub fn attach(
        &mut self,
        sources: &[Source],
        steps: usize,
        receivers: &Receivers,
    ) -> Result<()> {
        // Forces, then moment tensors: each kind in the sources' order, as a step adds them.
        let of_kind = |forces: bool| -> Vec<(&PointOp, usize)> {
            (sources.iter().enumerate())
                .filter(|(_, s)| s.is_force() == forces)
                .map(|(i, s)| (&s.op, i))
                .collect()
        };
        // What each source adds at each step, as Solver::step_with asks for it: a force with
        // the velocities at n·dt, a moment tensor with the stresses half a step later.
        let mut values = Vec::with_capacity(steps * sources.len());
        for n in 0..steps {
            let t = n as f64 * self.dt;
            for source in sources {
                let t = if source.is_force() {
                    t
                } else {
                    t + 0.5 * self.dt
                };
                values.push(source.value(n, t));
            }
        }
        let receiving: Vec<&PointOp> = receivers.ops.iter().flatten().collect();
        self.attach_ops(
            [&of_kind(true), &of_kind(false)],
            sources.len(),
            &values,
            steps,
            &receiving,
            receivers.samples,
        )
    }

    /// Puts on the device what a run of `steps` steps adds and reads: `added` the operators
    /// added with the velocities, then those added with the stresses, each with its column in
    /// the row of `columns` values that `values` holds for each step; and the operators of the
    /// receivers, which read `samples` samples.
    fn attach_ops(
        &mut self,
        added: [&[(&PointOp, usize)]; 2],
        columns: usize,
        values: &[f32],
        steps: usize,
        receivers: &[&PointOp],
        samples: usize,
    ) -> Result<()> {
        ensure!(
            values.len() == steps * columns
                && added
                    .iter()
                    .flat_map(|a| a.iter())
                    .all(|&(_, c)| c < columns),
            "{} values for {steps} steps of {columns}",
            values.len()
        );
        let mut terms = Terms::default();
        for (k, ops) in added.into_iter().enumerate() {
            let start = terms.len();
            for &(op, column) in ops {
                terms.extend(op, column);
            }
            let range = start..terms.len();
            if k == 0 {
                self.forces = range;
            } else {
                self.moments = range;
            }
        }
        let start = terms.len();
        for op in receivers {
            terms.extend(op, 0);
        }
        self.receivers = start..terms.len();
        ensure!(
            terms.len() < i32::MAX as usize
                && steps < i32::MAX as usize
                && values.len() < i32::MAX as usize,
            "too many terms"
        );
        self.steps = steps;
        self.args.sources = columns as i32;
        self.args.values = self.put(values)?;
        self.args.term_field = self.put(&terms.field)?;
        self.args.term_index = self.put(&terms.index)?;
        self.args.term_weight = self.put(&terms.weight)?;
        self.args.term_source = self.put(&terms.source)?;

        let row = 4 * self.receivers.len();
        self.chunk = match row {
            0 => samples.max(1),
            _ => (GATHERED_BYTES / row).clamp(1, samples.max(1)),
        };
        self.held = 0;
        self.args.gathered = self.reserve(self.chunk * row)?;
        Ok(())
    }

    /// Advances one step, with the sources attached: [`Solver::step_with`], in its order.
    pub fn step(&mut self) -> Result<()> {
        (self.args.sign, self.args.layers) = (1.0, 1);
        self.sweep(Kernel::Velocity)?;
        self.inject(self.forces.clone())?;
        if self.free_surface {
            self.mirror(0, 1.0)?;
        }
        self.sweep(Kernel::Stress)?;
        self.inject(self.moments.clone())?;
        if self.sponge {
            self.sweep(Kernel::Damp)?;
        }
        if self.free_surface {
            self.mirror(1, -1.0)?;
        }
        self.step += 1;
        Ok(())
    }

    /// Launches a kernel that takes a thread per node.
    fn sweep(&mut self, kernel: Kernel) -> Result<()> {
        let grid = std::array::from_fn(|a| (self.grid.n[a] as u32).div_ceil(self.block[a]));
        self.device.launch(kernel, grid, self.block, &self.args)
    }

    /// Adds the sources' terms `terms` for the step being taken, if they have values for it.
    fn inject(&mut self, terms: Range<usize>) -> Result<()> {
        if terms.is_empty() || self.step >= self.steps {
            return Ok(());
        }
        let args = Args {
            from: terms.start as i32,
            to: terms.end as i32,
            row: self.step as i32,
            ..self.args
        };
        self.device.launch(Kernel::Inject, [1; 3], [1; 3], &args)
    }

    /// Mirrors the velocities (`which` 0) or the stresses that vanish on the free surface (1)
    /// across it, with `sign`.
    fn mirror(&mut self, which: i32, sign: f32) -> Result<()> {
        let args = Args {
            which,
            sign,
            ..self.args
        };
        let grid = [
            (self.grid.plane_len() as u32).div_ceil(LINE),
            self.grid.halos()[2] as u32,
            1,
        ];
        self.device
            .launch(Kernel::Mirror, grid, [LINE, 1, 1], &args)
    }

    /// Reads the receivers attached, as their next sample; `receivers` gets it, and those
    /// before it, when the device has kept as many as it has room for
    /// ([`GpuSolver::fetch`]). Nothing if they are full.
    pub fn record(&mut self, receivers: &mut Receivers) -> Result<()> {
        if receivers.recorded() + self.held == receivers.samples {
            return Ok(());
        }
        if !self.receivers.is_empty() {
            let args = Args {
                from: self.receivers.start as i32,
                to: self.receivers.end as i32,
                row: self.held as i32,
                ..self.args
            };
            let grid = [(self.receivers.len() as u32).div_ceil(LINE), 1, 1];
            self.device
                .launch(Kernel::Gather, grid, [LINE, 1, 1], &args)?;
        }
        self.held += 1;
        if self.held == self.chunk {
            self.fetch(receivers)?;
        }
        Ok(())
    }

    /// Hands `receivers` the samples the device still keeps: each trace's is the weighted sum
    /// of what its terms read, taken here exactly as [`PointOp::sample`] takes it.
    pub fn fetch(&mut self, receivers: &mut Receivers) -> Result<()> {
        let row = self.receivers.len();
        let mut gathered = vec![0.0f32; self.held * row];
        self.device
            .download(self.args.gathered, bytes_mut(&mut gathered))?;
        for m in 0..self.held {
            let mut read = &gathered[m * row..(m + 1) * row];
            receivers.record_with(|op| {
                let (values, rest) = read.split_at(op.terms.len());
                read = rest;
                op.weigh(values)
            });
        }
        self.held = 0;
        Ok(())
    }

    /// The largest particle-velocity component anywhere, or infinity if any is not finite
    /// ([`Solver::peak_velocity`]).
    pub fn peak_velocity(&mut self) -> Result<f32> {
        let [_, ny, nz] = self.grid.n;
        let grid = [(ny as u32).div_ceil(16), (nz as u32).div_ceil(16), 1];
        self.device
            .launch(Kernel::Peak, grid, [16, 16, 1], &self.args)?;
        let mut rows = vec![0.0f32; ny * nz];
        self.device.download(self.args.peak, bytes_mut(&mut rows))?;
        Ok(rows.into_iter().fold(0.0, f32::max))
    }

    /// Returns once the device has taken every step asked of it.
    pub fn wait(&mut self) -> Result<()> {
        self.device.wait()
    }
}

/// The terms of point operators, as the device holds them.
#[derive(Default)]
struct Terms {
    field: Vec<i32>,
    index: Vec<i32>,
    weight: Vec<f32>,
    /// The source whose value each takes.
    source: Vec<i32>,
}

impl Terms {
    fn len(&self) -> usize {
        self.field.len()
    }

    fn extend(&mut self, op: &PointOp, source: usize) {
        for &(field, index, weight) in &op.terms {
            let f = Field::ALL
                .iter()
                .position(|&f| f == field)
                .expect("a field");
            self.field.push(f as i32);
            self.index.push(index as i32);
            self.weight.push(weight);
            self.source.push(source as i32);
        }
    }
}

/// Threads along x, y and z of a block of the sweeps over the nodes of `g`: rows of 32 along
/// x, the direction memory runs, and a few rows and planes of them (one row on a 2D grid);
/// [`BLOCK_VARIABLE`] sets another.
fn block_shape(g: &Grid) -> Result<[u32; 3]> {
    if let Ok(text) = std::env::var(BLOCK_VARIABLE) {
        let parts: Vec<u32> = text
            .split(',')
            .map(|p| p.trim().parse())
            .collect::<Result<_, _>>()
            .with_context(|| format!("{BLOCK_VARIABLE}={text} isn't three numbers"))?;
        let block: [u32; 3] = parts
            .try_into()
            .ok()
            .with_context(|| format!("{BLOCK_VARIABLE}={text} isn't three numbers"))?;
        ensure!(
            block.iter().all(|&b| b > 0) && block.iter().product::<u32>() <= 1024,
            "{BLOCK_VARIABLE}={text}: a block is 1 to 1024 threads"
        );
        return Ok(block);
    }
    Ok(if g.is_2d() { [64, 1, 4] } else { [32, 4, 4] })
}

/// Runs every step of `sim` on `device` as [`Simulation::run`] does on the CPU, and leaves
/// it as that would: the traces, the fields and the snapshots written into `out` are the
/// same, bit for bit unless the device says otherwise ([`Device::name`]).
pub fn run<D: Device>(
    sim: &mut Simulation,
    device: D,
    out: Option<&Path>,
    progress: bool,
) -> Result<()> {
    let snap_dir = out.map(|d| d.join("snapshots"));
    if let Some(dir) = &snap_dir
        && !sim.snapshots.is_empty()
    {
        std::fs::create_dir_all(dir).with_context(|| format!("creating {}", dir.display()))?;
    }
    let mut gpu = GpuSolver::new(device, &sim.solver)?;
    gpu.attach(&sim.sources, sim.steps, &sim.receivers)?;
    let name = gpu.device().name();
    if progress {
        eprintln!("device    {name}");
    }
    sim.info.device = Some(name);
    let start = Instant::now();
    let mut last_report = start;
    let cells = sim.solver.grid.interior_len() as f64;
    let first = gpu.step;
    while gpu.step < sim.steps {
        let n = gpu.step;
        gpu.step()?;
        if n.is_multiple_of(sim.record_every) {
            gpu.record(&mut sim.receivers)?;
        }
        let done = gpu.step;
        if let Some(dir) = &snap_dir
            && sim.snapshots.iter().any(|s| s.due(done))
        {
            gpu.download(&mut sim.solver)?;
            for s in sim.snapshots.iter().filter(|s| s.due(done)) {
                s.write(dir, &sim.solver)?;
            }
        }
        // The steps are queued, not waited for, so the clock is read only now and then; the
        // check below waits for them.
        let report = progress
            && (done == sim.steps
                || (done % 20 == 0 && last_report.elapsed() >= Duration::from_secs(2)));
        if report || done % 100 == 0 || done == sim.steps {
            let peak = gpu.peak_velocity()?;
            if !peak.is_finite() {
                bail!(
                    "the wavefield blew up by step {done}; check the model, or lower dt or courant"
                );
            }
            if report {
                last_report = Instant::now();
                let elapsed = start.elapsed().as_secs_f64();
                let eta = elapsed / (done - first) as f64 * (sim.steps - done) as f64;
                eprintln!(
                    "step {done:>6}/{}  t {:.4} s  peak |v| {peak:.3e} m/s  {:.0} Mcells/s  {} left",
                    sim.steps,
                    done as f64 * sim.solver.dt,
                    cells * (done - first) as f64 / elapsed / 1e6,
                    human_time(eta),
                );
            }
        }
    }
    gpu.fetch(&mut sim.receivers)?;
    gpu.download(&mut sim.solver)?;
    sim.finish();
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::config::Config;

    /// A solid under water with a free surface, uneven spacing, a force and a moment tensor
    /// between nodes, and receivers of every component on and off the nodes: every kernel and
    /// every branch of them has work.
    const SURFACE: &str = r#"
        [grid]
        points = [26, 22, 24]
        spacing = [10.0, 12.0, 8.0]
        [time]
        duration = 0.12
        [boundary]
        pml_width = 6
        [[model.layer]]
        top = 0.0
        vp = 1500.0
        vs = 0.0
        rho = 1000.0
        [[model.layer]]
        top = 40.0
        vp = 2400.0
        vs = 1200.0
        rho = 2100.0
        [[source]]
        position = [120.0, 130.0, 64.0]
        double_couple = { strike = 30.0, dip = 70.0, rake = 20.0, m0 = 1e12 }
        wavelet = { kind = "ricker", frequency = 12.0 }
        [[source]]
        position = [83.0, 121.0, 100.0]
        force = [1e9, -2e9, 3e9]
        wavelet = { kind = "ricker", frequency = 10.0 }
        [receivers]
        components = ["vx", "vy", "vz", "p"]
        [[receivers.point]]
        position = [100.0, 120.0, 0.0]
        [[receivers.point]]
        position = [147.0, 101.0, 21.0]
        [[receivers.line]]
        from = [80.0, 96.0, 48.0]
        to = [170.0, 168.0, 120.0]
        count = 7
    "#;

    /// Absorbing layers on every side, fourth-order Taylor stencils, the time error kept.
    const FULL_SPACE: &str = r#"
        [grid]
        points = [24, 26, 22]
        spacing = 10.0
        [time]
        duration = 0.1
        [boundary]
        pml_width = 5
        free_surface = false
        [scheme]
        order = 4
        coefficients = "taylor"
        remove_time_dispersion = false
        [[model.layer]]
        top = 0.0
        vp = 3000.0
        vs = 1700.0
        rho = 2500.0
        [[source]]
        position = [110.0, 120.0, 100.0]
        explosion = 1e12
        wavelet = { kind = "ricker", frequency = 15.0 }
        [receivers]
        components = ["vz", "p"]
        [[receivers.line]]
        from = [60.0, 60.0, 60.0]
        to = [170.0, 180.0, 150.0]
        count = 5
    "#;

    /// A 2D grid, sixteenth order.
    const PLANE: &str = r#"
        [grid]
        points = [60, 1, 50]
        spacing = 10.0
        [time]
        duration = 0.15
        [boundary]
        pml_width = 10
        [scheme]
        order = 16
        [[model.layer]]
        top = 0.0
        vp = 2000.0
        vs = 900.0
        rho = 2000.0
        [[model.layer]]
        top = 200.0
        vp = 3200.0
        vs = 1800.0
        rho = 2400.0
        [[source]]
        position = [300.0, 0.0, 0.0]
        force = [0.0, 0.0, 1e9]
        wavelet = { kind = "ricker", frequency = 8.0 }
        [receivers]
        components = ["vx", "vz", "p"]
        [[receivers.line]]
        from = [150.0, 0.0, 0.0]
        to = [450.0, 0.0, 0.0]
        count = 11
    "#;

    fn simulation(text: &str) -> Simulation {
        Simulation::new(&Config::from_toml(text).unwrap(), Path::new(".")).unwrap()
    }

    /// Every bit of what a run leaves: the fields, the layers' memory and the traces.
    fn bits(sim: &Simulation) -> Vec<(String, Vec<u32>)> {
        let [psi_v, psi_s] = sim.solver.layers();
        let fields = Field::ALL
            .iter()
            .zip(sim.solver.fields.arrays())
            .map(|(f, v)| (format!("{f:?}"), v));
        let memory = psi_v
            .iter()
            .chain(psi_s)
            .enumerate()
            .map(|(m, v)| (format!("layer memory {m}"), v));
        let traces = sim
            .receivers
            .components
            .iter()
            .zip(&sim.receivers.traces)
            .map(|(c, v)| (format!("traces of {}", c.name()), v));
        fields
            .chain(memory)
            .chain(traces)
            .map(|(name, v)| (name, v.iter().map(|x| x.to_bits()).collect()))
            .collect()
    }

    /// Runs `text` on the CPU and on `device`, and wants the same bits from both.
    fn runs_alike<D: Device>(name: &str, text: &str, device: D) {
        let mut cpu = simulation(text);
        cpu.run(None, false).unwrap();
        let mut gpu = simulation(text);
        let on = device.name();
        run(&mut gpu, device, None, false).unwrap();
        assert_eq!(gpu.solver.step, cpu.solver.step);
        let peak = cpu.solver.peak_velocity();
        assert!(peak > 0.0 && peak.is_finite(), "{name}: nothing to compare");
        for ((what, a), (_, b)) in bits(&cpu).into_iter().zip(bits(&gpu)) {
            let differ = a.iter().zip(&b).filter(|(x, y)| x != y).count();
            assert!(
                a.len() == b.len() && differ == 0,
                "{name}, {what}: {differ} of {} values differ on {on}",
                a.len()
            );
        }
    }

    const CASES: [(&str, &str); 3] = [
        ("water over rock under a free surface", SURFACE),
        ("a full space", FULL_SPACE),
        ("a 2D grid", PLANE),
    ];

    /// The kernels' source, compiled for this CPU and run a thread at a time, leaves the
    /// fields, the layers' memory and the traces the solver leaves, bit for bit.
    #[test]
    fn the_kernels_step_as_the_solver_bit_for_bit() {
        for (name, text) in CASES {
            runs_alike(name, text, Host::build().unwrap());
        }
    }

    /// The environment variable that makes the GPU's tests fail, instead of pass idly, where
    /// there is no GPU: set it where there should be one.
    const REQUIRE: &str = "VIBE_RATIONS_CUDA_REQUIRED";

    /// A GPU to test on, or None where there is none (unless [`REQUIRE`] is set).
    fn gpu() -> Option<Cuda> {
        match Cuda::open() {
            Ok(device) => Some(device.fuse(false)),
            Err(e) if std::env::var_os(REQUIRE).is_some() => panic!("no GPU: {e:#}"),
            Err(e) => {
                eprintln!("no GPU to test on, skipped: {e:#}");
                None
            }
        }
    }

    /// A GPU, where there is one, leaves the same bits as the CPU too.
    #[test]
    fn a_gpu_steps_as_the_solver_bit_for_bit() {
        for (name, text) in CASES {
            let Some(device) = gpu() else { return };
            eprintln!("{name} on {}", device.name());
            runs_alike(name, text, device);
        }
    }

    /// Stepped by hand: receivers fetched in chunks smaller than the run, the fields taken
    /// back and handed over again midway, and the largest velocity, all as the solver's.
    #[test]
    fn fetches_in_chunks_and_resumes_from_a_download() {
        let mut cpu = simulation(SURFACE);
        let mut sim = simulation(SURFACE);
        let mut gpu = GpuSolver::new(Host::build().unwrap(), &sim.solver).unwrap();
        gpu.attach(&sim.sources, sim.steps, &sim.receivers).unwrap();
        assert!(gpu.chunk > 3 && sim.steps > 20);
        gpu.chunk = 3;
        for n in 0..sim.steps {
            cpu.advance();
            gpu.step().unwrap();
            if n.is_multiple_of(sim.record_every) {
                gpu.record(&mut sim.receivers).unwrap();
            }
            if n == 17 {
                assert_eq!(gpu.peak_velocity().unwrap(), cpu.solver.peak_velocity());
                gpu.download(&mut sim.solver).unwrap();
                assert_eq!(sim.solver.step, 18);
                gpu.upload(&sim.solver).unwrap();
            }
        }
        gpu.fetch(&mut sim.receivers).unwrap();
        gpu.download(&mut sim.solver).unwrap();
        assert_eq!(sim.receivers.recorded(), cpu.receivers.recorded());
        assert_eq!(bits(&sim), bits(&cpu));
        assert!(gpu.memory_bytes() > 17 * 4 * sim.solver.grid.padded_len());
    }

    /// The L2 misfit of the traces against no data: residuals the traces themselves, so
    /// that every receiver drives the adjoint run.
    fn no_data(r: &Receivers) -> Result<(f64, Vec<Vec<f32>>)> {
        let misfit = r.traces.iter().flatten().map(|&x| 0.5 * (x as f64).powi(2));
        Ok((misfit.sum(), r.traces.clone()))
    }

    /// A gradient of `text` on `device`, with another checkpoint interval than the CPU's,
    /// leaves the CPU's bits: the gradient, the illumination, the misfit, the traces and the
    /// forward run's last state.
    fn gradients_alike<D: Device>(name: &str, text: &str, device: D) {
        let mut cpu = simulation(text);
        let a = crate::gradient::gradient(&mut cpu, no_data, None).unwrap();
        let mut gpu = simulation(text);
        let on = device.name();
        let b = gradient(&mut gpu, no_data, Some(7), device).unwrap();
        assert_ne!(a.checkpoint_every, b.checkpoint_every);
        assert!(
            a.gradient.mu.iter().any(|&g| g != 0.0),
            "{name}: nothing to compare"
        );
        assert_eq!(a.misfit.to_bits(), b.misfit.to_bits(), "{name}: the misfit");
        let parts = [
            ("ρ", &a.gradient.rho, &b.gradient.rho),
            ("λ", &a.gradient.lambda, &b.gradient.lambda),
            ("μ", &a.gradient.mu, &b.gradient.mu),
            ("the illumination", &a.illumination, &b.illumination),
        ];
        for (what, x, y) in parts {
            let differ = x.iter().zip(y).filter(|(p, q)| p.to_bits() != q.to_bits());
            assert!(
                x.len() == y.len() && differ.count() == 0,
                "{name}, {what}: values differ on {on}"
            );
        }
        let traces = |s: &Simulation| bits(s).into_iter().filter(|(w, _)| w.starts_with("traces"));
        assert!(traces(&cpu).eq(traces(&gpu)), "{name}: the traces");
    }

    /// Water over rock under a free surface with pressure receivers, a full space at order
    /// 4, and a 2D grid at order 16.
    #[test]
    fn a_gradient_on_the_kernels_is_the_solvers_bit_for_bit() {
        for (name, text) in CASES {
            gradients_alike(name, text, Host::build().unwrap());
        }
    }

    /// A GPU, where there is one, makes the CPU's gradients too.
    #[test]
    fn a_gpu_makes_the_solvers_gradient_bit_for_bit() {
        for (name, text) in CASES {
            let Some(device) = gpu() else { return };
            eprintln!("{name} on {}", device.name());
            gradients_alike(name, text, device);
        }
    }

    /// NVRTC, where there is one, compiles the kernels without fusing any product into a sum.
    #[test]
    fn nvrtc_compiles_the_kernels_without_fused_multiply_adds() {
        let options = |fused: bool| {
            [
                "--gpu-architecture=compute_61".to_string(),
                format!("--fmad={fused}"),
                "-DSTENCIL_L=4".to_string(),
            ]
        };
        let (ptx, version) = match compile(KERNELS, &options(false)) {
            Ok(done) => done,
            Err(e) if std::env::var_os(REQUIRE).is_some() => panic!("no NVRTC: {e:#}"),
            Err(e) => return eprintln!("no NVRTC to test with, skipped: {e:#}"),
        };
        let ptx = String::from_utf8_lossy(&ptx).into_owned();
        eprintln!("NVRTC {version}: {} bytes of PTX", ptx.len());
        for kernel in Kernel::ALL {
            assert!(
                ptx.contains(&format!(".entry {}(", kernel.name())),
                "{kernel:?}"
            );
        }
        // The structure arrives whole, by value.
        assert!(ptx.contains(&format!("velocity_param_0[{}]", size_of::<Args>())));
        assert!(!ptx.contains("fma."));
        let (fused, _) = compile(KERNELS, &options(true)).unwrap();
        assert!(String::from_utf8_lossy(&fused).contains("fma."));
        // An error comes with the compiler's words.
        let e = compile("__global__ void broken( {", &options(false)).unwrap_err();
        assert!(format!("{e:#}").contains("error"), "{e:#}");
    }

    #[test]
    fn the_arguments_are_laid_out_as_the_kernels_declare_them() {
        // 55 addresses, a double and 54 four-byte members, with no padding: cuda/kernels.cu.
        assert_eq!(size_of::<Args>(), 56 * 8 + 54 * 4);
        assert_eq!(align_of::<Args>(), 8);
        assert_eq!(std::mem::offset_of!(Args, n), 448);
        assert_eq!(std::mem::offset_of!(Args, c), size_of::<Args>() - 96);
        for (k, kernel) in Kernel::ALL.into_iter().enumerate() {
            assert_eq!(kernel as usize, k);
            assert!(KERNELS.contains(&format!("__global__ void {}(Args a)", kernel.name())));
        }
    }
}
