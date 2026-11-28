//! An NVIDIA GPU through the CUDA driver API (`libcuda.so.1`, which comes with the driver)
//! and the kernels compiled for it by NVRTC (`libnvrtc.so`, CUDA's run-time compiler), both
//! found when a run asks for them: nothing is linked when this is built, and no crate stands
//! between. The handful of calls used are declared here from CUDA 12.4's `cuda.h` and
//! `nvrtc.h`; the symbols are the ones those headers' macros name (`cuCtxCreate_v2` and so
//! on).

use std::ffi::{CStr, CString, c_char, c_int, c_uint, c_void};
use std::ptr::null_mut;

use anyhow::{Context, Result, anyhow, bail, ensure};

use super::{Args, Device, KERNELS, Kernel};

unsafe extern "C" {
    fn dlopen(filename: *const c_char, flags: c_int) -> *mut c_void;
    fn dlsym(handle: *mut c_void, symbol: *const c_char) -> *mut c_void;
}

/// Resolve every symbol when the library is opened (the same value on Linux and macOS).
const RTLD_NOW: c_int = 2;

/// The environment variable that names NVRTC's library, where the loader wouldn't find it.
pub const NVRTC_VARIABLE: &str = "VIBE_RATIONS_NVRTC";

/// A shared library opened at run time. It stays open for the life of the process: closing
/// the driver's library under it is not worth the risk.
pub(super) struct Library {
    handle: *mut c_void,
    name: String,
}

impl Library {
    /// The first of `names` the system's loader can open.
    pub fn open<S: AsRef<str>>(names: &[S]) -> Result<Self> {
        for name in names {
            let c = CString::new(name.as_ref())?;
            // SAFETY: a NUL-terminated name; the loader runs the library's initialisers, as
            // for any library a program links.
            let handle = unsafe { dlopen(c.as_ptr(), RTLD_NOW) };
            if !handle.is_null() {
                return Ok(Self {
                    handle,
                    name: name.as_ref().to_string(),
                });
            }
        }
        let names: Vec<&str> = names.iter().map(AsRef::as_ref).collect();
        bail!("the loader finds none of {}", names.join(", "))
    }

    /// The function called `name`.
    ///
    /// # Safety
    ///
    /// `T` must be the `extern "C"` function pointer type of that function.
    pub unsafe fn function<T: Copy>(&self, name: &str) -> Result<T> {
        assert_eq!(size_of::<T>(), size_of::<*mut c_void>());
        let c = CString::new(name)?;
        // SAFETY: an open library and a NUL-terminated name.
        let p = unsafe { dlsym(self.handle, c.as_ptr()) };
        ensure!(!p.is_null(), "{} has no {name}", self.name);
        // SAFETY: a function pointer is an address, and the caller vouches for its type.
        Ok(unsafe { std::mem::transmute_copy::<*mut c_void, T>(&p) })
    }
}

/// The driver's calls: `CUresult` is an int, `CUdevice` an int, `CUdeviceptr` 64 bits, and
/// contexts, modules, functions and streams are pointers.
#[derive(Clone, Copy)]
struct Api {
    init: unsafe extern "C" fn(c_uint) -> c_int,
    error_name: unsafe extern "C" fn(c_int, *mut *const c_char) -> c_int,
    error_string: unsafe extern "C" fn(c_int, *mut *const c_char) -> c_int,
    device_get: unsafe extern "C" fn(*mut c_int, c_int) -> c_int,
    device_name: unsafe extern "C" fn(*mut c_char, c_int, c_int) -> c_int,
    device_attribute: unsafe extern "C" fn(*mut c_int, c_int, c_int) -> c_int,
    context_create: unsafe extern "C" fn(*mut *mut c_void, c_uint, c_int) -> c_int,
    context_destroy: unsafe extern "C" fn(*mut c_void) -> c_int,
    synchronize: unsafe extern "C" fn() -> c_int,
    memory_info: unsafe extern "C" fn(*mut usize, *mut usize) -> c_int,
    module_load: unsafe extern "C" fn(*mut *mut c_void, *const c_void) -> c_int,
    module_function: unsafe extern "C" fn(*mut *mut c_void, *mut c_void, *const c_char) -> c_int,
    alloc: unsafe extern "C" fn(*mut u64, usize) -> c_int,
    to_device: unsafe extern "C" fn(u64, *const c_void, usize) -> c_int,
    to_host: unsafe extern "C" fn(*mut c_void, u64, usize) -> c_int,
    #[allow(clippy::type_complexity)]
    launch: unsafe extern "C" fn(
        *mut c_void,
        c_uint,
        c_uint,
        c_uint,
        c_uint,
        c_uint,
        c_uint,
        c_uint,
        *mut c_void,
        *mut *mut c_void,
        *mut *mut c_void,
    ) -> c_int,
}

impl Api {
    fn load(l: &Library) -> Result<Self> {
        // SAFETY (all): each type is the function's in cuda.h.
        unsafe {
            Ok(Self {
                init: l.function("cuInit")?,
                error_name: l.function("cuGetErrorName")?,
                error_string: l.function("cuGetErrorString")?,
                device_get: l.function("cuDeviceGet")?,
                device_name: l.function("cuDeviceGetName")?,
                device_attribute: l.function("cuDeviceGetAttribute")?,
                context_create: l.function("cuCtxCreate_v2")?,
                context_destroy: l.function("cuCtxDestroy_v2")?,
                synchronize: l.function("cuCtxSynchronize")?,
                memory_info: l.function("cuMemGetInfo_v2")?,
                module_load: l.function("cuModuleLoadData")?,
                module_function: l.function("cuModuleGetFunction")?,
                alloc: l.function("cuMemAlloc_v2")?,
                to_device: l.function("cuMemcpyHtoD_v2")?,
                to_host: l.function("cuMemcpyDtoH_v2")?,
                launch: l.function("cuLaunchKernel")?,
            })
        }
    }

    /// An error for a call that returned `code`, unless it is CUDA_SUCCESS.
    fn check(&self, code: c_int, what: &str) -> Result<()> {
        if code == 0 {
            return Ok(());
        }
        let text = |f: unsafe extern "C" fn(c_int, *mut *const c_char) -> c_int| {
            let mut p: *const c_char = std::ptr::null();
            // SAFETY: the driver sets `p` to a static string, or leaves it null.
            unsafe {
                f(code, &mut p);
                (!p.is_null()).then(|| CStr::from_ptr(p).to_string_lossy().into_owned())
            }
        };
        let name = text(self.error_name).unwrap_or_else(|| format!("error {code}"));
        let why = text(self.error_string).unwrap_or_default();
        Err(anyhow!("{what}: {name} ({why})"))
    }
}

/// CU_DEVICE_ATTRIBUTE_COMPUTE_CAPABILITY_MAJOR and _MINOR.
const CAPABILITY: [c_int; 2] = [75, 76];

/// CU_CTX_SCHED_BLOCKING_SYNC: this thread sleeps while it waits for the GPU instead of
/// spinning, so a run takes no core from whatever else the machine is doing.
const BLOCKING_SYNC: c_uint = 4;

/// The first NVIDIA GPU of this machine, with a context of its own.
pub struct Cuda {
    api: Api,
    context: *mut c_void,
    functions: Option<[*mut c_void; Kernel::ALL.len()]>,
    name: String,
    capability: [c_int; 2],
    fused: bool,
    compiler: String,
}

impl Cuda {
    /// Opens the driver and makes a context on the first GPU; an error says what is missing.
    pub fn open() -> Result<Self> {
        let library = Library::open(&["libcuda.so.1", "libcuda.so"])
            .context("no NVIDIA driver here (libcuda.so.1 comes with it)")?;
        let api = Api::load(&library)?;
        let mut device = 0;
        let mut name = [0 as c_char; 128];
        let mut capability = [0; 2];
        let mut context = null_mut();
        // SAFETY: the driver's calls as cuda.h declares them, with room for what they write.
        unsafe {
            api.check((api.init)(0), "cuInit")?;
            api.check((api.device_get)(&mut device, 0), "cuDeviceGet")
                .context("no CUDA device here")?;
            api.check(
                (api.device_name)(name.as_mut_ptr(), name.len() as c_int, device),
                "cuDeviceGetName",
            )?;
            for (value, attribute) in capability.iter_mut().zip(CAPABILITY) {
                api.check(
                    (api.device_attribute)(value, attribute, device),
                    "cuDeviceGetAttribute",
                )?;
            }
            api.check(
                (api.context_create)(&mut context, BLOCKING_SYNC, device),
                "cuCtxCreate",
            )?;
        }
        // SAFETY: the driver wrote a NUL-terminated name into the buffer.
        let name = unsafe { CStr::from_ptr(name.as_ptr()) }
            .to_string_lossy()
            .into_owned();
        Ok(Self {
            api,
            context,
            functions: None,
            name,
            capability,
            fused: false,
            compiler: String::new(),
        })
    }

    /// Whether the kernels may fuse products into sums ([`Cuda::fuse`]).
    pub fn fused(&self) -> bool {
        self.fused
    }

    /// Lets the compiler fuse the kernels' products into sums, each rounded once instead of
    /// twice (`run --gpu --fma`), or not, as a device is opened. Fused, the results are no
    /// longer the CPU's bit for bit. Before the kernels are compiled.
    pub fn fuse(mut self, on: bool) -> Self {
        self.fused = on;
        self
    }
}

impl Drop for Cuda {
    fn drop(&mut self) {
        // SAFETY: the context made in `open`; destroying it frees what was allocated in it.
        unsafe {
            (self.api.context_destroy)(self.context);
        }
    }
}

impl Device for Cuda {
    fn name(&self) -> String {
        let [major, minor] = self.capability;
        let arithmetic = if self.fused {
            "products fused into sums: not the CPU's results bit for bit"
        } else {
            "the CPU's arithmetic"
        };
        let compiler = if self.compiler.is_empty() {
            String::new()
        } else {
            format!(", kernels compiled by NVRTC {}", self.compiler)
        };
        format!(
            "{} (CUDA, compute capability {major}.{minor}{compiler}; {arithmetic})",
            self.name
        )
    }

    fn prepare(&mut self, half_length: usize) -> Result<()> {
        let [major, minor] = self.capability;
        let options = [
            format!("--gpu-architecture=compute_{major}{minor}"),
            format!("--fmad={}", self.fused),
            format!("-DSTENCIL_L={half_length}"),
        ];
        let (ptx, version) = compile(KERNELS, &options)?;
        self.compiler = version;
        let mut module = null_mut();
        let mut functions = [null_mut(); Kernel::ALL.len()];
        // SAFETY: NUL-terminated PTX and names; the module lives as long as the context.
        unsafe {
            self.api
                .check(
                    (self.api.module_load)(&mut module, ptx.as_ptr().cast()),
                    "cuModuleLoadData",
                )
                .context("the driver did not take the kernels NVRTC compiled")?;
            for (f, kernel) in functions.iter_mut().zip(Kernel::ALL) {
                let name = CString::new(kernel.name())?;
                self.api.check(
                    (self.api.module_function)(f, module, name.as_ptr()),
                    "cuModuleGetFunction",
                )?;
            }
        }
        self.functions = Some(functions);
        Ok(())
    }

    fn free_bytes(&mut self) -> Result<Option<usize>> {
        let (mut free, mut total) = (0, 0);
        // SAFETY: two places for the driver to write a size.
        let code = unsafe { (self.api.memory_info)(&mut free, &mut total) };
        self.api.check(code, "cuMemGetInfo")?;
        Ok(Some(free))
    }

    fn alloc(&mut self, bytes: usize) -> Result<u64> {
        let mut p = 0;
        // SAFETY: a place for the driver to write the address.
        let code = unsafe { (self.api.alloc)(&mut p, bytes.max(1)) };
        self.api
            .check(code, "cuMemAlloc")
            .with_context(|| format!("allocating {bytes} bytes on the GPU"))?;
        Ok(p)
    }

    fn upload(&mut self, to: u64, data: &[u8]) -> Result<()> {
        if data.is_empty() {
            return Ok(());
        }
        // SAFETY: `data` is readable for its length; the caller allocated as much at `to`.
        let code = unsafe { (self.api.to_device)(to, data.as_ptr().cast(), data.len()) };
        self.api.check(code, "cuMemcpyHtoD")
    }

    fn download(&mut self, from: u64, out: &mut [u8]) -> Result<()> {
        if out.is_empty() {
            return Ok(());
        }
        // SAFETY: `out` is writable for its length; the caller allocated as much at `from`.
        // The copy waits for the kernels launched before it.
        let code = unsafe { (self.api.to_host)(out.as_mut_ptr().cast(), from, out.len()) };
        self.api.check(code, "cuMemcpyDtoH")
    }

    fn launch(
        &mut self,
        kernel: Kernel,
        grid: [u32; 3],
        block: [u32; 3],
        args: &Args,
    ) -> Result<()> {
        let functions = self.functions.context("the kernels are not compiled yet")?;
        let mut args = *args;
        // The kernels take one parameter, the structure by value: the driver copies it from
        // here before the call returns.
        let mut parameters = [(&raw mut args).cast::<c_void>()];
        // SAFETY: a function of the module, with the one parameter it declares.
        let code = unsafe {
            (self.api.launch)(
                functions[kernel as usize],
                grid[0],
                grid[1],
                grid[2],
                block[0],
                block[1],
                block[2],
                0,
                null_mut(),
                parameters.as_mut_ptr(),
                null_mut(),
            )
        };
        // CUDA_ERROR_LAUNCH_OUT_OF_RESOURCES when a block's threads need more registers
        // between them than the card gives a block.
        self.api.check(code, kernel.name()).with_context(|| {
            let [x, y, z] = block;
            format!(
                "launching blocks of {x} × {y} × {z} threads ({} sets fewer)",
                super::BLOCK_VARIABLE
            )
        })
    }

    fn wait(&mut self) -> Result<()> {
        // SAFETY: no arguments; it blocks until the GPU has run what was launched.
        let code = unsafe { (self.api.synchronize)() };
        self.api.check(code, "cuCtxSynchronize")
    }
}

/// NVRTC's calls: `nvrtcResult` is an int and a program a pointer.
struct Compiler {
    version: unsafe extern "C" fn(*mut c_int, *mut c_int) -> c_int,
    error_string: unsafe extern "C" fn(c_int) -> *const c_char,
    #[allow(clippy::type_complexity)]
    create: unsafe extern "C" fn(
        *mut *mut c_void,
        *const c_char,
        *const c_char,
        c_int,
        *const *const c_char,
        *const *const c_char,
    ) -> c_int,
    destroy: unsafe extern "C" fn(*mut *mut c_void) -> c_int,
    compile: unsafe extern "C" fn(*mut c_void, c_int, *const *const c_char) -> c_int,
    ptx_size: unsafe extern "C" fn(*mut c_void, *mut usize) -> c_int,
    ptx: unsafe extern "C" fn(*mut c_void, *mut c_char) -> c_int,
    log_size: unsafe extern "C" fn(*mut c_void, *mut usize) -> c_int,
    log: unsafe extern "C" fn(*mut c_void, *mut c_char) -> c_int,
}

impl Compiler {
    fn check(&self, code: c_int, what: &str) -> Result<()> {
        if code == 0 {
            return Ok(());
        }
        // SAFETY: NVRTC returns a static string for any code.
        let why = unsafe { CStr::from_ptr((self.error_string)(code)) };
        Err(anyhow!("{what}: {}", why.to_string_lossy()))
    }
}

/// Compiles CUDA C `source` with NVRTC to PTX, the GPU's assembly, as text ending in a NUL,
/// with NVRTC's version: what the driver then compiles for its GPU. It needs no GPU.
pub fn compile(source: &str, options: &[String]) -> Result<(Vec<u8>, String)> {
    let names = match std::env::var(NVRTC_VARIABLE) {
        Ok(path) if !path.trim().is_empty() => vec![path],
        _ => [
            "libnvrtc.so.12",
            "libnvrtc.so.13",
            "libnvrtc.so.11.2",
            "libnvrtc.so",
        ]
        .map(String::from)
        .to_vec(),
    };
    let l = Library::open(&names).with_context(|| {
        format!(
            "no NVRTC here to compile the kernels (CUDA's libnvrtc; {NVRTC_VARIABLE} can name it)"
        )
    })?;
    // SAFETY (all): each type is the function's in nvrtc.h.
    let c = unsafe {
        Compiler {
            version: l.function("nvrtcVersion")?,
            error_string: l.function("nvrtcGetErrorString")?,
            create: l.function("nvrtcCreateProgram")?,
            destroy: l.function("nvrtcDestroyProgram")?,
            compile: l.function("nvrtcCompileProgram")?,
            ptx_size: l.function("nvrtcGetPTXSize")?,
            ptx: l.function("nvrtcGetPTX")?,
            log_size: l.function("nvrtcGetProgramLogSize")?,
            log: l.function("nvrtcGetProgramLog")?,
        }
    };
    let source = CString::new(source)?;
    let name = CString::new("kernels.cu")?;
    let options: Vec<CString> = options
        .iter()
        .map(|o| CString::new(o.as_str()))
        .collect::<Result<_, _>>()?;
    let pointers: Vec<*const c_char> = options.iter().map(|o| o.as_ptr()).collect();
    let (mut major, mut minor) = (0, 0);
    let mut program = null_mut();
    // SAFETY: NVRTC's calls as nvrtc.h declares them, on NUL-terminated strings that outlive
    // them, with buffers of the sizes it asked for.
    unsafe {
        c.check((c.version)(&mut major, &mut minor), "nvrtcVersion")?;
        c.check(
            (c.create)(
                &mut program,
                source.as_ptr(),
                name.as_ptr(),
                0,
                std::ptr::null(),
                std::ptr::null(),
            ),
            "nvrtcCreateProgram",
        )?;
        let compiled = (c.compile)(program, pointers.len() as c_int, pointers.as_ptr());
        let text = |size: unsafe extern "C" fn(*mut c_void, *mut usize) -> c_int,
                    get: unsafe extern "C" fn(*mut c_void, *mut c_char) -> c_int,
                    what: &str|
         -> Result<Vec<u8>> {
            let mut n = 0;
            c.check(size(program, &mut n), what)?;
            let mut out = vec![0u8; n.max(1)];
            c.check(get(program, out.as_mut_ptr().cast()), what)?;
            Ok(out)
        };
        let result = match c.check(compiled, "compiling the kernels") {
            Ok(()) => text(c.ptx_size, c.ptx, "nvrtcGetPTX"),
            Err(e) => {
                let log = text(c.log_size, c.log, "nvrtcGetProgramLog").unwrap_or_default();
                let log = String::from_utf8_lossy(&log);
                Err(anyhow!("{e}\n{}", log.trim_end_matches('\0')))
            }
        };
        (c.destroy)(&mut program);
        Ok((result?, format!("{major}.{minor}")))
    }
}
