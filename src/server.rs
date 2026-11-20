//! `vibe-rations serve`: one machine's simulations, run for clients that come and go.
//!
//! The server holds runs by id and keeps them going whether or not anyone is connected; a run
//! lives in its own directory under the data directory (its configuration, model files,
//! status and outputs), so finished runs' results outlast the server. `vibe-rations run
//! --remote` and the web page both use its HTTP API. It listens on localhost only by default:
//! reach it through an SSH tunnel, which authenticates and encrypts, rather than exposing it.
//!
//! ```text
//! GET    /api/info                        version, threads, memory available
//! POST   /api/estimate                    body: TOML → the memory it needs (memory::Memory)
//! GET    /api/runs                        every run's status, newest first
//! POST   /api/runs?name=N                 body: TOML → a new run; `needs` lists files to upload
//! GET    /api/runs/ID                     its status
//! PUT    /api/runs/ID/files/PATH          body: a model file the configuration names
//! POST   /api/runs/ID/start?paused=B      sets it up (allocating), then runs it or holds it
//! POST   /api/runs/ID/pause, …/resume
//! DELETE /api/runs/ID                     stops it and deletes it with its files
//! GET    /api/runs/ID/section?field=F&plane=P&index=I&gain=G   PNG, with an X-Clip header
//! GET    /api/runs/ID/gather?component=C&rows=R                PNG
//! GET    /api/runs/ID/traces/C.npy        the traces of component C so far
//! GET    /api/runs/ID/config              its configuration
//! GET    /api/runs/ID/outputs             the files it wrote: [{path, bytes}]
//! GET    /api/runs/ID/outputs/PATH        one of them
//! GET    /…                               the web page, when --page names its files
//! ```

use std::collections::BTreeMap;
use std::collections::hash_map::RandomState;
use std::fs::{self, File};
use std::hash::{BuildHasher, Hasher};
use std::io;
use std::path::{Component as PathPart, Path, PathBuf};
use std::sync::atomic::{AtomicUsize, Ordering};
use std::sync::{Arc, Condvar, Mutex, MutexGuard};
use std::time::{Duration, Instant, SystemTime, UNIX_EPOCH};

use anyhow::{Context, Result, anyhow};
use serde::{Deserialize, Serialize};
use tiny_http::{Header, Method, Request, Response};

use crate::config::Config;
use crate::memory;
use crate::simulation::Simulation;
use crate::view::View;

/// Where a run is.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "lowercase")]
pub enum Phase {
    /// Waiting for its model files, or to be started.
    Created,
    /// Set up, not stepping.
    Paused,
    Running,
    /// Finished; its outputs are written.
    Done,
    Failed,
    /// The server stopped while it was set up; its wavefield is gone.
    Interrupted,
}

impl Phase {
    /// Whether the run will take no more steps.
    pub fn finished(self) -> bool {
        matches!(self, Phase::Done | Phase::Failed | Phase::Interrupted)
    }
}

/// A run, as clients see it.
#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct Status {
    pub id: String,
    /// The configuration's file name.
    pub name: String,
    pub phase: Phase,
    /// When it was created, in seconds since 1970.
    pub created: u64,
    pub step: usize,
    pub total_steps: usize,
    /// Simulated time reached, in seconds.
    pub time: f64,
    /// Mcell updates per second so far.
    pub speed: f64,
    /// Model files the configuration names that are still to be uploaded.
    pub needs: Vec<String>,
    pub error: Option<String>,
    /// The run's settings once it is set up: `run.json`'s fields and the receivers' positions.
    pub info: Option<serde_json::Value>,
    /// Whether its wavefield is in memory, so that sections and gathers can be drawn.
    pub live: bool,
    /// Nobody has asked for a picture of it in the last few seconds, so it steps in longer
    /// chunks, which saves the per-chunk work (unless a watched run shares the machine).
    #[serde(default)]
    pub headless: bool,
}

/// How to serve.
pub struct Options {
    /// Address to listen on; 127.0.0.1:8088 by default.
    pub listen: String,
    /// Where runs keep their files.
    pub data: PathBuf,
    /// The web page's files (web/dist), to serve at /.
    pub page: Option<PathBuf>,
    /// Print every request and its answer on stderr.
    pub log: bool,
}

struct Run {
    dir: PathBuf,
    /// What clients read: locked briefly, never while stepping, so status requests don't wait.
    shared: Mutex<Shared>,
    /// The wavefield and its stepping, held for a whole chunk. Lock it with [`Run::state`],
    /// never while holding `shared`.
    state: Mutex<State>,
    /// Requests waiting for `state`: the scheduler lets them in before its next chunk, since
    /// the mutex isn't fair and would otherwise go straight back to the scheduler.
    waiting: AtomicUsize,
}

struct Shared {
    status: Status,
    /// When a client last asked for a picture of it.
    pictured: Option<Instant>,
}

struct State {
    view: Option<View>,
    /// Steps per scheduling round, adapted to the chunk length wanted.
    chunk: usize,
    /// Cell updates and seconds spent, for the speed.
    work: (f64, f64),
}

/// A run nobody has asked for a picture of for this long is headless.
const WATCH_WINDOW: Duration = Duration::from_secs(5);

impl Shared {
    fn watched(&self) -> bool {
        self.pictured.is_some_and(|t| t.elapsed() < WATCH_WINDOW)
    }

    /// The status as clients see it now.
    fn status(&self) -> Status {
        Status {
            headless: !self.watched(),
            ..self.status.clone()
        }
    }
}

impl Run {
    fn new(dir: PathBuf, status: Status) -> Self {
        Run {
            dir,
            shared: Mutex::new(Shared {
                status,
                pictured: None,
            }),
            state: Mutex::new(State {
                view: None,
                chunk: 1,
                work: (0.0, 0.0),
            }),
            waiting: AtomicUsize::new(0),
        }
    }

    fn shared(&self) -> MutexGuard<'_, Shared> {
        self.shared.lock().unwrap()
    }

    /// The wavefield, ahead of the scheduler's next chunk.
    fn state(&self) -> MutexGuard<'_, State> {
        self.waiting.fetch_add(1, Ordering::SeqCst);
        let guard = self.state.lock().unwrap();
        self.waiting.fetch_sub(1, Ordering::SeqCst);
        guard
    }

    fn status(&self) -> Status {
        self.shared().status()
    }
}

struct Server {
    data: PathBuf,
    page: Option<PathBuf>,
    log: bool,
    runs: Mutex<BTreeMap<String, Arc<Run>>>,
    /// Signalled when a run starts or resumes.
    wake: Condvar,
}

/// Serves until the process is stopped.
pub fn serve(opts: Options) -> Result<()> {
    fs::create_dir_all(&opts.data).with_context(|| format!("creating {}", opts.data.display()))?;
    let server = Arc::new(Server {
        runs: Mutex::new(load_runs(&opts.data)?),
        data: opts.data,
        page: opts.page,
        log: opts.log,
        wake: Condvar::new(),
    });
    let http = tiny_http::Server::http(&opts.listen)
        .map_err(|e| anyhow!("listening on {}: {e}", opts.listen))?;
    // The address bound, which differs from the one asked for with port 0.
    let bound = http
        .server_addr()
        .to_ip()
        .map_or(opts.listen.clone(), |a| a.to_string());
    eprintln!(
        "serving on http://{bound}: {} run(s) in {}",
        server.runs.lock().unwrap().len(),
        server.data.display()
    );
    let s = server.clone();
    std::thread::spawn(move || schedule(&s));
    for request in http.incoming_requests() {
        let s = server.clone();
        std::thread::spawn(move || handle(&s, request));
    }
    Ok(())
}

/// The runs a previous server left: finished ones keep their outputs, set-up ones lost their
/// wavefield.
fn load_runs(data: &Path) -> Result<BTreeMap<String, Arc<Run>>> {
    let mut runs = BTreeMap::new();
    for entry in fs::read_dir(data)? {
        let dir = entry?.path();
        let Ok(text) = fs::read_to_string(dir.join("status.json")) else {
            continue;
        };
        let Ok(mut status) = serde_json::from_str::<Status>(&text) else {
            continue;
        };
        if matches!(status.phase, Phase::Paused | Phase::Running) {
            status.phase = Phase::Interrupted;
            status.error = Some("the server stopped during the run".into());
        }
        status.live = false;
        let id = status.id.clone();
        let run = Run::new(dir, status);
        run.save(&run.shared().status);
        runs.insert(id, Arc::new(run));
    }
    Ok(runs)
}

/// Steps a watched run's chunk takes, in seconds: short, so its pictures move smoothly and
/// requests wait little for the run's lock.
const WATCHED_CHUNK: f64 = 0.1;
/// And a headless run's, when no watched run shares the machine: long, so the work done once
/// per chunk (the blow-up check, the status) costs next to nothing.
const HEADLESS_CHUNK: f64 = 2.0;

/// What a run did with its turn.
enum Turn {
    Stepped,
    /// Left its turn to a request waiting for it.
    Yielded,
    /// Nothing to do.
    Idle,
}

/// Steps every running run a chunk at a time, in turn, sleeping while none is.
fn schedule(server: &Server) {
    loop {
        let all: Vec<Arc<Run>> = server.runs.lock().unwrap().values().cloned().collect();
        let watched = all.iter().any(|r| {
            let s = r.shared();
            s.status.phase == Phase::Running && s.watched()
        });
        let (mut stepped, mut yielded) = (false, false);
        for run in all {
            match run.step(watched) {
                Turn::Stepped => stepped = true,
                Turn::Yielded => yielded = true,
                Turn::Idle => {}
            }
        }
        if !stepped && yielded {
            std::thread::sleep(Duration::from_millis(1));
        } else if !stepped {
            let runs = server.runs.lock().unwrap();
            let _ = server.wake.wait_timeout(runs, Duration::from_millis(500));
        }
    }
}

impl Run {
    /// Takes one chunk of steps if the run is running and no request is waiting for it.
    /// `watched`: whether any running run is being watched, which keeps every chunk short.
    fn step(&self, watched: bool) -> Turn {
        let headless = {
            let sh = self.shared();
            if sh.status.phase != Phase::Running {
                return Turn::Idle;
            }
            !sh.watched()
        };
        if self.waiting.load(Ordering::SeqCst) > 0 {
            return Turn::Yielded;
        }
        let Ok(mut guard) = self.state.try_lock() else {
            return Turn::Yielded;
        };
        let State { view, chunk, work } = &mut *guard;
        let Some(v) = view.as_mut() else {
            let mut sh = self.shared();
            sh.status.phase = Phase::Failed;
            sh.status.error = Some("the run lost its wavefield".into());
            self.save(&sh.status);
            return Turn::Idle;
        };
        let snapshots = self.dir.join("out/snapshots");
        let snapshots = (!v.sim.snapshots.is_empty()).then_some(snapshots.as_path());
        let start = Instant::now();
        let result = v.advance(*chunk, snapshots);
        let seconds = start.elapsed().as_secs_f64();
        if let Ok(taken) = result {
            work.0 += taken as f64 * v.sim.solver.grid.interior_len() as f64;
            work.1 += seconds;
        }
        let target = if headless && !watched {
            HEADLESS_CHUNK
        } else {
            WATCHED_CHUNK
        };
        *chunk =
            ((*chunk as f64 * target / seconds.max(1e-3)).round() as usize).clamp(1, 1_000_000);
        let outcome = match result {
            Err(e) => Err(e),
            Ok(_) if v.sim.remaining() == 0 => {
                v.sim.write_outputs(&self.dir.join("out")).map(|_| true)
            }
            Ok(_) => Ok(false),
        };
        let mut sh = self.shared();
        let status = &mut sh.status;
        status.step = v.sim.solver.step;
        status.time = v.sim.solver.step as f64 * v.sim.solver.dt;
        status.speed = if work.1 > 0.0 {
            work.0 / work.1 / 1e6
        } else {
            0.0
        };
        match outcome {
            Err(e) => {
                status.phase = Phase::Failed;
                status.error = Some(format!("{e:#}"));
                self.save(status);
            }
            Ok(true) => {
                status.phase = Phase::Done;
                self.save(status);
            }
            Ok(false) => {}
        }
        Turn::Stepped
    }

    /// Writes the status where the next server will find it.
    fn save(&self, status: &Status) {
        if let Ok(text) = serde_json::to_string_pretty(status) {
            let _ = fs::write(self.dir.join("status.json"), text);
        }
    }
}

/// An HTTP error: status code and message.
struct Error(u16, String);

impl From<anyhow::Error> for Error {
    fn from(e: anyhow::Error) -> Self {
        Error(400, format!("{e:#}"))
    }
}

fn not_found(what: &str) -> Error {
    Error(404, format!("no {what}"))
}

fn conflict(message: impl Into<String>) -> Error {
    Error(409, message.into())
}

/// What to send back.
enum Reply {
    Json(serde_json::Value),
    Bytes(&'static str, Vec<u8>, Vec<(&'static str, String)>),
    File(File, &'static str),
    Empty,
}

fn handle(server: &Server, mut request: Request) {
    let url = request.url().to_string();
    let (path, query) = url.split_once('?').unwrap_or((&url, ""));
    let parts: Vec<String> = path
        .split('/')
        .filter(|p| !p.is_empty())
        .map(decode)
        .collect();
    let parts: Vec<&str> = parts.iter().map(String::as_str).collect();
    let query = parse_query(query);
    // HEAD is GET without the body, which tiny_http leaves out.
    let method = match request.method() {
        Method::Head => Method::Get,
        m => m.clone(),
    };
    let reply = route(server, &method, &parts, &query, &mut request);
    if server.log {
        let code = match &reply {
            Ok(Reply::Empty) => 204,
            Ok(_) => 200,
            Err(Error(code, _)) => *code,
        };
        eprintln!("{} {url} -> {code}", request.method());
    }
    let _ = match reply {
        Ok(Reply::Json(v)) => respond(
            request,
            200,
            "application/json",
            v.to_string().into_bytes(),
            vec![],
        ),
        Ok(Reply::Bytes(kind, body, headers)) => respond(request, 200, kind, body, headers),
        Ok(Reply::File(file, kind)) => {
            let response = with_headers(Response::from_file(file), kind, vec![]);
            request.respond(response)
        }
        Ok(Reply::Empty) => respond(request, 204, "text/plain", vec![], vec![]),
        Err(Error(code, message)) => respond(
            request,
            code,
            "application/json",
            serde_json::json!({ "error": message })
                .to_string()
                .into_bytes(),
            vec![],
        ),
    };
}

fn respond(
    request: Request,
    code: u16,
    kind: &'static str,
    body: Vec<u8>,
    headers: Vec<(&'static str, String)>,
) -> io::Result<()> {
    let response = with_headers(Response::from_data(body), kind, headers).with_status_code(code);
    request.respond(response)
}

/// The content type, the extra headers, and those every response gets: cross-origin isolation
/// (the page's threads need it), no caching, no sniffing.
fn with_headers<R: io::Read>(
    mut response: Response<R>,
    kind: &str,
    extra: Vec<(&'static str, String)>,
) -> Response<R> {
    let fixed = [
        ("Content-Type", kind.to_string()),
        ("Cross-Origin-Opener-Policy", "same-origin".into()),
        ("Cross-Origin-Embedder-Policy", "require-corp".into()),
        ("Cache-Control", "no-cache".into()),
        ("X-Content-Type-Options", "nosniff".into()),
    ];
    for (name, value) in fixed.into_iter().chain(extra) {
        if let Ok(h) = Header::from_bytes(name.as_bytes(), value.as_bytes()) {
            response.add_header(h);
        }
    }
    response
}

fn route(
    server: &Server,
    method: &Method,
    parts: &[&str],
    query: &BTreeMap<String, String>,
    request: &mut Request,
) -> Result<Reply, Error> {
    use Method::{Delete, Get, Post, Put};
    match (method, parts) {
        (Get, ["api", "info"]) => Ok(Reply::Json(serde_json::json!({
            "version": env!("CARGO_PKG_VERSION"),
            "threads": rayon::current_num_threads(),
            "available": memory::available(),
        }))),
        (Post, ["api", "estimate"]) => {
            let cfg = Config::from_toml(&body_text(request)?)?;
            let m = memory::estimate(&cfg, &server.data)?;
            let mut value = serde_json::to_value(m).expect("Memory serialises");
            value["limit"] = serde_json::json!(memory::available());
            Ok(Reply::Json(value))
        }
        (Get, ["api", "runs"]) => {
            let runs: Vec<Arc<Run>> = server.runs.lock().unwrap().values().cloned().collect();
            let mut list: Vec<Status> = runs.iter().map(|r| status(r)).collect();
            list.sort_by(|a, b| b.created.cmp(&a.created).then(b.id.cmp(&a.id)));
            Ok(Reply::Json(serde_json::to_value(list).expect("serialises")))
        }
        (Post, ["api", "runs"]) => {
            let name = query.get("name").map_or("run.toml", String::as_str);
            create(server, name, &body_text(request)?).map(|s| json(&s))
        }
        (Get, ["api", "runs", id]) => Ok(json(&status(&*find(server, id)?))),
        (Delete, ["api", "runs", id]) => {
            let run = server
                .runs
                .lock()
                .unwrap()
                .remove(*id)
                .ok_or_else(|| not_found("such run"))?;
            // Waits for the scheduler to let go of it, then frees its memory.
            run.state().view = None;
            fs::remove_dir_all(&run.dir).map_err(|e| Error(500, e.to_string()))?;
            Ok(Reply::Empty)
        }
        (Put, ["api", "runs", id, "files", rest @ ..]) => {
            upload(&*find(server, id)?, &rest.join("/"), request).map(|s| json(&s))
        }
        (Post, ["api", "runs", id, "start"]) => {
            start(server, &find(server, id)?, flag(query, "paused")).map(|s| json(&s))
        }
        (Post, ["api", "runs", id, action @ ("pause" | "resume")]) => {
            let run = find(server, id)?;
            let mut sh = run.shared();
            let (from, to) = match *action {
                "pause" => (Phase::Running, Phase::Paused),
                _ => (Phase::Paused, Phase::Running),
            };
            if sh.status.phase != from {
                return Err(conflict(
                    format!("the run is {:?}, not {from:?}", sh.status.phase).to_lowercase(),
                ));
            }
            sh.status.phase = to;
            run.save(&sh.status);
            let s = sh.status();
            drop(sh);
            server.wake.notify_all();
            Ok(json(&s))
        }
        (Get, ["api", "runs", id, "section"]) => {
            let run = find(server, id)?;
            run.shared().pictured = Some(Instant::now());
            let mut st = run.state();
            let view = live(&mut st)?;
            let index: usize = number(query, "index", 0.0)? as usize;
            let gain = number(query, "gain", 1.0)? as f32;
            let field = query.get("field").map_or("vz", String::as_str);
            let plane = query.get("plane").map_or("xz", String::as_str);
            let frame = view.section(field, plane, index, gain)?;
            Ok(Reply::Bytes(
                "image/png",
                frame.png()?,
                vec![("X-Clip", frame.clip.to_string())],
            ))
        }
        (Get, ["api", "runs", id, "gather"]) => {
            let run = find(server, id)?;
            run.shared().pictured = Some(Instant::now());
            let mut st = run.state();
            let view = live(&mut st)?;
            let rows = number(query, "rows", 600.0)? as usize;
            let component = query.get("component").map_or("vz", String::as_str);
            Ok(Reply::Bytes(
                "image/png",
                view.gather(component, rows)?.png()?,
                vec![],
            ))
        }
        (Get, ["api", "runs", id, "traces", file]) => {
            let run = find(server, id)?;
            let mut st = run.state();
            let view = live(&mut st)?;
            let component = file.strip_suffix(".npy").unwrap_or(file);
            Ok(Reply::Bytes(
                "application/octet-stream",
                view.traces_npy(component)?,
                vec![],
            ))
        }
        (Get, ["api", "runs", id, "config"]) => {
            let run = find(server, id)?;
            let text =
                fs::read(run.dir.join("config.toml")).map_err(|e| Error(500, e.to_string()))?;
            Ok(Reply::Bytes("text/plain; charset=utf-8", text, vec![]))
        }
        (Get, ["api", "runs", id, "outputs"]) => {
            let run = find(server, id)?;
            let out = run.dir.join("out");
            let mut files = Vec::new();
            list_files(&out, &out, &mut files);
            Ok(Reply::Json(
                serde_json::to_value(files).expect("serialises"),
            ))
        }
        (Get, ["api", "runs", id, "outputs", rest @ ..]) => {
            let run = find(server, id)?;
            let path = inside(&run.dir.join("out"), &rest.join("/"))?;
            let file = File::open(&path).map_err(|_| not_found("such output"))?;
            Ok(Reply::File(file, "application/octet-stream"))
        }
        (Get, ["api", ..]) | (_, ["api", ..]) => Err(not_found("such API call")),
        (Get, rest) => {
            let page = server.page.as_ref().ok_or_else(|| not_found("page here"))?;
            let rel = if rest.is_empty() {
                "index.html".to_string()
            } else {
                rest.join("/")
            };
            let path = inside(page, &rel)?;
            let file = File::open(&path).map_err(|_| not_found("such file"))?;
            Ok(Reply::File(file, content_type(&path)))
        }
        _ => Err(Error(405, "method not allowed".into())),
    }
}

fn json(status: &Status) -> Reply {
    Reply::Json(serde_json::to_value(status).expect("Status serialises"))
}

fn find(server: &Server, id: &str) -> Result<Arc<Run>, Error> {
    server
        .runs
        .lock()
        .unwrap()
        .get(id)
        .cloned()
        .ok_or_else(|| not_found("such run"))
}

fn status(run: &Run) -> Status {
    run.status()
}

fn live(state: &mut State) -> Result<&mut View, Error> {
    state.view.as_mut().ok_or_else(|| {
        conflict("this run's wavefield is no longer in memory; its outputs can still be fetched")
    })
}

fn body_text(request: &mut Request) -> Result<String, Error> {
    let mut text = String::new();
    io::Read::read_to_string(request.as_reader(), &mut text)
        .map_err(|e| Error(400, format!("reading the request: {e}")))?;
    Ok(text)
}

/// A query flag: true for `true` or `1`.
fn flag(query: &BTreeMap<String, String>, key: &str) -> bool {
    query.get(key).is_some_and(|v| v == "true" || v == "1")
}

fn number(query: &BTreeMap<String, String>, key: &str, default: f64) -> Result<f64, Error> {
    match query.get(key) {
        None => Ok(default),
        Some(v) => v
            .parse()
            .map_err(|_| Error(400, format!("{key} must be a number"))),
    }
}

/// A new run from a configuration: its directory, configuration and status.
fn create(server: &Server, name: &str, toml: &str) -> Result<Status, Error> {
    let cfg = Config::from_toml(toml)?;
    let needs: Vec<String> = match &cfg.model.files {
        Some(f) => [&f.vp, &f.vs, &f.rho]
            .iter()
            .map(|p| p.to_string_lossy().into_owned())
            .collect(),
        None => vec![],
    };
    for n in &needs {
        relative(n)?;
    }
    let stem: String = Path::new(name)
        .file_stem()
        .map(|s| s.to_string_lossy())
        .unwrap_or_default()
        .chars()
        .filter(|c| c.is_ascii_alphanumeric() || *c == '-' || *c == '_')
        .take(32)
        .collect();
    let created = SystemTime::now()
        .duration_since(UNIX_EPOCH)
        .map_or(0, |d| d.as_secs());
    let random = RandomState::new().build_hasher().finish() & 0xffff;
    let id = format!(
        "{}-{created}-{random:04x}",
        if stem.is_empty() { "run" } else { &stem }
    );
    let dir = server.data.join(&id);
    fs::create_dir_all(&dir).map_err(|e| Error(500, e.to_string()))?;
    fs::write(dir.join("config.toml"), toml).map_err(|e| Error(500, e.to_string()))?;
    let status = Status {
        id: id.clone(),
        name: name.to_string(),
        phase: Phase::Created,
        created,
        step: 0,
        total_steps: 0,
        time: 0.0,
        speed: 0.0,
        needs,
        error: None,
        info: None,
        live: false,
        headless: false,
    };
    let run = Run::new(dir, status.clone());
    run.save(&status);
    server.runs.lock().unwrap().insert(id, Arc::new(run));
    Ok(status)
}

/// Stores a model file the configuration names.
fn upload(run: &Run, name: &str, request: &mut Request) -> Result<Status, Error> {
    if !run.shared().status.needs.iter().any(|n| n == name) {
        return Err(conflict(format!("the run doesn't need {name}")));
    }
    let path = run.dir.join(relative(name)?);
    if let Some(parent) = path.parent() {
        fs::create_dir_all(parent).map_err(|e| Error(500, e.to_string()))?;
    }
    let mut file = File::create(&path).map_err(|e| Error(500, e.to_string()))?;
    io::copy(request.as_reader(), &mut file)
        .map_err(|e| Error(400, format!("receiving {name}: {e}")))?;
    let mut sh = run.shared();
    sh.status.needs.retain(|n| n != name);
    run.save(&sh.status);
    Ok(sh.status())
}

/// Sets a run up, if it fits in the memory available, and runs it or holds it at step 0.
fn start(server: &Server, run: &Arc<Run>, paused: bool) -> Result<Status, Error> {
    // Held throughout, so two starts can't both set it up.
    let mut st = run.state();
    {
        let sh = run.shared();
        if sh.status.phase != Phase::Created {
            return Err(conflict(format!(
                "the run has been started already ({:?})",
                sh.status.phase
            )));
        }
        if !sh.status.needs.is_empty() {
            return Err(conflict(format!(
                "upload {} first",
                sh.status.needs.join(", ")
            )));
        }
    }
    let text =
        fs::read_to_string(run.dir.join("config.toml")).map_err(|e| Error(500, e.to_string()))?;
    let cfg = Config::from_toml(&text)?;
    let need = memory::estimate(&cfg, &run.dir)?.peak;
    if let Some(mut free) = memory::available() {
        if need > free {
            free += evict(server, run, need - free);
        }
        if need > free {
            return Err(Error(
                507,
                format!(
                    "the run needs {} and only {} is available",
                    memory::human(need),
                    memory::human(free)
                ),
            ));
        }
    }
    let sim = Simulation::new(&cfg, &run.dir)?;
    let out = run.dir.join("out");
    fs::create_dir_all(out.join("snapshots")).map_err(|e| Error(500, e.to_string()))?;
    let view = View::new(sim);
    let mut sh = run.shared();
    sh.status.info = Some(view.info());
    sh.status.total_steps = view.sim.steps;
    sh.status.live = true;
    sh.status.phase = if paused {
        Phase::Paused
    } else {
        Phase::Running
    };
    st.view = Some(view);
    run.save(&sh.status);
    let s = sh.status();
    drop(sh);
    drop(st);
    server.wake.notify_all();
    Ok(s)
}

/// Frees the wavefields of finished runs, oldest first, until about `wanted` bytes are free
/// or none is left; returns what it freed.
fn evict(server: &Server, keep: &Arc<Run>, wanted: u64) -> u64 {
    let mut runs: Vec<Arc<Run>> = server.runs.lock().unwrap().values().cloned().collect();
    runs.sort_by_key(|r| r.shared().status.created);
    let mut freed = 0;
    for r in runs {
        if freed >= wanted || Arc::ptr_eq(&r, keep) {
            continue;
        }
        if !r.shared().status.phase.finished() {
            continue;
        }
        let Ok(mut st) = r.state.try_lock() else {
            continue;
        };
        if let Some(v) = st.view.take() {
            freed += v.sim.info.memory.resident;
            let mut sh = r.shared();
            sh.status.live = false;
            r.save(&sh.status);
        }
    }
    freed
}

/// A relative path with no way out of the directory it is joined to.
fn relative(path: &str) -> Result<PathBuf, Error> {
    let p = Path::new(path);
    let ok = !path.is_empty() && p.components().all(|c| matches!(c, PathPart::Normal(_)));
    if ok {
        Ok(p.to_path_buf())
    } else {
        Err(Error(
            400,
            format!("{path} must be a relative path inside the run"),
        ))
    }
}

/// `rel` under `root`, refusing anything that would leave it.
fn inside(root: &Path, rel: &str) -> Result<PathBuf, Error> {
    Ok(root.join(relative(rel)?))
}

fn list_files(root: &Path, dir: &Path, out: &mut Vec<serde_json::Value>) {
    let Ok(entries) = fs::read_dir(dir) else {
        return;
    };
    let mut entries: Vec<_> = entries.flatten().collect();
    entries.sort_by_key(|e| e.path());
    for e in entries {
        let path = e.path();
        if path.is_dir() {
            list_files(root, &path, out);
        } else if let (Ok(rel), Ok(meta)) = (path.strip_prefix(root), e.metadata()) {
            out.push(serde_json::json!({
                "path": rel.to_string_lossy(),
                "bytes": meta.len(),
            }));
        }
    }
}

fn content_type(path: &Path) -> &'static str {
    match path.extension().and_then(|e| e.to_str()) {
        Some("html") => "text/html; charset=utf-8",
        Some("js") => "text/javascript",
        Some("css") => "text/css; charset=utf-8",
        Some("wasm") => "application/wasm",
        Some("json") => "application/json",
        Some("png") => "image/png",
        Some("svg") => "image/svg+xml",
        Some("toml") | Some("txt") => "text/plain; charset=utf-8",
        _ => "application/octet-stream",
    }
}

fn parse_query(query: &str) -> BTreeMap<String, String> {
    query
        .split('&')
        .filter(|p| !p.is_empty())
        .map(|p| {
            let (k, v) = p.split_once('=').unwrap_or((p, ""));
            (decode(k), decode(v))
        })
        .collect()
}

/// Percent-decoding, with + as a space in query strings.
fn decode(s: &str) -> String {
    let b = s.as_bytes();
    let mut out = Vec::with_capacity(b.len());
    let mut i = 0;
    while i < b.len() {
        let hex = (b[i] == b'%' && i + 3 <= b.len())
            .then(|| std::str::from_utf8(&b[i + 1..i + 3]).ok())
            .flatten()
            .and_then(|h| u8::from_str_radix(h, 16).ok());
        match (hex, b[i]) {
            (Some(byte), _) => {
                out.push(byte);
                i += 3;
                continue;
            }
            (None, b'+') => out.push(b' '),
            (None, c) => out.push(c),
        }
        i += 1;
    }
    String::from_utf8_lossy(&out).into_owned()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn decodes_percent_escapes() {
        assert_eq!(decode("a%20b+c"), "a b c");
        assert_eq!(decode("100%"), "100%");
        assert_eq!(decode("%zz"), "%zz");
        assert_eq!(parse_query("a=1&b=x%2Fy")["b"], "x/y");
    }

    #[test]
    fn keeps_paths_inside() {
        assert!(relative("models/vp.npy").is_ok());
        for bad in ["../etc/passwd", "/etc/passwd", "a/../../b", ""] {
            assert!(relative(bad).is_err(), "{bad}");
        }
    }

    #[test]
    fn serves_the_page_with_isolation_headers() {
        assert_eq!(content_type(Path::new("x.wasm")), "application/wasm");
        assert_eq!(content_type(Path::new("app.js")), "text/javascript");
    }
}
