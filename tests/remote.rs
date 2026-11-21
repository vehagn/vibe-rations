//! `vibe-rations serve` and its client, as separate processes talking HTTP: remote runs give
//! the same results as local ones, outlast their clients, and outlast the server's restarts.

use std::fs;
use std::io::{BufRead, BufReader};
use std::path::{Path, PathBuf};
use std::process::{Child, Command, Stdio};
use std::time::{Duration, Instant};

use vibe_rations::npy;
use vibe_rations::remote::{self, Remote};
use vibe_rations::server::Phase;

const BIN: &str = env!("CARGO_BIN_EXE_vibe-rations");

/// A fresh, empty directory for one test.
fn scratch(name: &str) -> PathBuf {
    let dir =
        std::env::temp_dir().join(format!("vibe-rations-remote-{name}-{}", std::process::id()));
    let _ = fs::remove_dir_all(&dir);
    fs::create_dir_all(&dir).unwrap();
    dir
}

/// A server process, stopped when dropped.
struct Server {
    child: Child,
    url: String,
}

impl Server {
    /// Starts a server on a port of its own choosing, which it names on its first line.
    fn start(data: &Path, page: Option<&Path>) -> Server {
        let mut cmd = Command::new(BIN);
        cmd.args(["serve", "--listen", "127.0.0.1:0", "--data"])
            .arg(data)
            .stderr(Stdio::piped());
        if let Some(p) = page {
            cmd.arg("--page").arg(p);
        }
        let mut child = cmd.spawn().unwrap();
        let mut line = String::new();
        BufReader::new(child.stderr.take().unwrap())
            .read_line(&mut line)
            .unwrap();
        let url = line
            .strip_prefix("serving on ")
            .and_then(|rest| rest.split(':').take(3).collect::<Vec<_>>().join(":").into())
            .unwrap_or_else(|| panic!("the server said {line:?}"));
        Server { child, url }
    }

    fn client(&self) -> Remote {
        remote::connect(&self.url, 0).unwrap()
    }
}

impl Drop for Server {
    fn drop(&mut self) {
        let _ = self.child.kill();
        let _ = self.child.wait();
    }
}

const CONFIG: &str = r#"
[grid]
points = [24, 24, 20]
spacing = 10.0
[time]
duration = 0.05
[boundary]
pml_width = 6
[[model.layer]]
top = 0.0
vp = 2000.0
vs = 1000.0
rho = 2000.0
[[source]]
position = [120.0, 120.0, 60.0]
explosion = 1e12
wavelet = { kind = "ricker", frequency = 10.0 }
[receivers]
components = ["vz", "p"]
[[receivers.line]]
from = [90.0, 120.0, 0.0]
to = [150.0, 120.0, 0.0]
count = 4
[[snapshot]]
field = "vz"
plane = "xz"
at = 120.0
interval = 0.025
formats = ["png", "npy"]
"#;

/// The traces a local run of `config` (in `dir`) writes.
fn local_traces(dir: &Path, config: &Path) -> Vec<f32> {
    let out = dir.join("local-out");
    let status = Command::new(BIN)
        .args(["run", "--quiet", "--out"])
        .arg(&out)
        .arg(config)
        .stderr(Stdio::null())
        .status()
        .unwrap();
    assert!(status.success());
    npy::read(&out.join("vz.npy")).unwrap().1
}

#[test]
fn a_remote_run_gives_what_a_local_one_does() {
    let dir = scratch("same");
    let config = dir.join("tiny.toml");
    fs::write(&config, CONFIG).unwrap();
    let server = Server::start(&dir.join("data"), None);
    let r = server.client();
    let s = r.submit(&config).unwrap();
    assert!(s.info.is_some());
    let end = r.follow(&s.id, false).unwrap();
    assert_eq!(end.phase, Phase::Done, "{:?}", end.error);
    let out = dir.join("remote-out");
    let n = r.download(&s.id, &out).unwrap();
    assert!(n >= 6, "{n} files");
    assert!(out.join("snapshots").read_dir().unwrap().count() >= 2);
    assert_eq!(
        npy::read(&out.join("vz.npy")).unwrap().1,
        local_traces(&dir, &config)
    );
    // The CLI does the same in one go.
    let cli_out = dir.join("cli-out");
    let ok = Command::new(BIN)
        .args(["run", "--quiet", "--remote", &server.url, "--out"])
        .arg(&cli_out)
        .arg(&config)
        .stderr(Stdio::null())
        .status()
        .unwrap();
    assert!(ok.success());
    assert_eq!(
        fs::read(cli_out.join("vz.npy")).unwrap(),
        fs::read(out.join("vz.npy")).unwrap()
    );
    fs::remove_dir_all(dir).unwrap();
}

#[test]
fn uploads_model_files() {
    let dir = scratch("files");
    let (nx, ny, nz) = (24, 24, 20);
    let n = nx * ny * nz;
    fs::create_dir_all(dir.join("model")).unwrap();
    for (name, v) in [("vp", 2000.0f32), ("vs", 1000.0), ("rho", 2000.0)] {
        npy::write(
            &dir.join(format!("model/{name}.npy")),
            &[nz, ny, nx],
            &vec![v; n],
        )
        .unwrap();
    }
    let text = CONFIG.replace(
        "[[model.layer]]\ntop = 0.0\nvp = 2000.0\nvs = 1000.0\nrho = 2000.0",
        "[model.files]\nvp = \"model/vp.npy\"\nvs = \"model/vs.npy\"\nrho = \"model/rho.npy\"",
    );
    assert!(text.contains("model/vp.npy"));
    let config = dir.join("files.toml");
    fs::write(&config, &text).unwrap();
    let server = Server::start(&dir.join("data"), None);
    let r = server.client();
    let s = r.submit(&config).unwrap();
    assert!(s.needs.is_empty());
    assert_eq!(r.follow(&s.id, false).unwrap().phase, Phase::Done);
    let out = dir.join("remote-out");
    r.download(&s.id, &out).unwrap();
    assert_eq!(
        npy::read(&out.join("vz.npy")).unwrap().1,
        local_traces(&dir, &config)
    );
    fs::remove_dir_all(dir).unwrap();
}

/// What the page does: create a run held at step 0, look at it, run it, pause and resume it,
/// and pick it up again from another connection; and what a restart leaves.
#[test]
fn runs_outlast_their_clients_and_the_server() {
    let dir = scratch("outlast");
    let page = dir.join("page");
    fs::create_dir_all(&page).unwrap();
    fs::write(page.join("index.html"), "<!doctype html><title>t</title>").unwrap();
    let data = dir.join("data");
    let agent: ureq::Agent = ureq::Agent::config_builder()
        .http_status_as_error(false)
        .build()
        .into();
    let call = |method: &str, url: &str, body: &str| -> (u16, String) {
        let mut r = match method {
            "GET" => agent.get(url).call(),
            "DELETE" => agent.delete(url).call(),
            "PUT" => agent.put(url).send(body),
            _ => agent.post(url).send(body),
        }
        .unwrap();
        let code = r.status().as_u16();
        (code, r.body_mut().read_to_string().unwrap_or_default())
    };
    let (held, deleted);
    {
        let server = Server::start(&data, Some(&page));
        let u = &server.url;
        // The page, with the headers its threads need.
        let page_reply = agent.get(u).call().unwrap();
        assert_eq!(
            page_reply.headers()["cross-origin-embedder-policy"],
            "require-corp"
        );
        assert_eq!(call("GET", &format!("{u}/../etc/passwd"), "").0, 400);
        // Created, held at step 0, pictured.
        let (code, body) = call("POST", &format!("{u}/api/runs?name=tiny.toml"), CONFIG);
        assert_eq!(code, 200, "{body}");
        let id = serde_json::from_str::<serde_json::Value>(&body).unwrap()["id"]
            .as_str()
            .unwrap()
            .to_string();
        let (code, body) = call("POST", &format!("{u}/api/runs/{id}/start?paused=true"), "");
        assert_eq!(code, 200, "{body}");
        assert!(body.contains("\"phase\":\"paused\""), "{body}");
        // Nobody has asked for a picture yet: headless.
        assert!(body.contains("\"headless\":true"), "{body}");
        let png = agent
            .get(&format!(
                "{u}/api/runs/{id}/section?field=vz&plane=xz&index=12&gain=1"
            ))
            .call()
            .unwrap();
        assert!(png.headers().contains_key("x-clip"));
        // Asked for a picture: watched.
        let (_, body) = call("GET", &format!("{u}/api/runs/{id}"), "");
        assert!(body.contains("\"headless\":false"), "{body}");
        assert_eq!(call("POST", &format!("{u}/api/runs/{id}/start"), "").0, 409);
        assert_eq!(call("POST", &format!("{u}/api/runs/{id}/pause"), "").0, 409);
        // Run, then follow from a second connection to the end.
        assert_eq!(
            call("POST", &format!("{u}/api/runs/{id}/resume"), "").0,
            200
        );
        let other = server.client();
        assert_eq!(other.follow(&id, false).unwrap().phase, Phase::Done);
        let (code, gather) = call(
            "GET",
            &format!("{u}/api/runs/{id}/gather?component=p&rows=50"),
            "",
        );
        assert_eq!(code, 200, "{gather}");
        // Another run held at step 0, and one deleted.
        let (_, body) = call("POST", &format!("{u}/api/runs?name=held.toml"), CONFIG);
        held = serde_json::from_str::<serde_json::Value>(&body).unwrap()["id"]
            .as_str()
            .unwrap()
            .to_string();
        call(
            "POST",
            &format!("{u}/api/runs/{held}/start?paused=true"),
            "",
        );
        let (_, body) = call("POST", &format!("{u}/api/runs?name=gone.toml"), CONFIG);
        deleted = serde_json::from_str::<serde_json::Value>(&body).unwrap()["id"]
            .as_str()
            .unwrap()
            .to_string();
        assert_eq!(
            call("DELETE", &format!("{u}/api/runs/{deleted}"), "").0,
            204
        );
        assert!(!data.join(&deleted).exists());
        assert_eq!(call("GET", &format!("{u}/api/runs/nonesuch"), "").0, 404);
        // A file the configuration doesn't name, or one outside the run, is refused.
        let (_, body) = call("POST", &format!("{u}/api/runs?name=x.toml"), CONFIG);
        let x = serde_json::from_str::<serde_json::Value>(&body).unwrap()["id"]
            .as_str()
            .unwrap()
            .to_string();
        assert_eq!(
            call(
                "PUT",
                &format!("{u}/api/runs/{x}/files/..%2F..%2Fevil"),
                "x"
            )
            .0,
            409
        );
    }
    // After a restart: the finished run's results are still there, the held one is gone.
    let server = Server::start(&data, None);
    let r = server.client();
    let runs = r.runs().unwrap();
    assert_eq!(runs.len(), 3, "{runs:?}");
    let done = runs.iter().find(|s| s.name == "tiny.toml").unwrap();
    assert_eq!(done.phase, Phase::Done);
    assert!(!done.live);
    assert!(r.download(&done.id, &dir.join("after")).unwrap() >= 6);
    assert_eq!(r.status(&held).unwrap().phase, Phase::Interrupted);
    fs::remove_dir_all(dir).unwrap();
}

/// A headless run steps in long chunks, but its status answers at once throughout: std's
/// Mutex isn't fair, and with the status under the wavefield's lock a request once waited out
/// a whole run.
#[test]
fn status_answers_promptly_during_a_headless_run() {
    let dir = scratch("prompt");
    let config = dir.join("long.toml");
    let text = CONFIG
        .replace("points = [24, 24, 20]", "points = [72, 72, 64]")
        .replace("duration = 0.05", "duration = 2.0")
        .replace("[[snapshot]]", "[[snapshot.ignored]]");
    // No snapshots: just a run long enough for several headless chunks.
    let text = text[..text.find("[[snapshot.ignored]]").unwrap()].to_string();
    fs::write(&config, text).unwrap();
    let server = Server::start(&dir.join("data"), None);
    let r = server.client();
    let id = r.submit(&config).unwrap().id;
    let mut steps = Vec::new();
    for _ in 0..6 {
        std::thread::sleep(Duration::from_millis(300));
        let start = Instant::now();
        let s = r.status(&id).unwrap();
        let waited = start.elapsed();
        assert!(
            waited < Duration::from_millis(300),
            "status took {waited:?}"
        );
        if s.phase == Phase::Running {
            assert!(s.headless, "nobody asked for pictures");
            steps.push(s.step);
        }
    }
    assert!(steps.len() >= 2, "the run should still be going: {steps:?}");
    assert_eq!(r.follow(&id, false).unwrap().phase, Phase::Done);
    fs::remove_dir_all(dir).unwrap();
}

/// The start of a small synthetic inversion (two shots, the top held fixed) and its truth, with
/// a faster sphere, sampled alike: what the web page sends.
fn inversion_pair() -> (String, String) {
    let config = |sphere: &str, inversion: &str| {
        format!(
            r#"
[grid]
points = [24, 20, 18]
spacing = 10.0
[time]
duration = 0.15
dt = 5e-4
[boundary]
pml_width = 4
pml_vp_max = 2700.0
[scheme]
coefficients = "taylor"
[[model.layer]]
top = 0.0
vp = 2400.0
vs = 1300.0
rho = 2100.0
{sphere}
[[source]]
position = [80.0, 90.0, 30.0]
force = [0.0, 0.0, 1e9]
wavelet = {{ kind = "ricker", frequency = 15.0 }}
[[source]]
position = [150.0, 100.0, 30.0]
force = [1e9, 0.0, 0.0]
wavelet = {{ kind = "ricker", frequency = 15.0 }}
[receivers]
components = ["vx", "vz"]
[[receivers.line]]
from = [60.0, 70.0, 0.0]
to = [170.0, 70.0, 0.0]
count = 12
{inversion}
"#
        )
    };
    let sphere = "[[model.sphere]]\ncenter = [115.0, 95.0, 100.0]\nradius = 25.0\n\
                  vp = 2600.0\nvs = 1450.0\nrho = 2100.0";
    (
        config(
            "",
            "[inversion]\nobserved = \"observed\"\niterations = 2\nfixed_above = 40.0",
        ),
        config(sphere, ""),
    )
}

/// An inversion on the server ends where the same inversion here does, bit for bit; it
/// pictures and writes what the page and `invert` expect, and lists as one.
#[test]
fn an_inversion_on_the_server_gives_what_one_here_does() {
    use vibe_rations::config::Config;
    use vibe_rations::server::Kind;
    use vibe_rations::view::InversionView;

    let dir = scratch("inversion");
    let (start, truth) = inversion_pair();
    let server = Server::start(&dir.join("data"), None);
    let u = &server.url;
    let agent: ureq::Agent = ureq::Agent::config_builder()
        .http_status_as_error(false)
        .build()
        .into();
    let post = |path: &str, body: &str| {
        let mut r = agent.post(&format!("{u}{path}")).send(body).unwrap();
        (r.status().as_u16(), r.body_mut().read_to_string().unwrap())
    };
    // Its memory, more than a forward run's; and no [inversion], no inversion.
    let peak = |body: &str| {
        serde_json::from_str::<serde_json::Value>(body).unwrap()["peak"]
            .as_u64()
            .unwrap()
    };
    let (code, inversion) = post("/api/estimate?kind=inversion", &start);
    assert_eq!(code, 200, "{inversion}");
    assert!(peak(&inversion) > peak(&post("/api/estimate", &start).1));
    assert_eq!(post("/api/runs?kind=inversion", &truth).0, 400);

    let r = server.client();
    let s = r.submit_inversion("sphere.toml", &start, &truth).unwrap();
    assert_eq!((s.kind, s.total_steps), (Kind::Inversion, 5));
    let end = r.follow(&s.id, false).unwrap();
    assert_eq!(end.phase, Phase::Done, "{:?}", end.error);
    let status = end.inversion.unwrap();
    assert_eq!(status["history"].as_array().unwrap().len(), 2);
    assert_eq!(end.step, 5);
    assert!(remote::progress_line(&r.status(&s.id).unwrap()).contains("5/5 units"));

    // The same here.
    let base = Path::new(".");
    let (start_cfg, truth_cfg) = (
        Config::from_toml(&start).unwrap(),
        Config::from_toml(&truth).unwrap(),
    );
    let mut here = InversionView::new(&start_cfg, &truth_cfg, base).unwrap();
    while here.advance().unwrap() {}
    assert_eq!(here.status()["misfit"], status["misfit"]);
    let out = dir.join("out");
    assert!(r.download(&s.id, &out).unwrap() >= 4);
    assert_eq!(
        fs::read(out.join("vs.npy")).unwrap(),
        here.pictures.npy("vs").unwrap()
    );
    let summary: serde_json::Value =
        serde_json::from_str(&fs::read_to_string(out.join("inversion.json")).unwrap()).unwrap();
    assert_eq!(summary["iterations"].as_array().unwrap().len(), 2);

    // Its pictures, and its true model back; nothing else of the run's.
    let picture = agent
        .get(&format!(
            "{u}/api/runs/{}/model?parameter=vs&show=update&plane=xz&index=9&gain=1",
            s.id
        ))
        .call()
        .unwrap();
    assert_eq!(picture.status().as_u16(), 200);
    assert!(picture.headers().contains_key("x-low"));
    let mut t = agent
        .get(&format!("{u}/api/runs/{}/files/truth.toml", s.id))
        .call()
        .unwrap();
    assert_eq!(t.body_mut().read_to_string().unwrap(), truth);
    for other in ["status.json", "config.toml", "out/vs.npy"] {
        let code = agent
            .get(&format!("{u}/api/runs/{}/files/{other}", s.id))
            .call()
            .unwrap()
            .status()
            .as_u16();
        assert_eq!(code, 404, "{other}");
    }
    let list = Command::new(BIN)
        .args(["remote", u, "list"])
        .output()
        .unwrap();
    assert!(String::from_utf8_lossy(&list.stdout).contains("sphere.toml (inversion)"));
    fs::remove_dir_all(dir).unwrap();
}
