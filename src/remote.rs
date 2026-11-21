//! The client of `vibe-rations serve` (crate::server): submitting runs, following them and
//! fetching what they wrote, either at a URL or on a machine reached over SSH, through a
//! tunnel this opens with the system's `ssh` (so its keys, agent and config apply).

use std::fs::{self, File};
use std::io::{self, IsTerminal, Write as _};
use std::net::{TcpListener, TcpStream};
use std::path::Path;
use std::process::{Child, Command, Stdio};
use std::time::{Duration, Instant};

use anyhow::{Context, Result, bail, ensure};
use serde::de::DeserializeOwned;
use ureq::Agent;
use ureq::http::Response;

use crate::config::Config;
use crate::memory::human;
use crate::server::{Kind, Phase, Status, TRUTH};
use crate::simulation::human_time;

/// A connection to a run server.
pub struct Remote {
    base: String,
    agent: Agent,
    tunnel: Option<Child>,
}

impl Drop for Remote {
    fn drop(&mut self) {
        if let Some(mut t) = self.tunnel.take() {
            let _ = t.kill();
            let _ = t.wait();
        }
    }
}

/// Connects to the server at `target`: a URL (`http://localhost:8089`), or an SSH destination
/// (`user@host`, or a Host from ~/.ssh/config) whose server listens on `port` on that machine.
pub fn connect(target: &str, port: u16) -> Result<Remote> {
    let agent: Agent = Agent::config_builder()
        .http_status_as_error(false)
        .build()
        .into();
    if target.starts_with("http://") || target.starts_with("https://") {
        let r = Remote {
            base: target.trim_end_matches('/').to_string(),
            agent,
            tunnel: None,
        };
        r.info().with_context(|| format!("reaching {target}"))?;
        return Ok(r);
    }
    // A free port here, forwarded to the server's port there.
    let local = TcpListener::bind("127.0.0.1:0")?.local_addr()?.port();
    let mut child = Command::new("ssh")
        .args([
            "-N",
            "-o",
            "ExitOnForwardFailure=yes",
            "-L",
            &format!("{local}:127.0.0.1:{port}"),
            target,
        ])
        .stdout(Stdio::null())
        .stderr(Stdio::piped())
        .spawn()
        .context("starting ssh; is it installed?")?;
    let start = Instant::now();
    loop {
        if TcpStream::connect(("127.0.0.1", local)).is_ok() {
            break;
        }
        if let Some(exit) = child.try_wait()? {
            let mut stderr = String::new();
            if let Some(mut e) = child.stderr.take() {
                let _ = io::Read::read_to_string(&mut e, &mut stderr);
            }
            bail!("ssh to {target} ended ({exit}): {}", stderr.trim());
        }
        if start.elapsed() > Duration::from_secs(60) {
            let _ = child.kill();
            bail!("ssh to {target} didn't open its tunnel within a minute");
        }
        std::thread::sleep(Duration::from_millis(100));
    }
    let r = Remote {
        base: format!("http://127.0.0.1:{local}"),
        agent,
        tunnel: Some(child),
    };
    r.info().with_context(|| {
        format!("reaching vibe-rations serve on {target}, port {port}; is it running there?")
    })?;
    Ok(r)
}

impl Remote {
    fn url(&self, path: &str) -> String {
        format!("{}{path}", self.base)
    }

    /// The response, or the server's error message as an error.
    fn check(
        &self,
        response: Result<Response<ureq::Body>, ureq::Error>,
    ) -> Result<Response<ureq::Body>> {
        let mut response = response.context("talking to the server")?;
        let code = response.status().as_u16();
        if code >= 400 {
            let text = response.body_mut().read_to_string().unwrap_or_default();
            let message = serde_json::from_str::<serde_json::Value>(&text)
                .ok()
                .and_then(|v| v["error"].as_str().map(str::to_string))
                .unwrap_or(text);
            bail!("{message} (HTTP {code})");
        }
        Ok(response)
    }

    fn parse<T: DeserializeOwned>(&self, response: Response<ureq::Body>) -> Result<T> {
        let mut response = response;
        let text = response
            .body_mut()
            .with_config()
            .limit(64 << 20)
            .read_to_string()?;
        serde_json::from_str(&text).context("reading the server's answer")
    }

    fn get<T: DeserializeOwned>(&self, path: &str) -> Result<T> {
        let r = self.check(self.agent.get(&self.url(path)).call())?;
        self.parse(r)
    }

    fn post<T: DeserializeOwned>(&self, path: &str, body: &str) -> Result<T> {
        let r = self.check(self.agent.post(&self.url(path)).send(body))?;
        self.parse(r)
    }

    /// The server: version, threads, memory available.
    pub fn info(&self) -> Result<serde_json::Value> {
        self.get("/api/info")
    }

    pub fn runs(&self) -> Result<Vec<Status>> {
        self.get("/api/runs")
    }

    pub fn status(&self, id: &str) -> Result<Status> {
        self.get(&format!("/api/runs/{}", encode(id)))
    }

    /// Sends the configuration at `path` and the model files it names, then starts the run.
    pub fn submit(&self, path: &Path) -> Result<Status> {
        let text =
            fs::read_to_string(path).with_context(|| format!("reading {}", path.display()))?;
        let cfg = Config::from_toml(&text)?;
        let name = path
            .file_name()
            .map_or("run.toml".into(), |n| n.to_string_lossy());
        let status: Status = self.post(&format!("/api/runs?name={}", encode(&name)), &text)?;
        let base = path.parent().unwrap_or(Path::new("."));
        if let Some(files) = &cfg.model.files {
            for f in [&files.vp, &files.vs, &files.rho] {
                let local = base.join(f);
                let file =
                    File::open(&local).with_context(|| format!("opening {}", local.display()))?;
                let size = file.metadata()?.len();
                eprintln!("uploading {} ({})", f.display(), human(size));
                let url = self.url(&format!(
                    "/api/runs/{}/files/{}",
                    encode(&status.id),
                    encode(&f.to_string_lossy())
                ));
                self.check(self.agent.put(&url).send(file))?;
            }
        }
        self.post(&format!("/api/runs/{}/start", encode(&status.id)), "")
    }

    /// Starts an inversion of synthetic data named `name`: `start` is the configuration it
    /// starts from, with `[inversion]`, and `truth` the same survey with the true model, both
    /// sampled alike (as the web page makes them).
    pub fn submit_inversion(&self, name: &str, start: &str, truth: &str) -> Result<Status> {
        let status: Status = self.post(
            &format!("/api/runs?name={}&kind=inversion", encode(name)),
            start,
        )?;
        ensure!(
            status.needs == [TRUTH],
            "the inversion needs {:?} besides its true model",
            status.needs
        );
        let url = self.url(&format!("/api/runs/{}/files/{TRUTH}", encode(&status.id)));
        self.check(self.agent.put(&url).send(truth))?;
        self.post(&format!("/api/runs/{}/start", encode(&status.id)), "")
    }

    /// Stops a run and deletes it with its files.
    pub fn delete(&self, id: &str) -> Result<()> {
        self.check(
            self.agent
                .delete(&self.url(&format!("/api/runs/{}", encode(id))))
                .call(),
        )?;
        Ok(())
    }

    /// Waits until the run has finished, reporting its progress on stderr if `progress`: one
    /// line, rewritten in place on a terminal, otherwise a line every minute.
    pub fn follow(&self, id: &str, progress: bool) -> Result<Status> {
        let terminal = io::stderr().is_terminal();
        let every = if terminal {
            Duration::from_secs(2)
        } else {
            Duration::from_secs(60)
        };
        let mut last: Option<Instant> = None;
        let mut paused = false;
        let mut drawn = false;
        loop {
            let s = self.status(id)?;
            if s.phase.finished() {
                if progress && drawn {
                    eprintln!("{}", if terminal { "" } else { "finished" });
                }
                return Ok(s);
            }
            if s.phase == Phase::Paused && !paused && progress {
                eprintln!("{}paused at step {}", if drawn { "\n" } else { "" }, s.step);
                drawn = false;
            }
            paused = s.phase == Phase::Paused;
            if progress && s.phase == Phase::Running && last.is_none_or(|t| t.elapsed() >= every) {
                last = Some(Instant::now());
                let line = progress_line(&s);
                if terminal {
                    eprint!("\r{line}\x1b[K");
                } else {
                    eprintln!("{line}");
                }
                drawn = true;
            }
            std::thread::sleep(Duration::from_secs(2));
        }
    }

    /// The run's results now, if it has finished: into `out`, as [`Remote::download`].
    pub fn fetch(&self, id: &str, out: &Path) -> Result<usize> {
        let s = self.status(id)?;
        ensure!(
            s.phase == Phase::Done,
            "{id} is {}{}",
            format!("{:?}", s.phase).to_lowercase(),
            match s.phase {
                Phase::Running | Phase::Paused => format!(", {}", progress_line(&s)),
                _ => s.error.map(|e| format!(": {e}")).unwrap_or_default(),
            }
        );
        self.download(id, out)
    }

    /// Fetches everything the run wrote into `out`, keeping its layout; returns the count.
    pub fn download(&self, id: &str, out: &Path) -> Result<usize> {
        #[derive(serde::Deserialize)]
        struct Output {
            path: String,
        }
        let files: Vec<Output> = self.get(&format!("/api/runs/{}/outputs", encode(id)))?;
        for f in &files {
            let target = out.join(&f.path);
            ensure!(
                Path::new(&f.path)
                    .components()
                    .all(|c| matches!(c, std::path::Component::Normal(_))),
                "the server named an output outside the run: {}",
                f.path
            );
            if let Some(parent) = target.parent() {
                fs::create_dir_all(parent)?;
            }
            let url = self.url(&format!(
                "/api/runs/{}/outputs/{}",
                encode(id),
                encode(&f.path)
            ));
            let mut response = self.check(self.agent.get(&url).call())?;
            let mut file =
                File::create(&target).with_context(|| format!("creating {}", target.display()))?;
            io::copy(
                &mut response.body_mut().with_config().limit(u64::MAX).reader(),
                &mut file,
            )?;
            file.flush()?;
        }
        Ok(files.len())
    }
}

/// A run's progress in one line: share done, step, speed and time left; for an inversion, its
/// units of work (shots, the start, iterations) and its misfit.
pub fn progress_line(s: &Status) -> String {
    let share = if s.total_steps > 0 {
        100.0 * s.step as f64 / s.total_steps as f64
    } else {
        0.0
    };
    if s.kind == Kind::Inversion {
        let mut line = format!("{share:5.1}%  {}/{} units", s.step, s.total_steps);
        let i = s.inversion.as_ref();
        if let Some(stage) = i.and_then(|i| i["stage"].as_str()) {
            line += &format!("  {stage}");
        }
        let misfit = |k: &str| i.and_then(|i| i[k].as_f64());
        if let (Some(now), Some(first)) = (misfit("misfit"), misfit("starting_misfit")) {
            line += &format!("  misfit {:.1}% of the start", 100.0 * now / first);
        }
        return line;
    }
    let mut line = format!("{share:5.1}%  step {}/{}", s.step, s.total_steps);
    if let Some(left) = seconds_left(s) {
        line += &format!("  {:.0} Mcell/s  {} left", s.speed, human_time(left));
    }
    line
}

/// Seconds until the run is done at its speed so far, once it has one.
pub fn seconds_left(s: &Status) -> Option<f64> {
    let cells: f64 = s.info.as_ref()?["points"]
        .as_array()?
        .iter()
        .filter_map(|v| v.as_f64())
        .product();
    (s.speed > 0.0 && s.step > 0)
        .then(|| cells * s.total_steps.saturating_sub(s.step) as f64 / (s.speed * 1e6))
}

/// Percent-encoding for a URL path or query value; keeps / so paths stay paths.
fn encode(s: &str) -> String {
    let mut out = String::with_capacity(s.len());
    for b in s.bytes() {
        if b.is_ascii_alphanumeric() || b"-_.~/".contains(&b) {
            out.push(b as char);
        } else {
            out.push_str(&format!("%{b:02X}"));
        }
    }
    out
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn encodes_what_a_url_path_cannot_hold() {
        assert_eq!(encode("models/vp 1.npy"), "models/vp%201.npy");
        assert_eq!(encode("a&b=c"), "a%26b%3Dc");
    }
}
