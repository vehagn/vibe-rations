//! The client of `vibe-rations serve` (crate::server): submitting runs, following them and
//! fetching what they wrote, either at a URL or on a machine reached over SSH, through a
//! tunnel this opens with the system's `ssh` (so its keys, agent and config apply).

use std::fs::{self, File};
use std::io::{self, Write as _};
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
use crate::server::{Phase, Status};

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

    /// Stops a run and deletes it with its files.
    pub fn delete(&self, id: &str) -> Result<()> {
        self.check(
            self.agent
                .delete(&self.url(&format!("/api/runs/{}", encode(id))))
                .call(),
        )?;
        Ok(())
    }

    /// Waits until the run has finished, reporting progress on stderr if `progress`.
    pub fn follow(&self, id: &str, progress: bool) -> Result<Status> {
        let mut last = Instant::now() - Duration::from_secs(10);
        let mut paused = false;
        loop {
            let s = self.status(id)?;
            if s.phase.finished() {
                return Ok(s);
            }
            if s.phase == Phase::Paused && !paused && progress {
                eprintln!("paused at step {}; resume it from the web page", s.step);
            }
            paused = s.phase == Phase::Paused;
            if progress && s.phase == Phase::Running && last.elapsed() >= Duration::from_secs(2) {
                last = Instant::now();
                let left = if s.speed > 0.0 && s.step > 0 {
                    let cells = s.info.as_ref().map_or(0.0, |i| {
                        i["points"]
                            .as_array()
                            .map_or(0.0, |p| p.iter().filter_map(|v| v.as_f64()).product())
                    });
                    format!(
                        "  {:.0}s left",
                        cells * (s.total_steps - s.step) as f64 / (s.speed * 1e6)
                    )
                } else {
                    String::new()
                };
                eprintln!(
                    "step {:>6}/{}  t {:.4} s  {:.0} Mcells/s{left}",
                    s.step, s.total_steps, s.time, s.speed
                );
            }
            std::thread::sleep(Duration::from_secs(1));
        }
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
