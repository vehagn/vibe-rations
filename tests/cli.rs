//! The command-line interface, run as its own process.

use std::fs;
use std::path::{Path, PathBuf};
use std::process::{Command, Output};

fn vibe(args: &[&str]) -> Output {
    Command::new(env!("CARGO_BIN_EXE_vibe-rations"))
        .args(args)
        .current_dir(env!("CARGO_MANIFEST_DIR"))
        .output()
        .expect("the binary runs")
}

fn stderr(out: &Output) -> String {
    String::from_utf8_lossy(&out.stderr).into_owned()
}

/// A fresh, empty directory for one test.
fn scratch(name: &str) -> PathBuf {
    let dir = std::env::temp_dir().join(format!("vibe-rations-cli-{name}-{}", std::process::id()));
    let _ = fs::remove_dir_all(&dir);
    fs::create_dir_all(&dir).unwrap();
    dir
}

const TINY: &str = r#"
[grid]
points = [24, 24, 20]
spacing = 10.0
[time]
duration = 0.05
[boundary]
pml_width = 8
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

#[test]
fn check_describes_a_configuration() {
    let out = vibe(&["check", "examples/explosion.toml"]);
    let err = stderr(&out);
    assert!(out.status.success(), "{err}");
    for label in [
        "grid", "model", "time", "accuracy", "boundary", "output", "memory",
    ] {
        assert!(
            err.lines().any(|l| l.starts_with(label)),
            "no {label} line in\n{err}"
        );
    }
}

/// `memory` itemises a model's memory without setting it up, so even one far larger than this
/// machine is answered at once; `check` refuses such a model before allocating it.
#[test]
fn memory_itemises_a_model_without_building_it() {
    let out = vibe(&["memory", "examples/explosion.toml"]);
    let text = String::from_utf8_lossy(&out.stdout).into_owned();
    assert!(out.status.success(), "{}", stderr(&out));
    for part in ["wavefield", "material", "model", "traces", "at most"] {
        assert!(text.contains(part), "no {part} in\n{text}");
    }
    // A gradient run holds the adjoint's solver, checkpoints and a replayed segment too.
    let out = vibe(&["memory", "--gradient", "examples/explosion.toml"]);
    let gradient = String::from_utf8_lossy(&out.stdout).into_owned();
    assert!(out.status.success(), "{}", stderr(&out));
    for part in [
        "adjoint's solver",
        "checkpoints",
        "replayed segment",
        "strains",
    ] {
        assert!(gradient.contains(part), "no {part} in\n{gradient}");
    }
    let dir = scratch("memory");
    let huge = dir.join("huge.toml");
    fs::write(
        &huge,
        TINY.replace("points = [24, 24, 20]", "points = [4000, 4000, 2000]"),
    )
    .unwrap();
    let out = vibe(&["memory", huge.to_str().unwrap()]);
    let text = String::from_utf8_lossy(&out.stdout).into_owned();
    assert!(out.status.success(), "{}", stderr(&out));
    assert!(text.contains("GiB"), "{text}");
    if Path::new("/proc/meminfo").exists() {
        let out = vibe(&["check", huge.to_str().unwrap()]);
        assert!(!out.status.success());
        assert!(stderr(&out).contains("only"), "{}", stderr(&out));
    }
    fs::remove_dir_all(dir).unwrap();
}

#[test]
fn run_writes_traces_gathers_snapshots_and_settings() {
    let dir = scratch("run");
    let config = dir.join("tiny.toml");
    fs::write(&config, TINY).unwrap();
    let out_dir = dir.join("out");
    let out = vibe(&[
        "run",
        config.to_str().unwrap(),
        "--out",
        out_dir.to_str().unwrap(),
        "--quiet",
    ]);
    assert!(out.status.success(), "{}", stderr(&out));
    for file in [
        "vz.npy",
        "p.npy",
        "time.npy",
        "receivers.csv",
        "run.json",
        "vz_gather.png",
        "p_gather.png",
    ] {
        assert!(out_dir.join(file).is_file(), "no {file}");
    }
    let snapshots: Vec<String> = fs::read_dir(out_dir.join("snapshots"))
        .unwrap()
        .map(|e| e.unwrap().file_name().to_string_lossy().into_owned())
        .collect();
    assert!(
        snapshots.iter().any(|s| s.ends_with(".png")),
        "{snapshots:?}"
    );
    assert!(
        snapshots.iter().any(|s| s.ends_with(".npy")),
        "{snapshots:?}"
    );
    let csv = fs::read_to_string(out_dir.join("receivers.csv")).unwrap();
    assert_eq!(csv.lines().count(), 1 + 4);
    let info: serde_json::Value =
        serde_json::from_str(&fs::read_to_string(out_dir.join("run.json")).unwrap()).unwrap();
    assert_eq!(info["receivers"], 4);
    assert_eq!(info["points"], serde_json::json!([24, 24, 20]));
    fs::remove_dir_all(dir).unwrap();
}

/// The gradient of a model's misfit to traces another model recorded.
#[test]
fn gradient_writes_one_array_per_parameter() {
    let dir = scratch("gradient");
    // The same time step for both, or the samples wouldn't line up.
    let tiny = TINY.replace("duration = 0.05", "duration = 0.05\ndt = 6e-4");
    let truth = dir.join("truth.toml");
    fs::write(&truth, &tiny).unwrap();
    let observed = dir.join("observed");
    let out = vibe(&[
        "run",
        truth.to_str().unwrap(),
        "--out",
        observed.to_str().unwrap(),
        "--quiet",
    ]);
    assert!(out.status.success(), "{}", stderr(&out));
    let guess = dir.join("guess.toml");
    fs::write(&guess, tiny.replace("vp = 2000.0", "vp = 2100.0")).unwrap();
    let grad = dir.join("grad");
    let args = |observed: &Path| {
        vibe(&[
            "gradient",
            guess.to_str().unwrap(),
            "--observed",
            observed.to_str().unwrap(),
            "--out",
            grad.to_str().unwrap(),
            "--balance",
        ])
    };
    let out = args(&observed);
    assert!(out.status.success(), "{}", stderr(&out));
    assert!(
        stderr(&out).contains("memory --gradient"),
        "{}",
        stderr(&out)
    );
    for name in ["vp", "vs", "rho"] {
        let bytes = fs::read(grad.join(format!("grad_{name}.npy"))).unwrap();
        assert!(bytes.starts_with(b"\x93NUMPY"), "grad_{name}.npy");
        assert!(bytes.len() > 4 * 24 * 24 * 20, "grad_{name}.npy");
    }
    assert!(grad.join("vz.npy").is_file());
    let info: serde_json::Value =
        serde_json::from_str(&fs::read_to_string(grad.join("gradient.json")).unwrap()).unwrap();
    assert!(info["misfit"].as_f64().unwrap() > 0.0, "{info}");
    assert_eq!(info["weights"].as_array().unwrap().len(), 2);
    // TINY's receivers are clear of the layers.
    assert_eq!(info["receivers_next_to_layers"], 0, "vz is read along z");
    // The observed traces must be there, sampled as the receivers record.
    fs::write(&guess, TINY.replace("vp = 2000.0", "vp = 2100.0")).unwrap();
    let out = args(&observed);
    assert!(!out.status.success());
    assert!(stderr(&out).contains("same [time] dt"), "{}", stderr(&out));
    fs::write(&guess, tiny.replace("vp = 2000.0", "vp = 2100.0")).unwrap();
    fs::remove_file(observed.join("p.npy")).unwrap();
    let out = args(&observed);
    assert!(!out.status.success());
    assert!(stderr(&out).contains("p.npy"), "{}", stderr(&out));
    fs::remove_dir_all(dir).unwrap();
}

#[test]
fn reports_a_typo_and_fails() {
    let dir = scratch("typo");
    let config = dir.join("typo.toml");
    fs::write(&config, TINY.replace("duration", "duraton")).unwrap();
    let out = vibe(&["check", config.to_str().unwrap()]);
    let err = stderr(&out);
    assert!(!out.status.success());
    assert!(
        err.starts_with("error:") && err.contains("duraton"),
        "{err}"
    );
    fs::remove_dir_all(dir).unwrap();
}

#[test]
fn reports_a_missing_file_and_fails() {
    let out = vibe(&["check", Path::new("no/such.toml").to_str().unwrap()]);
    assert!(!out.status.success());
    assert!(
        stderr(&out).contains("reading no/such.toml"),
        "{}",
        stderr(&out)
    );
}
