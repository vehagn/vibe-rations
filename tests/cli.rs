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

/// With --shots each source fires on its own, into its own directory.
#[test]
fn run_fires_shots_one_at_a_time() {
    let dir = scratch("shots");
    let config = dir.join("two.toml");
    let two = TINY.replace(
        "[receivers]",
        "[[source]]\nposition = [60.0, 120.0, 60.0]\nforce = [0.0, 0.0, 1e9]\n\
         wavelet = { kind = \"ricker\", frequency = 10.0 }\n[receivers]",
    );
    fs::write(&config, two).unwrap();
    let out_dir = dir.join("out");
    let out = vibe(&[
        "run",
        config.to_str().unwrap(),
        "--out",
        out_dir.to_str().unwrap(),
        "--quiet",
        "--shots",
    ]);
    assert!(out.status.success(), "{}", stderr(&out));
    for (shot, at) in [
        ("shot-000", [120.0, 120.0, 60.0]),
        ("shot-001", [60.0, 120.0, 60.0]),
    ] {
        assert!(out_dir.join(shot).join("vz.npy").is_file(), "{shot}");
        let info: serde_json::Value =
            serde_json::from_str(&fs::read_to_string(out_dir.join(shot).join("run.json")).unwrap())
                .unwrap();
        let sources = info["sources"].as_array().unwrap();
        assert_eq!(sources.len(), 1, "{shot}");
        assert_eq!(sources[0]["position"], serde_json::json!(at), "{shot}");
    }
    assert!(!out_dir.join("shot-002").exists());
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

/// An inversion of two shots another model recorded, one iteration: the model and the
/// history come out, and a configuration without [inversion] is refused.
#[test]
fn invert_writes_the_model_and_its_history() {
    let dir = scratch("invert");
    let survey = TINY
        .replace("duration = 0.05", "duration = 0.05\ndt = 6e-4")
        .replace("pml_width = 8", "pml_width = 8\npml_vp_max = 2500.0")
        .replace(
            "[receivers]",
            "[[source]]\nposition = [150.0, 110.0, 60.0]\nforce = [0.0, 0.0, 1e9]\n\
             wavelet = { kind = \"ricker\", frequency = 10.0 }\n[receivers]",
        );
    let truth = dir.join("truth.toml");
    fs::write(&truth, &survey).unwrap();
    let observed = dir.join("observed");
    let out = vibe(&[
        "run",
        truth.to_str().unwrap(),
        "--out",
        observed.to_str().unwrap(),
        "--quiet",
        "--shots",
    ]);
    assert!(out.status.success(), "{}", stderr(&out));
    let guess = dir.join("guess.toml");
    let model_out = dir.join("inverted");
    fs::write(
        &guess,
        survey.replace("vp = 2000.0", "vp = 2080.0")
            + "\n[inversion]\nobserved = \"observed\"\niterations = 1\n",
    )
    .unwrap();
    let out = vibe(&[
        "invert",
        guess.to_str().unwrap(),
        "--out",
        model_out.to_str().unwrap(),
    ]);
    let err = stderr(&out);
    assert!(out.status.success(), "{err}");
    assert!(err.contains("iteration   1"), "{err}");
    // Shaped like model files, (nz, ny, nx), and changed from the start somewhere.
    let mut moved = false;
    for (name, start) in [("vp", 2080.0), ("vs", 1000.0), ("rho", 2000.0)] {
        let (shape, values) =
            vibe_rations::npy::read(&model_out.join(format!("{name}.npy"))).unwrap();
        assert_eq!(shape, vec![20, 24, 24], "{name}.npy");
        moved |= values.iter().any(|&v| v != start);
    }
    assert!(moved, "the inverted model is the starting one");
    let info: serde_json::Value =
        serde_json::from_str(&fs::read_to_string(model_out.join("inversion.json")).unwrap())
            .unwrap();
    let iterations = info["iterations"].as_array().unwrap();
    assert_eq!(iterations.len(), 1, "{info}");
    let (start, now) = (
        info["starting_misfit"].as_f64().unwrap(),
        iterations[0]["misfit"].as_f64().unwrap(),
    );
    assert!(now < start, "{now} vs {start}");
    // Without [inversion] there is nothing to do.
    let out = vibe(&["invert", truth.to_str().unwrap()]);
    assert!(!out.status.success());
    assert!(stderr(&out).contains("[inversion]"), "{}", stderr(&out));
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
