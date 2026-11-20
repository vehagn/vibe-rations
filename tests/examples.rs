//! Every configuration in `examples/` parses and sets up, without warnings.

use std::fs;
use std::path::Path;

use vibe_rations::config::Config;
use vibe_rations::simulation::Simulation;

#[test]
fn every_example_sets_up_cleanly() {
    let dir = Path::new(env!("CARGO_MANIFEST_DIR")).join("examples");
    let mut seen = 0;
    for entry in fs::read_dir(&dir).unwrap() {
        let path = entry.unwrap().path();
        if path.extension().is_none_or(|e| e != "toml") {
            continue;
        }
        let cfg = Config::load(&path).unwrap_or_else(|e| panic!("{}: {e:#}", path.display()));
        let sim =
            Simulation::new(&cfg, &dir).unwrap_or_else(|e| panic!("{}: {e:#}", path.display()));
        assert!(
            sim.info.warnings.is_empty(),
            "{}: {:?}",
            path.display(),
            sim.info.warnings
        );
        seen += 1;
    }
    assert!(seen >= 2, "found only {seen} examples");
}

/// The Marmousi2 configurations, whose model files `just marmousi2` downloads and makes,
/// parse; where the files are there, they set up cleanly too.
#[test]
fn the_marmousi2_examples_parse_and_set_up_once_made() {
    let dir = Path::new(env!("CARGO_MANIFEST_DIR")).join("examples/marmousi2");
    let mut seen = 0;
    for entry in fs::read_dir(&dir).unwrap() {
        let path = entry.unwrap().path();
        if path.extension().is_none_or(|e| e != "toml") {
            continue;
        }
        let cfg = Config::load(&path).unwrap_or_else(|e| panic!("{}: {e:#}", path.display()));
        seen += 1;
        let files = cfg
            .model
            .files
            .as_ref()
            .expect("Marmousi2 comes as model files");
        if [&files.vp, &files.vs, &files.rho]
            .iter()
            .any(|f| !dir.join(f).exists())
        {
            eprintln!(
                "{}: no model files; `just marmousi2` makes them",
                path.display()
            );
            continue;
        }
        let sim =
            Simulation::new(&cfg, &dir).unwrap_or_else(|e| panic!("{}: {e:#}", path.display()));
        assert!(
            sim.info.warnings.is_empty(),
            "{}: {:?}",
            path.display(),
            sim.info.warnings
        );
    }
    assert_eq!(seen, 3);
}
