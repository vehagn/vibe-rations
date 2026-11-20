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
