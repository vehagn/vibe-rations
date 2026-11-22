//! Every configuration in `examples/` and `examples/marmousi2/` parses and sets up, without
//! warnings. Marmousi2's model files are made by `just marmousi2`, not kept in the repository:
//! until they are, the configurations that need them are only parsed.

use std::fs;
use std::path::Path;

use vibe_rations::config::Config;
use vibe_rations::simulation::Simulation;

/// Checks every configuration in `dir`; returns how many there were and how many set up.
fn check(dir: &Path) -> (usize, usize) {
    let (mut seen, mut set_up) = (0, 0);
    for entry in fs::read_dir(dir).unwrap() {
        let path = entry.unwrap().path();
        if path.extension().is_none_or(|e| e != "toml") {
            continue;
        }
        let cfg = Config::load(&path).unwrap_or_else(|e| panic!("{}: {e:#}", path.display()));
        seen += 1;
        if let Some(files) = &cfg.model.files
            && [&files.vp, &files.vs, &files.rho]
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
            Simulation::new(&cfg, dir).unwrap_or_else(|e| panic!("{}: {e:#}", path.display()));
        assert!(
            sim.info.warnings.is_empty(),
            "{}: {:?}",
            path.display(),
            sim.info.warnings
        );
        set_up += 1;
    }
    (seen, set_up)
}

#[test]
fn every_example_sets_up_cleanly() {
    let dir = Path::new(env!("CARGO_MANIFEST_DIR")).join("examples");
    let (seen, set_up) = check(&dir);
    // explosion and layered need no files, the two Marmousi2 ones do.
    assert!(seen >= 4 && set_up >= 2, "{seen} examples, {set_up} set up");
}

#[test]
fn the_marmousi2_examples_parse_and_set_up_once_made() {
    let dir = Path::new(env!("CARGO_MANIFEST_DIR")).join("examples/marmousi2");
    assert_eq!(check(&dir).0, 5);
}
