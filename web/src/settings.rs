//! The page's settings form and its commands. The configuration's text stays the source of
//! truth: [`settings`] reads what it sets, defaults filled in, for the form to show, and
//! [`edit`] writes the form's changes back into the text, keeping its comments and layout.
//! [`commands`] are the command lines that run the configuration without the page, and
//! [`synthetic`] makes the pair of configurations a synthetic inversion needs.

use anyhow::{Context, Result, bail, ensure};
use serde_json::{Value as Json, json};
use std::path::Path;
use toml_edit::{Array, ArrayOfTables, DocumentMut, InlineTable, Item, Table, TableLike, Value};
use vibe_rations::config::{
    Coefficients, Config, InversionConfig, Spacing, WaveletConfig, default_courant,
};
use vibe_rations::inversion::shared_operator;

/// What a configuration sets, with every default filled in, as JSON for the form.
pub fn settings(text: &str) -> Result<Json> {
    let c = Config::from_toml(text)?;
    let spacing = match c.grid.spacing {
        Spacing::Uniform(h) => [h; 3],
        Spacing::PerAxis(h) => h,
    };
    let sources: Vec<Json> = c
        .sources
        .iter()
        .map(|s| {
            let (kind, frequency, delay) = match s.wavelet {
                WaveletConfig::Ricker { frequency, delay } => ("ricker", frequency, delay),
                WaveletConfig::Gaussian { frequency, delay } => ("gaussian", frequency, delay),
            };
            let mechanism = if s.explosion.is_some() {
                "explosion"
            } else if s.moment.is_some() {
                "moment"
            } else if s.double_couple.is_some() {
                "double couple"
            } else if s.force.is_some() {
                "force"
            } else {
                "none"
            };
            json!({
                "position": s.position,
                "wavelet": kind,
                "frequency": frequency,
                "delay": delay,
                "mechanism": mechanism,
            })
        })
        .collect();
    let receivers = c.receivers.positions()?.len();
    let inversion = match &c.inversion {
        Some(i) => inversion_json(i),
        None => Json::Null,
    };
    Ok(json!({
        "grid": { "points": c.grid.points, "spacing": spacing },
        "time": {
            "duration": c.time.duration,
            "dt": c.time.dt,
            "courant": c.time.courant,
            "default_courant": default_courant(c.scheme.order, c.scheme.remove_time_dispersion),
        },
        "scheme": {
            "order": c.scheme.order,
            "coefficients": match c.scheme.coefficients {
                Coefficients::Taylor => "taylor",
                Coefficients::Optimised => "optimised",
            },
            "remove_time_dispersion": c.scheme.remove_time_dispersion,
            "optimised_for": c.scheme.optimised_for,
        },
        "boundary": {
            "pml_width": c.boundary.pml_width,
            "free_surface": c.boundary.free_surface,
            "pml_sponge": c.boundary.pml_sponge,
            "pml_vp_max": c.boundary.pml_vp_max,
        },
        "output": { "record_every": c.output.record_every, "dir": c.output.dir },
        "sources": sources,
        "receivers": {
            "components": c.receivers.components.iter().map(|c| c.name()).collect::<Vec<_>>(),
            "count": receivers,
            "lines": c.receivers.lines.len(),
            "points": c.receivers.points.len(),
        },
        "model": {
            "files": c.model.files.is_some(),
            "layers": c.model.layers.iter().map(|l| json!({
                "top": l.top, "vp": l.vp, "vs": l.vs, "rho": l.rho,
            })).collect::<Vec<_>>(),
            "spheres": c.model.spheres.iter().map(|s| json!({
                "center": s.center, "radius": s.radius, "vp": s.vp, "vs": s.vs, "rho": s.rho,
            })).collect::<Vec<_>>(),
        },
        "inversion": inversion,
        "inversion_defaults": inversion_json(&default_inversion()),
    }))
}

/// `[inversion]` as `invert` reads it when it says nothing but where the data are.
fn default_inversion() -> InversionConfig {
    toml::from_str("observed = \"observed\"").expect("the defaults parse")
}

fn inversion_json(i: &InversionConfig) -> Json {
    serde_json::to_value(i).expect("InversionConfig serialises")
}

/// Keys whose values are whole numbers; every other number is written as a float, so that
/// `frequency = 12` from the form comes out as `12.0`.
const INTEGER_KEYS: [&str; 7] = [
    "points",
    "order",
    "pml_width",
    "record_every",
    "count",
    "iterations",
    "memory",
];

/// The configuration `text` with `edits` applied in turn: a JSON array of `{path, value}`.
///
/// A path is dotted keys, with numbers indexing arrays of tables: `scheme.order`,
/// `source.0.wavelet.frequency`, `model.layer.2.vp`. A `null` value removes the key (or that
/// table of the array, and the array with its last table); a path one past the end of an
/// array of tables appends the table its value describes. Tables and keys that don't exist yet
/// are added; everything else in the text is left as it was, comments included.
pub fn edit(text: &str, edits: &Json) -> Result<String> {
    let mut doc: DocumentMut = text.parse().context("the configuration isn't valid TOML")?;
    let edits = edits.as_array().context("edits must be a list")?;
    for e in edits {
        let path = e["path"].as_str().context("an edit needs a path")?;
        match e["copy"].as_u64() {
            Some(from) => copy(doc.as_table_mut(), path, from as usize),
            None => set(doc.as_table_mut(), path, &e["value"]),
        }
        .with_context(|| format!("setting {path}"))?;
    }
    Ok(doc.to_string())
}

/// Appends a copy of entry `from` of the array of tables at `path` (a key of the root), its
/// comments and formatting included.
fn copy(root: &mut dyn TableLike, path: &str, from: usize) -> Result<()> {
    let array = root
        .get_mut(path)
        .and_then(Item::as_array_of_tables_mut)
        .with_context(|| format!("{path} isn't an array of tables"))?;
    let mut table = array
        .get(from)
        .with_context(|| format!("{path} has no entry {}", from + 1))?
        .clone();
    table.set_position(None);
    array.push(table);
    Ok(())
}

/// Sets `path` under `table` to `v` (or removes it, for null).
fn set(table: &mut dyn TableLike, path: &str, v: &Json) -> Result<()> {
    let segments: Vec<&str> = path.split('.').collect();
    ensure!(segments.iter().all(|s| !s.is_empty()), "an empty key");
    let create = !v.is_null();
    let mut t: &mut dyn TableLike = table;
    // Whether `t` is an inline table, whose children must be values too.
    let mut inline = false;
    let mut i = 0;
    // Down to the table holding the last key (or the array holding the last index).
    while i + 1 < segments.len() {
        let key = segments[i];
        ensure!(
            key.parse::<usize>().is_err(),
            "{key} is an index, not a key"
        );
        if let Some(index) = segments.get(i + 1).and_then(|s| s.parse::<usize>().ok()) {
            if i + 2 == segments.len() {
                return set_element(t, key, index, v);
            }
            let Some(item) = t.get_mut(key) else {
                bail!("there is no {key}");
            };
            let array = item
                .as_array_of_tables_mut()
                .with_context(|| format!("{key} isn't an array of tables"))?;
            let len = array.len();
            t = array
                .get_mut(index)
                .with_context(|| format!("{key} has {len} entries, not {}", index + 1))?;
            inline = false;
            i += 2;
        } else {
            if t.get(key).is_none() {
                if !create {
                    return Ok(());
                }
                // Inside an inline table, an inline table; else a table with a header once it
                // has values of its own.
                let new = if inline {
                    Item::Value(Value::InlineTable(InlineTable::new()))
                } else {
                    let mut table = Table::new();
                    table.set_implicit(true);
                    Item::Table(table)
                };
                t.insert(key, new);
            }
            let item = t.get_mut(key).expect("just made");
            inline = item.is_inline_table();
            t = item
                .as_table_like_mut()
                .with_context(|| format!("{key} isn't a table"))?;
            i += 1;
        }
    }
    let key = segments[i];
    if v.is_null() {
        t.remove(key);
        return Ok(());
    }
    let mut new = to_value(key, v)?;
    // Keep the comment after the old value, and its spacing.
    if let Some(old) = t.get(key).and_then(Item::as_value) {
        *new.decor_mut() = old.decor().clone();
    }
    t.insert(key, Item::Value(new));
    Ok(())
}

/// Sets entry `index` of the array of tables `key` in `t`: removes it for null, appends it if
/// it is one past the end.
fn set_element(t: &mut dyn TableLike, key: &str, index: usize, v: &Json) -> Result<()> {
    if v.is_null() {
        if let Some(array) = t.get_mut(key).and_then(Item::as_array_of_tables_mut) {
            ensure!(index < array.len(), "{key} has no entry {}", index + 1);
            array.remove(index);
            if array.is_empty() {
                t.remove(key);
            }
        }
        return Ok(());
    }
    let object = v
        .as_object()
        .with_context(|| format!("an entry of {key} must be a table"))?;
    let mut table = Table::new();
    for (k, x) in ordered(object) {
        table.insert(k, Item::Value(to_value(k, x)?));
    }
    if t.get(key).is_none() {
        t.insert(key, Item::ArrayOfTables(ArrayOfTables::new()));
    }
    let array = t
        .get_mut(key)
        .and_then(Item::as_array_of_tables_mut)
        .with_context(|| format!("{key} isn't an array of tables"))?;
    ensure!(
        index == array.len(),
        "{key} has {} entries; only entry {} can be added",
        array.len(),
        array.len() + 1
    );
    array.push(table);
    Ok(())
}

/// The order keys are written in, as the examples have them; others follow, alphabetically.
const KEY_ORDER: [&str; 16] = [
    "top",
    "center",
    "radius",
    "position",
    "vp",
    "vs",
    "rho",
    "explosion",
    "moment",
    "double_couple",
    "force",
    "wavelet",
    "kind",
    "frequency",
    "delay",
    "observed",
];

/// An object's entries in [`KEY_ORDER`] (serde_json's own order is alphabetical).
fn ordered(map: &serde_json::Map<String, Json>) -> Vec<(&String, &Json)> {
    let mut entries: Vec<_> = map.iter().collect();
    let rank = |k: &str| {
        KEY_ORDER
            .iter()
            .position(|o| *o == k)
            .unwrap_or(KEY_ORDER.len())
    };
    entries.sort_by(|a, b| rank(a.0).cmp(&rank(b.0)).then(a.0.cmp(b.0)));
    entries
}

/// A JSON value as a TOML one, numbers typed by the key they are for.
fn to_value(key: &str, v: &Json) -> Result<Value> {
    let integer = INTEGER_KEYS.contains(&key);
    Ok(match v {
        Json::Bool(b) => Value::from(*b),
        Json::String(s) => Value::from(s.as_str()),
        Json::Number(n) => {
            let x = n.as_f64().context("a number out of range")?;
            if integer {
                ensure!(
                    x >= 0.0 && x.fract() == 0.0,
                    "{key} must be a whole number, not {x}"
                );
                Value::from(x as i64)
            } else {
                Value::from(x)
            }
        }
        Json::Array(items) => {
            let mut a = Array::new();
            for x in items {
                a.push(to_value(key, x)?);
            }
            a.fmt();
            Value::Array(a)
        }
        Json::Object(map) => {
            let mut t = InlineTable::new();
            for (k, x) in ordered(map) {
                t.insert(k, to_value(k, x)?);
            }
            t.fmt();
            Value::InlineTable(t)
        }
        Json::Null => bail!("null can't be written"),
    })
}

/// The command lines that run a configuration saved as `file` without the page, as a list of
/// `{what, command}`. `input` is JSON: `file`, `target` (an SSH destination or a server's URL,
/// for remote runs), `run` (the id of the run on the page's server, if any), `sources`, and
/// for an inversion `start` and `truth` (the files [`synthetic`] makes), `observed` and
/// `balance`.
pub fn commands(input: &Json) -> Result<Json> {
    let file = input["file"]
        .as_str()
        .filter(|f| !f.is_empty())
        .unwrap_or("run.toml");
    let target = input["target"]
        .as_str()
        .filter(|t| !t.is_empty())
        .unwrap_or("user@host");
    let (f, t) = (quote(file), quote(target));
    let sources = input["sources"].as_u64().unwrap_or(1);
    let mut out = vec![
        json!({ "what": "Run it here", "command": format!("vibe-rations run {f}") }),
        json!({ "what": "Set it up only: grid, time step, accuracy, memory",
                "command": format!("vibe-rations check {f}") }),
        json!({ "what": "In a checkout of the repository (the dev container)",
                "command": format!("just run {f}") }),
        json!({ "what": "On another machine, headless; the run carries on there",
                "command": format!("vibe-rations run --remote {t} --detach {f}") }),
        json!({ "what": "Its progress there", "command": format!("vibe-rations remote {t} list") }),
    ];
    if sources > 1 {
        out.push(json!({
            "what": "Each source a shot, as an inversion takes them",
            "command": format!("vibe-rations run --shots {f}"),
        }));
    }
    if let Some(id) = input["run"].as_str().filter(|r| !r.is_empty()) {
        let id = quote(id);
        out.push(
            json!({ "what": "Follow this run on the server, then fetch its results",
                         "command": format!("vibe-rations remote {t} attach {id}") }),
        );
        out.push(json!({ "what": "Fetch its results now",
                         "command": format!("vibe-rations remote {t} fetch {id}") }));
    }
    if let (Some(start), Some(truth)) = (input["start"].as_str(), input["truth"].as_str()) {
        let (s, tr) = (quote(start), quote(truth));
        out.push(
            json!({ "what": "Inversion: the observed traces, from the true model",
                         "command": format!("vibe-rations run --shots {tr}") }),
        );
        out.push(json!({ "what": "Inversion: its memory",
                         "command": format!("vibe-rations memory --invert {s}") }));
        out.push(json!({ "what": "Inversion: L-BFGS as [inversion] says",
                         "command": format!("vibe-rations invert {s}") }));
        if sources == 1 {
            let observed = input["observed"].as_str().unwrap_or("observed");
            let balance = if input["balance"].as_bool().unwrap_or(true) {
                " --balance"
            } else {
                ""
            };
            out.push(json!({
                "what": "Inversion: the gradient at the start only",
                "command": format!(
                    "vibe-rations gradient {s} --observed {}{balance}",
                    quote(&format!("{observed}/shot-000"))
                ),
            }));
        }
    }
    Ok(Json::Array(out))
}

/// A word for a POSIX shell: as it is if it is safe, else in single quotes.
fn quote(word: &str) -> String {
    let safe = !word.is_empty()
        && word
            .chars()
            .all(|c| c.is_ascii_alphanumeric() || "._-/@:+=,%".contains(c));
    if safe {
        word.to_string()
    } else {
        format!("'{}'", word.replace('\'', r"'\''"))
    }
}

/// A true model for a synthetic inversion of `text`: its own `[model]`, as TOML, with a sphere
/// added in the part the inversion may change, below the deepest source if there is room
/// (else in the middle), a tenth faster than what is there.
pub fn default_truth(text: &str) -> Result<String> {
    let c = Config::from_toml(text)?;
    ensure!(
        c.model.files.is_none(),
        "a true model is written as layers and spheres, not files"
    );
    let mut truth: DocumentMut = model_of(text)?.parse()?;
    let g = c.grid.grid()?;
    let w = c.boundary.pml_width as f64;
    let e = g.extent();
    let fixed = c
        .inversion
        .as_ref()
        .and_then(|i| i.fixed_above)
        .unwrap_or(0.0);
    let top = if c.boundary.free_surface {
        0.0
    } else {
        w * g.h[2]
    };
    let top = top.max(fixed);
    let bottom = e[2] - w * g.h[2];
    let x = [w * g.h[0], e[0] - w * g.h[0]];
    let y = [w * g.h[1], e[1] - w * g.h[1]];
    // Below the deepest source, if there is room there, so it doesn't sit on one.
    let h = g.h.iter().copied().fold(0.0, f64::max);
    let deepest = c
        .sources
        .iter()
        .map(|s| s.position[2])
        .fold(f64::NEG_INFINITY, f64::max);
    let below = top.max(deepest + 3.0 * h);
    let top = if bottom - below >= 6.0 * h {
        below
    } else {
        top
    };
    let center = [
        round((x[0] + x[1]) / 2.0),
        round((y[0] + y[1]) / 2.0),
        round((top.max(0.0) + bottom) / 2.0),
    ];
    let radius = round(
        [
            (x[1] - x[0]) / 6.0,
            (y[1] - y[0]) / 6.0,
            (bottom - top) / 4.0,
        ]
        .into_iter()
        .fold(f64::INFINITY, f64::min)
        .max(2.0 * h),
    );
    let layer = c
        .model
        .layers
        .iter()
        .rev()
        .find(|l| l.top <= center[2])
        .or(c.model.layers.first())
        .context("the model has no layers")?;
    let sphere = json!({
        "center": center,
        "radius": radius,
        "vp": round(1.1 * layer.vp),
        "vs": round(1.1 * layer.vs),
        "rho": layer.rho,
    });
    let n = c.model.spheres.len();
    set(truth.as_table_mut(), &format!("model.sphere.{n}"), &sphere)?;
    let text = truth.to_string();
    Ok(format!(
        "# The true model: the observed traces come from it. The configuration's own model is\n\
         # where the inversion starts.\n{}",
        text.trim_start()
    ))
}

/// Lets tables follow whatever they are put after, instead of where they were in the text
/// they came from.
fn forget_positions(item: &mut Item) {
    match item {
        Item::Table(t) => {
            t.set_position(None);
            for (_, child) in t.iter_mut() {
                forget_positions(child);
            }
        }
        Item::ArrayOfTables(a) => {
            for t in a.iter_mut() {
                t.set_position(None);
                for (_, child) in t.iter_mut() {
                    forget_positions(child);
                }
            }
        }
        _ => {}
    }
}

/// The `[model]` of the configuration `text` alone, as TOML, comments kept: how the page
/// shows a true model.
pub fn model_of(text: &str) -> Result<String> {
    let doc: DocumentMut = text.parse().context("the configuration isn't valid TOML")?;
    let mut model = DocumentMut::new();
    model.insert(
        "model",
        doc.get("model")
            .context("the configuration has no model")?
            .clone(),
    );
    Ok(model.to_string())
}

fn round(x: f64) -> f64 {
    (x * 10.0).round() / 10.0
}

/// The two configurations of a synthetic inversion, as TOML: the start (`text`, with an
/// `[inversion]` that reads its data from `observed`) and the truth (`text` with its model
/// replaced by the `[model]` in `truth`, writing its traces to `observed`, without snapshots).
/// Both get the same time step, absorbing layers and optimised coefficients, set for both
/// models (`shared_operator`), so the observed and the synthetic traces share their sampling.
pub fn synthetic(text: &str, truth: &str, observed: &str) -> Result<(String, String)> {
    let mut start: DocumentMut = text.parse().context("the configuration isn't valid TOML")?;
    let fragment: DocumentMut = truth.parse().context("the true model isn't valid TOML")?;
    let model = fragment
        .get("model")
        .context("the true model needs [[model.layer]]s (and [[model.sphere]]s)")?;
    let extra: Vec<&str> = fragment
        .iter()
        .map(|(k, _)| k)
        .filter(|k| *k != "model")
        .collect();
    ensure!(
        extra.is_empty(),
        "the true model is only [model]; it also has {}",
        extra.join(", ")
    );
    set(start.as_table_mut(), "inversion.observed", &json!(observed))?;
    let mut truth_doc = start.clone();
    let mut model = model.clone();
    // Where [model] is in the start, not where it was in the fragment.
    forget_positions(&mut model);
    truth_doc.insert("model", model);
    truth_doc.remove("inversion");
    truth_doc.remove("snapshot");
    set(truth_doc.as_table_mut(), "output.dir", &json!(observed))?;

    let start_cfg = Config::from_toml(&start.to_string())?;
    let truth_cfg = Config::from_toml(&truth_doc.to_string()).context("in the true model")?;
    let op = shared_operator(&start_cfg, &truth_cfg, Path::new("."))?;
    for doc in [&mut start, &mut truth_doc] {
        let t = doc.as_table_mut();
        set(t, "time.dt", &json!(op.dt))?;
        set(t, "time.courant", &Json::Null)?;
        set(t, "boundary.pml_vp_max", &json!(op.pml_vp_max))?;
        if let Some(p) = op.optimised_for {
            set(t, "scheme.optimised_for", &json!(p))?;
        }
    }
    Ok((start.to_string(), truth_doc.to_string()))
}

#[cfg(test)]
mod tests {
    use super::*;

    const CONFIG: &str = r#"# A test.

[grid]
points = [36, 32, 26]  # nodes
spacing = 10.0

[time]
duration = 0.3

[boundary]
pml_width = 8

[[model.layer]]  # the only layer
top = 0.0
vp = 2400.0
vs = 1300.0
rho = 2100.0

[[source]]
position = [110.0, 155.0, 30.0]
force = [0.0, 0.0, 1e9]
wavelet = { kind = "ricker", frequency = 12.0 }

[receivers]
components = ["vx", "vz"]

[[receivers.line]]
from = [90.0, 120.0, 0.0]
to = [260.0, 120.0, 0.0]
count = 18

[output]
dir = "out/test"
"#;

    fn edited(edits: Json) -> String {
        edit(CONFIG, &edits).unwrap()
    }

    #[test]
    fn reads_settings_with_the_defaults() {
        let s = settings(CONFIG).unwrap();
        assert_eq!(s["grid"]["points"], json!([36, 32, 26]));
        assert_eq!(s["grid"]["spacing"], json!([10.0, 10.0, 10.0]));
        assert_eq!(s["scheme"]["order"], 8);
        assert_eq!(s["scheme"]["coefficients"], "optimised");
        assert_eq!(s["scheme"]["remove_time_dispersion"], true);
        assert_eq!(s["time"]["dt"], Json::Null);
        assert_eq!(s["time"]["default_courant"], 0.85);
        assert_eq!(s["boundary"]["free_surface"], true);
        assert_eq!(s["output"]["record_every"], 1);
        assert_eq!(s["sources"][0]["wavelet"], "ricker");
        assert_eq!(s["sources"][0]["frequency"], 12.0);
        assert_eq!(s["sources"][0]["mechanism"], "force");
        assert_eq!(s["receivers"]["components"], json!(["vx", "vz"]));
        assert_eq!(s["receivers"]["count"], 18);
        assert_eq!(s["model"]["layers"][0]["vs"], 1300.0);
        assert_eq!(s["inversion"], Json::Null);
        assert_eq!(s["inversion_defaults"]["iterations"], 10);
        assert_eq!(s["inversion_defaults"]["parameters"], json!(["vp", "vs"]));
        assert!(settings("[grid]\npoints = 3").is_err());
    }

    #[test]
    fn edits_keep_comments_and_layout() {
        let text = edited(json!([
            { "path": "grid.points", "value": [40, 32, 26] },
            { "path": "time.duration", "value": 0.5 },
        ]));
        assert!(text.contains("points = [40, 32, 26]  # nodes\n"), "{text}");
        assert!(text.contains("duration = 0.5\n"), "{text}");
        assert!(text.starts_with("# A test.\n"));
        assert!(text.contains("[[model.layer]]  # the only layer"));
        // Nothing else moved.
        assert_eq!(text.lines().count(), CONFIG.lines().count());
        let c = Config::from_toml(&text).unwrap();
        assert_eq!(c.grid.points, [40, 32, 26]);
    }

    #[test]
    fn adds_missing_tables_and_keys_and_removes_them() {
        let text = edited(json!([
            { "path": "scheme.order", "value": 4 },
            { "path": "scheme.coefficients", "value": "taylor" },
            { "path": "time.courant", "value": 0.5 },
            { "path": "boundary.free_surface", "value": false },
            { "path": "output.record_every", "value": 2 },
        ]));
        let c = Config::from_toml(&text).unwrap();
        assert_eq!(c.scheme.order, 4);
        assert_eq!(c.scheme.coefficients, Coefficients::Taylor);
        assert_eq!(c.time.courant, Some(0.5));
        assert!(!c.boundary.free_surface);
        assert_eq!(c.output.record_every, 2);
        assert!(text.contains("\n[scheme]\norder = 4\n"), "{text}");
        // A Courant fraction, then a time step instead.
        let text = edit(
            &text,
            &json!([{ "path": "time.dt", "value": 1e-3 }, { "path": "time.courant", "value": null }]),
        )
        .unwrap();
        let c = Config::from_toml(&text).unwrap();
        assert_eq!((c.time.dt, c.time.courant), (Some(1e-3), None));
        // Removing what isn't there is nothing.
        assert_eq!(
            edit(CONFIG, &json!([{ "path": "scheme.order", "value": null }])).unwrap(),
            CONFIG
        );
    }

    #[test]
    fn numbers_are_typed_by_their_key() {
        let text = edited(json!([
            { "path": "source.0.wavelet.frequency", "value": 15 },
            { "path": "grid.spacing", "value": 12 },
        ]));
        assert!(text.contains("frequency = 15.0 }"), "{text}");
        assert!(text.contains("spacing = 12.0\n"), "{text}");
        let e = edit(CONFIG, &json!([{ "path": "scheme.order", "value": 4.5 }]));
        assert!(format!("{:#}", e.unwrap_err()).contains("whole number"));
    }

    #[test]
    fn edits_inline_tables_and_arrays_of_tables() {
        let text = edited(json!([
            { "path": "source.0.wavelet.kind", "value": "gaussian" },
            { "path": "source.0.position", "value": [100.0, 150.0, 40.0] },
            { "path": "model.layer.0.vs", "value": 1200.0 },
            { "path": "model.layer.1", "value": { "top": 100.0, "vp": 3000.0, "vs": 1700.0, "rho": 2300.0 } },
            { "path": "model.sphere.0", "value": { "center": [180.0, 160.0, 130.0], "radius": 30.0, "vp": 2600.0, "vs": 1400.0, "rho": 2100.0 } },
            { "path": "source", "copy": 0 },
            { "path": "source.1.wavelet.frequency", "value": 10.0 },
        ]));
        let c = Config::from_toml(&text).unwrap();
        assert!(
            matches!(c.sources[0].wavelet, WaveletConfig::Gaussian { frequency, .. } if frequency == 12.0)
        );
        assert_eq!(c.sources[0].position, [100.0, 150.0, 40.0]);
        assert_eq!(c.model.layers.len(), 2);
        assert_eq!(c.model.layers[0].vs, 1200.0);
        assert_eq!(c.model.layers[1].top, 100.0);
        assert_eq!(c.model.spheres[0].radius, 30.0);
        assert_eq!(c.sources.len(), 2);
        assert!(
            text.contains("wavelet = { kind = \"gaussian\", frequency = 10.0 }"),
            "{text}"
        );
        assert!(
            text.contains("[[model.layer]]\ntop = 100.0\nvp = 3000.0\nvs = 1700.0\nrho = 2300.0\n"),
            "{text}"
        );
        assert_eq!(c.sources[1].force, Some([0.0, 0.0, 1e9]));
        // New layers follow the old ones, not the end of the file.
        let layers: Vec<usize> = text.match_indices("[[model.").map(|(i, _)| i).collect();
        assert!(
            layers.iter().all(|&i| i < text.find("[[source]]").unwrap()),
            "{text}"
        );
        // Removing the last of an array removes the array.
        let text = edit(&text, &json!([{ "path": "model.sphere.0", "value": null }])).unwrap();
        assert!(!text.contains("sphere"), "{text}");
        let e = edit(
            &text,
            &json!([{ "path": "model.layer.5", "value": { "top": 1.0 } }]),
        );
        assert!(e.is_err());
    }

    #[test]
    fn adds_an_inversion_section() {
        let text = edited(json!([
            { "path": "inversion.observed", "value": "observed" },
            { "path": "inversion.parameters", "value": ["vs"] },
            { "path": "inversion.iterations", "value": 3 },
            { "path": "inversion.vs_range", "value": [600.0, 2000.0] },
        ]));
        let c = Config::from_toml(&text).unwrap();
        let i = c.inversion.unwrap();
        assert_eq!(i.iterations, 3);
        assert_eq!(i.vs_range, Some([600.0, 2000.0]));
        assert!(
            text.contains("\n[inversion]\nobserved = \"observed\"\n"),
            "{text}"
        );
        let s = settings(&text).unwrap();
        assert_eq!(s["inversion"]["parameters"], json!(["vs"]));
    }

    #[test]
    fn writes_the_commands() {
        let c = commands(&json!({ "file": "explosion.toml", "target": "me@box", "sources": 1 }))
            .unwrap();
        let all: Vec<&str> = c
            .as_array()
            .unwrap()
            .iter()
            .map(|c| c["command"].as_str().unwrap())
            .collect();
        assert!(all.contains(&"vibe-rations run explosion.toml"));
        assert!(all.contains(&"vibe-rations run --remote me@box --detach explosion.toml"));
        assert!(all.contains(&"just run explosion.toml"));
        assert!(
            !all.iter()
                .any(|c| c.contains("--shots") || c.contains("attach"))
        );
        let c = commands(&json!({
            "file": "my run.toml", "target": "http://localhost:8089", "run": "explosion-1-ab",
            "sources": 2, "start": "my-start.toml", "truth": "my-truth.toml",
        }))
        .unwrap();
        let all: Vec<&str> = c
            .as_array()
            .unwrap()
            .iter()
            .map(|c| c["command"].as_str().unwrap())
            .collect();
        assert!(all.contains(&"vibe-rations run 'my run.toml'"), "{all:?}");
        assert!(all.contains(&"vibe-rations remote http://localhost:8089 attach explosion-1-ab"));
        assert!(all.contains(&"vibe-rations run --shots 'my run.toml'"));
        assert!(all.contains(&"vibe-rations invert my-start.toml"));
        assert!(all.contains(&"vibe-rations run --shots my-truth.toml"));
        // The gradient command takes one source's traces.
        assert!(!all.iter().any(|c| c.contains("gradient")));
        assert_eq!(quote("it's"), r"'it'\''s'");
    }

    #[test]
    fn makes_a_true_model_and_the_synthetic_pair() {
        let truth = default_truth(CONFIG).unwrap();
        assert!(
            truth.contains("[[model.layer]]  # the only layer"),
            "{truth}"
        );
        assert!(truth.contains("[[model.sphere]]"), "{truth}");
        let (start, true_text) = synthetic(CONFIG, &truth, "test-observed").unwrap();
        let s = Config::from_toml(&start).unwrap();
        let t = Config::from_toml(&true_text).unwrap();
        assert_eq!(s.model.spheres.len(), 0);
        assert_eq!(t.model.spheres.len(), 1);
        assert_eq!(t.model.spheres[0].vs, 1430.0);
        // Below the source (at 30 m) and clear of it, above the bottom layer (from 170 m).
        let sphere = &t.model.spheres[0];
        assert!(sphere.center[2] - sphere.radius >= 50.0, "{sphere:?}");
        assert!(sphere.center[2] + sphere.radius <= 170.0, "{sphere:?}");
        assert!(t.inversion.is_none());
        assert_eq!(t.output.dir, Path::new("test-observed"));
        // The true model where the start has its own.
        let at = |k: &str| true_text.find(k).unwrap();
        assert!(at("[boundary]") < at("[[model.layer]]"), "{true_text}");
        assert!(at("[[model.sphere]]") < at("[[source]]"), "{true_text}");
        let i = s.inversion.unwrap();
        assert_eq!(i.observed, Path::new("test-observed"));
        // The same operator for both, set for the faster model.
        assert!(s.time.dt.is_some() && s.time.dt == t.time.dt && s.time.courant.is_none());
        assert_eq!(s.boundary.pml_vp_max, Some(2640.0));
        assert_eq!(t.boundary.pml_vp_max, s.boundary.pml_vp_max);
        assert!(s.scheme.optimised_for.is_some());
        assert_eq!(s.scheme.optimised_for, t.scheme.optimised_for);
        // The start's own [inversion] stays.
        let with = edit(
            CONFIG,
            &json!([
                { "path": "inversion.observed", "value": "x" },
                { "path": "inversion.iterations", "value": 2 },
            ]),
        )
        .unwrap();
        let (start, _) = synthetic(&with, &truth, "obs").unwrap();
        let i = Config::from_toml(&start).unwrap().inversion.unwrap();
        assert_eq!((i.iterations, i.observed.to_str()), (2, Some("obs")));
        assert!(synthetic(CONFIG, "[grid]\npoints = [1, 2, 3]", "obs").is_err());
        // The truth's model alone, as the page shows it, makes the same pair again.
        let shown = model_of(&true_text).unwrap();
        assert!(
            shown.contains("[[model.layer]]  # the only layer"),
            "{shown}"
        );
        assert!(!shown.contains("[grid]") && shown.contains("[[model.sphere]]"));
        assert_eq!(
            synthetic(CONFIG, &shown, "test-observed").unwrap().1,
            true_text
        );
    }
}
