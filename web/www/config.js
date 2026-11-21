// The configuration's text, read and edited on the page's own thread by the single-threaded
// WebAssembly build (web/src/settings.rs), whichever solver runs it: the form, the commands
// and the inversion's pair of configurations all go through here.

import init, * as wasm from "./pkg/vibe_rations_web.js";

export const ready = init();

/// What a configuration sets, defaults filled in; throws if it doesn't parse.
export const readSettings = (text) => JSON.parse(wasm.readSettings(text));

/// The configuration with `edits` ([{path, value}] or [{path, copy}]) applied, comments kept.
export const editConfig = (text, edits) => wasm.editConfig(text, JSON.stringify(edits));

/// [{what, command}] that run the configuration without the page.
export const commands = (input) => JSON.parse(wasm.commands(JSON.stringify(input)));

/// The configuration's model with a faster sphere added, as TOML: a true model to start from.
export const defaultTruth = (text) => wasm.defaultTruth(text);

/// The [model] of a configuration alone, as TOML.
export const modelOf = (text) => wasm.modelOf(text);

/// {start, truth}: the two configurations of an inversion of synthetic data.
export const synthetic = (text, truth, observed) => JSON.parse(wasm.synthetic(text, truth, observed));
