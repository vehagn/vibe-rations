// The Command tab: the command lines that run the configuration being edited without the page
// (here, in the repository, on another machine, or the run on the page's server), and the files
// they refer to.

import { commands, synthetic } from "./config.js";
import { $, el, copy, download } from "./util.js";

const TARGET_KEY = "vibe-rations.command-target";
let hooks = { text: () => "", run: () => null, inversionRun: () => null, server: () => false, truth: () => "", sources: () => 1, balance: () => true };

export function initCommands(h) {
  hooks = h;
  let saved = null;
  try {
    saved = localStorage.getItem(TARGET_KEY);
  } catch {}
  $("#command-target").value = saved ?? "";
  $("#command-target").addEventListener("input", (e) => {
    try {
      localStorage.setItem(TARGET_KEY, e.target.value);
    } catch {}
    refreshCommands();
  });
  $("#command-file").addEventListener("input", refreshCommands);
  $("#download-toml").addEventListener("click", () =>
    download(new Blob([hooks.text()], { type: "text/plain" }), file()),
  );
  $("#download-fwi").addEventListener("click", () => {
    const pair = fwiPair();
    if (!pair) return;
    download(new Blob([pair.start], { type: "text/plain" }), pair.files.start);
    // Two downloads at once are sometimes refused as one too many; space them.
    setTimeout(() => download(new Blob([pair.truth], { type: "text/plain" }), pair.files.truth), 400);
  });
}

/// A new configuration's file name: the example's, until the visitor names it.
export function nameCommands(name) {
  $("#command-file").value = name.endsWith(".toml") ? name : `${name}.toml`;
}

const file = () => $("#command-file").value.trim() || "run.toml";
const stem = () => file().replace(/\.toml$/, "");

/// The machine for remote commands: what the visitor typed, or this server when it served the
/// page (its URL works for the command line on this machine), or a placeholder.
function target() {
  const typed = $("#command-target").value.trim();
  if (typed) return typed;
  return hooks.server() ? location.origin : "user@host";
}

/// The inversion's two files and their contents, or null (and why, shown) if they can't be made.
function fwiPair() {
  const files = { start: `${stem()}-start.toml`, truth: `${stem()}-truth.toml`, observed: `${stem()}-observed` };
  try {
    const pair = synthetic(hooks.text(), hooks.truth(), files.observed);
    $("#fwi-commands-error").hidden = true;
    return { ...pair, files };
  } catch (e) {
    $("#fwi-commands-error").textContent = `No inversion files yet: ${e.message ?? e}`;
    $("#fwi-commands-error").hidden = false;
    return null;
  }
}

export function refreshCommands() {
  if ($('[data-panel="command"]').hidden) return;
  const input = { file: file(), target: target(), run: hooks.run(), sources: hooks.sources() };
  list("#commands", commands(input));
  const pair = fwiPair();
  $("#download-fwi").disabled = !pair;
  const items = pair
    ? commands({ ...input, run: null, start: pair.files.start, truth: pair.files.truth, observed: pair.files.observed, balance: hooks.balance() }).filter((c) =>
        c.what.startsWith("Inversion"),
      )
    : [];
  // The inversion on the page's server, if there is one: follow it or fetch its model.
  const id = hooks.inversionRun?.();
  if (id) {
    for (const c of commands({ ...input, run: id })) {
      if (/ (attach|fetch) /.test(c.command)) items.push({ ...c, what: c.what.replace("this run", "this inversion") });
    }
  }
  list("#fwi-commands", items);
}

function list(selector, items) {
  $(selector).replaceChildren(
    ...items.map(({ what, command }) => {
      const button = el("button", { className: "copy", textContent: "Copy", title: "Copy the command" });
      button.addEventListener("click", () => copy(command, button));
      return el("li", {}, [el("span", { className: "what", textContent: capital(what.replace(/^Inversion: /, "")) }), el("div", { className: "line" }, [el("code", { textContent: command }), button])]);
    }),
  );
}

const capital = (s) => s.charAt(0).toUpperCase() + s.slice(1);
