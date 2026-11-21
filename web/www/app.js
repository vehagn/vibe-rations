// The page: picks a configuration, hands it to a solver, and draws what comes back. The solver
// is WebAssembly in worker.js, on this machine, or, when `vibe-rations serve` served the page,
// the server (remote.js); both answer the same messages. The configuration's text is the source
// of truth: the settings form (form.js) and the inversion (fwi.js) edit it through config.js,
// and the Command tab (commands.js) shows how to run it without the page.

import { RemoteSolver, findServer, listRuns } from "./remote.js";
import { ready as configReady, defaultTruth, modelOf } from "./config.js";
import { initForm, refreshForm, settings } from "./form.js";
import { initCommands, nameCommands, refreshCommands } from "./commands.js";
import * as fwi from "./fwi.js";
import { plot, fitHeight, clearCanvas } from "./plot.js";
import {
  $,
  el,
  download,
  formatBytes,
  formatDuration,
  formatGain,
  formatLength,
  formatValue,
  showMemory,
} from "./util.js";

const UNITS = { vx: "m/s", vy: "m/s", vz: "m/s", p: "Pa" };
const NORMAL_AXIS = { xy: 2, xz: 1, yz: 0 };

const state = {
  info: null,
  running: false,
  view: { field: "vz", plane: "xz", index: 0, gain: 1, gather: null, gatherRows: 600 },
  section: null,
  gather: null,
  clip: 0,
  step: 0,
  totalSteps: 0,
  // The configuration the solver has set up, to tell edits from it.
  loaded: null,
  // The run on the server, and the inversion there, for the commands.
  runId: null,
  inversionId: null,
};

// ?threads=N picks the solver's thread count, e.g. 1 to compare.
const workerUrl = new URL("./worker.js", import.meta.url);
const askedThreads = new URLSearchParams(location.search).get("threads");
if (askedThreads) workerUrl.searchParams.set("threads", askedThreads);
let worker;

/// Starts the solver: "browser" (WebAssembly here) or "server" (the one that served the page).
function useSolver(kind) {
  worker?.terminate();
  worker = kind === "server" ? new RemoteSolver() : new Worker(workerUrl, { type: "module" });
  worker.onmessage = ({ data }) => on[data.type]?.(data);
  worker.onerror = (e) => showError(`The solver failed to start: ${e.message}`);
  state.compute = kind;
  state.runId = state.inversionId = null;
  $("#server-runs-field").hidden = kind !== "server";
}

// The example list; set at the end of this module, before any message can arrive.
let examplesLoaded;

const on = {
  async ready({ threads, isolated, server }) {
    await examplesLoaded;
    $("#stat-threads").textContent =
      (threads > 1 ? `${threads}` : isolated ? "1" : "1 (page not isolated)") + (server ? " on the server" : "");
    setStatus("ready", "Solver ready");
    worker.postMessage({ type: "live", on: $("#live").checked });
    if (server) refreshRuns();
    // Switching solvers keeps the configuration being edited.
    if (state.switched) load(editor().value);
    else loadExample($("#example").value);
    state.switched = true;
  },
  loaded({ info, toml, name, id }) {
    if (toml !== undefined) {
      // An attached run: show its configuration, and its example if it is one.
      if ([...$("#example").options].some((o) => o.value === name)) $("#example").value = name;
      setText(toml);
      state.loaded = toml;
      $("#example-description").textContent = describe(toml);
    }
    state.runId = id ?? null;
    if (state.compute === "server") refreshRuns();
    state.info = info;
    clearError();
    setupView(info);
    showInfo(info);
    state.step = 0;
    state.totalSteps = info.steps;
    updateProgress({ step: 0, totalSteps: info.steps, time: 0, speed: 0 });
    setStatus("ready", "Ready to run");
    markEdited();
    refreshCommands();
    sendView();
  },
  state({ running, done }) {
    if (state.compute === "server" && state.running && !running) refreshRuns();
    state.running = running;
    $("#run").textContent = running ? "Pause" : state.step > 0 && !done ? "Resume" : "Run";
    $("#run").disabled = !state.info || (done ?? false) || state.step >= state.totalSteps;
    if (running) setStatus("running", "Running");
    else if (done || (state.totalSteps && state.step >= state.totalSteps)) setStatus("done", "Done");
    else if (state.info) setStatus("ready", state.step > 0 ? "Paused" : "Ready to run");
  },
  frame(f) {
    // Without live pictures, progress only: the last pictures stay.
    if (f.section) {
      state.section?.bitmap.close();
      state.section = f.section;
      state.clip = f.clip;
    }
    if (f.gather) {
      state.gather?.bitmap.close();
      state.gather = f.gather;
    }
    state.step = f.step;
    updateProgress(f);
    draw();
  },
  estimate(message) {
    if (String(message.id).startsWith("fwi-")) return fwi.showEstimate(message);
    if (message.id === state.estimate) showMemory($("#memory"), message.memory, message.error);
  },
  inversion(message) {
    if (message.id) state.inversionId = message.id;
    fwi.showStatus(message);
    if (message.running) setStatus("running", "Inverting");
    else if (message.status?.stage === "done") setStatus("done", "Inversion done");
    else if (message.status) setStatus("ready", "Inversion stopped");
    refreshCommands();
  },
  inversionFrame(f) {
    fwi.showFrame(f);
  },
  inversionAttached(message) {
    state.inversionId = message.id;
    setText(message.start);
    let truth = message.truth;
    try {
      truth = modelOf(truth);
    } catch {}
    fwi.attached({ ...message, truth });
    showTab("fwi");
  },
  modelNpy(message) {
    fwi.savedModel(message);
  },
  npy({ component, bytes }) {
    download(new Blob([bytes], { type: "application/octet-stream" }), `${component}.npy`);
  },
  error({ during, message }) {
    const what = {
      load: "The configuration has a problem",
      run: "The run stopped",
      invert: "The inversion stopped",
    }[during];
    showError(what ? `${what}: ${message}` : message);
    if (during === "load") setStatus("error", "Configuration error");
    if (during === "run" || during === "invert") setStatus("error", "Stopped");
  },
};

const editor = () => $("#editor");

// ---------------------------------------------------------------------------- scenarios

async function loadExamples() {
  const select = $("#example");
  const files = await fetch("examples/index.json").then((r) => r.json());
  for (const file of files) {
    select.append(el("option", { value: file, textContent: file.replace(/\.toml$/, "").replace(/[-_]/g, " ") }));
  }
  const wanted = new URLSearchParams(location.search).get("example");
  select.value = files.includes(`${wanted}.toml`) ? `${wanted}.toml` : files.includes("explosion.toml") ? "explosion.toml" : files[0];
}

async function loadExample(file) {
  const toml = await fetch(`examples/${file}`).then((r) => r.text());
  $("#example-description").textContent = describe(toml);
  nameCommands(file);
  editor().value = toml;
  fwi.resetTruth();
  setText(toml);
  load(toml);
}

/// The first paragraph of the comment block at the top of a configuration.
function describe(toml) {
  const lines = [];
  for (const line of toml.split("\n")) {
    if (!line.startsWith("#")) break;
    const text = line.replace(/^#\s?/, "").trim();
    if (!text && lines.length) break;
    if (text) lines.push(text);
  }
  return lines.join(" ");
}

/// Puts a configuration in the editor and has everything that follows it catch up.
function setText(text) {
  editor().value = text;
  edited();
}

/// The configuration has changed: the form, the memory, the commands and the inversion's memory
/// follow it, at once for small changes and after a pause in typing.
let following;
function edited({ wait = 0 } = {}) {
  clearTimeout(following);
  const run = () => {
    const text = editor().value;
    estimate(text);
    refreshForm(text);
    markEdited();
    refreshCommands();
    fwi.estimateFwi();
  };
  if (wait) following = setTimeout(run, wait);
  else run();
}

/// Shows Apply and Revert while the configuration differs from the one set up.
function markEdited() {
  $("#changes").hidden = state.loaded === null || editor().value === state.loaded;
}

/// Asks the solver what a configuration would need in memory; the latest answer is shown.
function estimate(toml) {
  state.estimate = (state.estimate ?? 0) + 1;
  worker.postMessage({ type: "estimate", toml, id: state.estimate });
}

function load(toml) {
  estimate(toml);
  state.info = null;
  state.loaded = toml;
  markEdited();
  state.section?.bitmap.close();
  state.gather?.bitmap.close();
  state.section = state.gather = null;
  $("#run").disabled = true;
  setStatus("busy", "Setting up");
  worker.postMessage({ type: "load", toml, name: $("#example").value || "page.toml" });
}

/// Fills the list of runs on the server, to attach to.
async function refreshRuns() {
  const select = $("#server-runs");
  try {
    const runs = await listRuns();
    const progress = (r) => (r.total_steps ? ` ${Math.round((100 * r.step) / r.total_steps)}%` : "");
    select.replaceChildren(
      el("option", { value: "", textContent: `${runs.length} run(s) on the server…` }),
      ...runs
        .filter((r) => r.info)
        .map((r) =>
          el("option", {
            value: r.id,
            textContent:
              `${r.name.replace(/\.toml$/, "")}${r.kind === "inversion" ? " · inversion" : ""} · ${r.phase}` +
              `${r.phase === "done" ? "" : progress(r)}${r.live ? "" : " (results only)"}`,
            disabled: !r.live,
          }),
        ),
    );
  } catch (e) {
    showError(`Listing the server's runs: ${e.message}`);
  }
}

// ---------------------------------------------------------------------------- tabs

function showTab(name) {
  for (const tab of document.querySelectorAll('[role="tab"]')) tab.setAttribute("aria-selected", tab.dataset.tab === name);
  for (const panel of document.querySelectorAll(".tab-panel")) panel.hidden = panel.dataset.panel !== name;
  if (name === "command") refreshCommands();
  if (name === "fwi") fwi.estimateFwi();
  try {
    localStorage.setItem("vibe-rations.tab", name);
  } catch {}
}

for (const tab of document.querySelectorAll('[role="tab"]')) {
  tab.addEventListener("click", () => showTab(tab.dataset.tab));
}

// ---------------------------------------------------------------------------- view

function setupView(info) {
  const source = info.sources[0]?.position ?? [0, 0, 0];
  // Through the source along the plane's normal.
  state.view.index = Math.round(source[NORMAL_AXIS[state.view.plane]] / info.spacing[NORMAL_AXIS[state.view.plane]]);
  const gather = $("#gather-component");
  gather.replaceChildren(...info.components.map((c) => el("option", { value: c, textContent: c })));
  state.view.gather = info.components.includes(state.view.field) ? state.view.field : info.components[0] ?? null;
  gather.value = state.view.gather ?? "";
  gather.disabled = !state.view.gather;
  $("#gather-panel").hidden = info.receivers < 2;
  updatePositionSlider();
}

function updatePositionSlider() {
  const info = state.info;
  if (!info) return;
  const axis = NORMAL_AXIS[state.view.plane];
  const slider = $("#position");
  slider.max = info.points[axis] - 1;
  state.view.index = Math.min(state.view.index, info.points[axis] - 1);
  slider.value = state.view.index;
  $("#position-label").textContent = `${["x", "y", "z"][axis]} = ${formatLength(state.view.index * info.spacing[axis])}`;
}

function sendView() {
  worker.postMessage({ type: "view", view: state.view });
}

for (const input of document.querySelectorAll('input[name="field"]')) {
  input.addEventListener("change", () => {
    state.view.field = input.value;
    sendView();
  });
}

for (const input of document.querySelectorAll('input[name="plane"]')) {
  input.addEventListener("change", () => {
    state.view.plane = input.value;
    if (state.info) {
      const axis = NORMAL_AXIS[input.value];
      const source = state.info.sources[0]?.position ?? [0, 0, 0];
      state.view.index = Math.round(source[axis] / state.info.spacing[axis]);
    }
    updatePositionSlider();
    sendView();
  });
}

$("#position").addEventListener("input", (e) => {
  state.view.index = Number(e.target.value);
  updatePositionSlider();
  sendView();
});

$("#gain").addEventListener("input", (e) => {
  state.view.gain = 2 ** Number(e.target.value);
  $("#gain-label").textContent = `×${formatGain(state.view.gain)}`;
  sendView();
  fwi.sendView();
});

$("#gather-component").addEventListener("change", (e) => {
  state.view.gather = e.target.value;
  sendView();
});

// ---------------------------------------------------------------------------- controls

$("#example").addEventListener("change", (e) => loadExample(e.target.value));
$("#apply").addEventListener("click", () => {
  $("#example-description").textContent = describe(editor().value);
  load(editor().value);
});
$("#revert").addEventListener("click", () => setText(state.loaded));
$("#reset").addEventListener("click", () => load(editor().value));
for (const input of document.querySelectorAll('input[name="compute"]')) {
  input.addEventListener("change", () => {
    if (!input.checked) return;
    state.running = false;
    setStatus("busy", "Switching the solver");
    useSolver(input.value);
  });
}
$("#server-runs").addEventListener("focus", refreshRuns);
$("#live").addEventListener("change", (e) => worker.postMessage({ type: "live", on: e.target.checked }));
$("#server-runs").addEventListener("change", (e) => {
  if (!e.target.value) return;
  state.section?.bitmap.close();
  state.gather?.bitmap.close();
  state.section = state.gather = null;
  setStatus("busy", "Attaching");
  worker.postMessage({ type: "attach", id: e.target.value });
  e.target.value = "";
});
// The memory estimate, the form and the commands follow the configuration as it is typed.
editor().addEventListener("input", () => edited({ wait: 300 }));
$("#run").addEventListener("click", () => worker.postMessage({ type: state.running ? "pause" : "run" }));
document.addEventListener("keydown", (e) => {
  if (e.key === " " && !["TEXTAREA", "INPUT", "SELECT", "BUTTON"].includes(document.activeElement?.tagName)) {
    e.preventDefault();
    if (!$("#run").disabled) $("#run").click();
  }
});
$("#download-section").addEventListener("click", () => {
  const { field, plane, index } = state.view;
  $("#section").toBlob((blob) => download(blob, `${field}_${plane}${index}_${String(state.step).padStart(6, "0")}.png`));
});
$("#download-traces").addEventListener("click", () => {
  if (state.view.gather) worker.postMessage({ type: "npy", component: state.view.gather });
});
new ResizeObserver(() => draw()).observe($("#views"));

// ---------------------------------------------------------------------------- drawing

function draw() {
  drawSection();
  drawGather();
}

function drawSection() {
  const { info, section, view } = state;
  const canvas = $("#section");
  if (!info || !section) return clearCanvas(canvas);
  const [nx, ny, nz] = info.points;
  const [dx, dy, dz] = info.spacing;
  const extents = { x: (nx - 1) * dx, y: (ny - 1) * dy, z: (nz - 1) * dz };
  const [h, v] = { xz: ["x", "z"], xy: ["x", "y"], yz: ["y", "z"] }[view.plane];
  const axis = NORMAL_AXIS[view.plane];
  const near = (p) => Math.abs(p[axis] / info.spacing[axis] - view.index) <= 1;
  const coord = { x: 0, y: 1, z: 2 };
  const markers = [
    ...info.receiver_positions.filter(near).map((p) => ({ x: p[coord[h]], y: p[coord[v]], kind: "receiver" })),
    ...info.sources.map((s) => s.position).filter(near).map((p) => ({ x: p[coord[h]], y: p[coord[v]], kind: "source" })),
  ];
  fitHeight(canvas, extents[h] / extents[v]);
  plot(canvas, section.bitmap, {
    x: [0, extents[h], h === "z" ? "depth" : h],
    y: [0, extents[v], v === "z" ? "depth" : v],
    aspect: extents[h] / extents[v],
    smooth: true,
    markers,
    format: formatLength,
  });
  const unit = UNITS[view.field];
  const clip = formatValue(state.clip);
  $("#legend-min").textContent = state.clip > 0 ? `−${clip} ${unit}` : `0 ${unit}`;
  $("#legend-max").textContent = state.clip > 0 ? `+${clip} ${unit}` : `0 ${unit}`;
  $("#section-title").textContent = `${view.field} on the ${view.plane} plane at ${["x", "y", "z"][axis]} = ${formatLength(view.index * info.spacing[axis])}`;
}

function drawGather() {
  const { info, gather } = state;
  const canvas = $("#gather");
  if (!info || !gather || info.receivers < 2) return clearCanvas(canvas);
  plot(canvas, gather.bitmap, {
    x: [1, info.receivers, "receiver", (v) => String(Math.round(v))],
    y: [0, info.duration, "time", (v) => `${+v.toFixed(3)} s`],
    smooth: false,
    markers: [],
    format: formatLength,
  });
}

// ---------------------------------------------------------------------------- panels

function showInfo(info) {
  const [nx, ny, nz] = info.points;
  const rows = [
    ["Grid", `${nx} × ${ny} × ${nz} nodes`],
    ["Spacing", info.spacing.map((h) => `${h}`).join(" × ") + " m"],
    ["Extent", info.extent.map(formatLength).join(" × ")],
    ["vp", `${info.vp_range.map(Math.round).join("–")} m/s`],
    ["vs", `${info.vs_range.map(Math.round).join("–")} m/s`],
    ["Time step", `${(info.dt * 1e3).toFixed(3)} ms`],
    ["Steps", `${info.steps} to ${info.duration.toFixed(3)} s`],
    [
      "Accuracy",
      `${info.points_per_wavelength.toFixed(1)} points per wavelength, order ${info.order}` +
        (info.band ? ", optimised coefficients" : "") +
        (info.remove_time_dispersion ? ", time dispersion removed" : ""),
    ],
    ["Boundary", `${info.pml_width}-node layers${info.free_surface ? ", free surface" : ""}`],
    ["Receivers", `${info.receivers} × ${info.samples} samples`],
    ["Memory", formatBytes(info.memory_bytes)],
  ];
  $("#info").replaceChildren(...rows.flatMap(([k, v]) => [el("dt", { textContent: k }), el("dd", { textContent: v })]));
  $("#warnings").replaceChildren(...info.warnings.map((w) => el("li", { textContent: w })));
}

function updateProgress({ step, totalSteps, time, speed }) {
  $("#progress").value = totalSteps ? step / totalSteps : 0;
  $("#stat-time").textContent = `${time.toFixed(3)} / ${(state.info?.duration ?? 0).toFixed(3)} s`;
  $("#stat-step").textContent = `${step} / ${totalSteps}`;
  $("#stat-speed").textContent = speed > 0 ? `${speed.toFixed(0)} Mcells/s` : "–";
  const cells = state.info ? state.info.points.reduce((a, b) => a * b, 1) : 0;
  const left = speed > 0 ? ((totalSteps - step) * cells) / (speed * 1e6) : NaN;
  $("#stat-left").textContent = Number.isFinite(left) && step < totalSteps ? formatDuration(left) : "–";
}

function setStatus(kind, text) {
  const status = $("#status");
  status.dataset.kind = kind;
  status.textContent = text;
}

function showError(message) {
  const box = $("#error");
  box.textContent = message;
  box.hidden = false;
}

function clearError() {
  $("#error").hidden = true;
}

// ---------------------------------------------------------------------------- start

const stem = () => ($("#command-file").value.trim() || "run.toml").replace(/\.toml$/, "");

examplesLoaded = loadExamples();
await configReady;
initForm({ text: () => editor().value, change: setText, info: () => state.info });
initCommands({
  text: () => editor().value,
  run: () => state.runId,
  inversionRun: () => state.inversionId,
  server: () => state.compute === "server",
  truth: () => $("#truth").value,
  sources: () => settings()?.sources.length ?? 1,
  balance: () => (settings()?.inversion ?? settings()?.inversion_defaults)?.balance ?? true,
});
fwi.initFwi({
  solver: () => worker,
  text: () => editor().value,
  truthDefault: () => {
    try {
      return defaultTruth(editor().value);
    } catch (e) {
      return `# ${e.message ?? e}\n`;
    }
  },
  name: stem,
  gain: () => state.view.gain,
  markers: () => [
    ...(settings()?.sources ?? []).map((s) => ({ position: s.position, kind: "source" })),
    ...(state.info?.receiver_positions ?? []).map((p) => ({ position: p, kind: "receiver" })),
  ],
});
try {
  const tab = localStorage.getItem("vibe-rations.tab");
  if (tab && document.querySelector(`[data-panel="${tab}"]`)) showTab(tab);
} catch {}

// Served by `vibe-rations serve`: offer its solver, and use it.
const server = await findServer();
if (server) {
  $("#compute").hidden = false;
  $('input[name="compute"][value="server"]').checked = true;
  useSolver("server");
} else {
  useSolver("browser");
}
await examplesLoaded;
