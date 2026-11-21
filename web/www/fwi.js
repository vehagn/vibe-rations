// The FWI tab and the inversion panel: an inversion of synthetic data (the configuration's model
// the start, the true model's traces the data), run by the solver (worker.js here, remote.js on
// the server) a unit at a time, with its misfit per iteration and pictures of the model, its
// update, the true update and the gradient.

import { editConfig, readSettings, synthetic } from "./config.js";
import { misfitChart, plot, fitHeight, clearCanvas } from "./plot.js";
import { $, el, download, formatLength, formatSigned, showMemory } from "./util.js";

const NORMAL_AXIS = { xy: 2, xz: 1, yz: 0 };
const UNITS = { vp: "m/s", vs: "m/s", rho: "kg/m³" };
const LABEL = { vp: "vp", vs: "vs", rho: "ρ" };

let hooks = { solver: () => null, text: () => "", truthDefault: () => "", name: () => "run", gain: () => 1, markers: () => [] };
const state = {
  status: null,
  running: false,
  grid: null,
  frame: null,
  view: { parameter: "vs", show: "update", plane: "xz", index: 0 },
  estimate: 0,
};

export function initFwi(h) {
  hooks = h;
  $("#fwi-run").addEventListener("click", () => start(null));
  $("#fwi-gradient").addEventListener("click", () => start(0));
  $("#fwi-stop").addEventListener("click", () => {
    const solver = hooks.solver();
    if (state.running) solver.postMessage({ type: "stopInversion" });
    else solver.postMessage({ type: "resumeInversion" });
  });
  $("#truth-reset").addEventListener("click", () => {
    $("#truth").value = hooks.truthDefault();
    estimateFwi();
  });
  let typing;
  $("#truth").addEventListener("input", () => {
    clearTimeout(typing);
    typing = setTimeout(estimateFwi, 400);
  });
  for (const input of document.querySelectorAll('input[name="fwi-parameter"]')) {
    input.addEventListener("change", () => setView({ parameter: input.value }));
  }
  $("#fwi-show").addEventListener("change", (e) => setView({ show: e.target.value }));
  for (const input of document.querySelectorAll('input[name="fwi-plane"]')) {
    input.addEventListener("change", () => {
      state.view.plane = input.value;
      state.view.index = middle(input.value);
      updateSlider();
      sendView();
    });
  }
  $("#fwi-position").addEventListener("input", (e) => {
    state.view.index = Number(e.target.value);
    updateSlider();
    sendView();
  });
  $("#fwi-download").addEventListener("click", () =>
    hooks.solver().postMessage({ type: "modelNpy", parameter: state.view.parameter }),
  );
  new ResizeObserver(() => draw()).observe($("#inversion-panel"));
}

/// The true model for a new configuration: its own model with a sphere.
export function resetTruth() {
  $("#truth").value = hooks.truthDefault();
}

/// The inversion's two configurations, with the gradient's iterations if `iterations` is set.
function pair(iterations) {
  const observed = `${hooks.name()}-observed`;
  const { start, truth } = synthetic(hooks.text(), $("#truth").value, observed);
  return {
    start: iterations == null ? start : editConfig(start, [{ path: "inversion.iterations", value: iterations }]),
    truth,
  };
}

function start(iterations) {
  let configs;
  try {
    configs = pair(iterations);
  } catch (e) {
    return showError(e.message ?? String(e));
  }
  const s = readSettings(configs.start);
  state.grid = s.grid;
  state.view.index = Math.min(state.view.index || middle(state.view.plane), s.grid.points[NORMAL_AXIS[state.view.plane]] - 1);
  if (!state.view.index) state.view.index = middle(state.view.plane);
  state.view.parameter = s.inversion.parameters.includes(state.view.parameter) ? state.view.parameter : s.inversion.parameters[0];
  document.querySelector(`input[name="fwi-parameter"][value="${state.view.parameter}"]`).checked = true;
  updateSlider();
  showError(null);
  // A gradient shows the gradient; an inversion, its update.
  state.view.show = iterations === 0 ? "descent" : "update";
  $("#fwi-show").value = state.view.show;
  state.status = null;
  state.frame?.section?.bitmap.close();
  state.frame = null;
  $("#inversion-panel").hidden = false;
  sendView();
  // The server lists it as an inversion; its name is the configuration's.
  hooks.solver().postMessage({ type: "invert", ...configs, name: `${hooks.name()}.toml` });
  showStatus({ status: null, running: true });
}

/// The grid an attached inversion runs on, from its start.
export function attached({ start, truth }) {
  const s = readSettings(start);
  state.grid = s.grid;
  state.view.index = middle(state.view.plane);
  $("#truth").value = truth;
  $("#inversion-panel").hidden = false;
  updateSlider();
  sendView();
}

function middle(plane) {
  const n = state.grid?.points[NORMAL_AXIS[plane]] ?? 1;
  return Math.floor((n - 1) / 2);
}

function setView(change) {
  Object.assign(state.view, change);
  sendView();
}

export function sendView() {
  hooks.solver()?.postMessage({ type: "inversionView", view: { ...state.view, gain: hooks.gain() } });
}

function updateSlider() {
  if (!state.grid) return;
  const axis = NORMAL_AXIS[state.view.plane];
  const slider = $("#fwi-position");
  slider.max = state.grid.points[axis] - 1;
  slider.value = state.view.index;
  slider.title = `${["x", "y", "z"][axis]} = ${formatLength(state.view.index * state.grid.spacing[axis])}`;
}

/// The solver's report: { status, running } (status: InversionView::status, or null).
export function showStatus({ status, running }) {
  state.running = running;
  if (status) state.status = status;
  const s = state.status;
  $("#fwi-run").disabled = running;
  $("#fwi-gradient").disabled = running;
  const done = s?.stage === "done";
  const resumable = !running && s && !done;
  $("#fwi-stop").disabled = !running && !resumable;
  $("#fwi-stop").textContent = resumable ? "Continue" : "Stop";
  $("#fwi-download").disabled = !s?.unknowns;
  $("#fwi-progress").value = s ? s.done / s.units : 0;
  $("#fwi-stage").textContent = s ? describe(s, running) : running ? "Setting up" : "";
  const history = s?.history ?? [];
  $("#fwi-history").hidden = s?.starting_misfit == null;
  if (s?.starting_misfit != null) {
    const first = s.starting_misfit;
    const row = (it, misfit, step, trials, change) =>
      el("tr", {}, [it, misfit, step, trials, change].map((v) => el("td", { textContent: v })));
    $("#fwi-history tbody").replaceChildren(
      row("0", "100%", "", "", ""),
      ...history.map((h) =>
        row(
          String(h.iteration),
          `${+((100 * h.misfit) / first).toPrecision(3)}%`,
          h.step.toFixed(2) + (h.restarted ? "*" : ""),
          String(h.trials),
          `${(100 * h.largest_change).toFixed(1)}%`,
        ),
      ),
    );
  }
  drawMisfit();
}

function describe(s, running) {
  const what = {
    observing: `the true model's traces, shot ${s.shot + 1} of ${s.shots}`,
    starting: `the misfit and gradient at the start (${s.shots} shot(s))`,
    iterating: `iteration ${s.iteration} of ${s.iterations}`,
    done: s.stopped
      ? `stopped after ${s.history.length} iteration(s): ${s.stopped}`
      : s.iterations
        ? `done: ${s.history.length} iteration(s)`
        : "done: the gradient at the start",
  }[s.stage];
  const unknowns = s.unknowns ? `; ${s.unknowns.toLocaleString()} unknowns` : "";
  return s.stage === "done" ? what + unknowns : `${running ? "Computing" : "Stopped before"} ${what}${unknowns}`;
}

/// A picture from the solver: { section, clip, low, note }.
export function showFrame(f) {
  if (f.section) {
    state.frame?.section?.bitmap.close();
    state.frame = f;
    $("#fwi-picture-note").hidden = true;
  } else if (f.note) {
    $("#fwi-picture-note").textContent = f.note;
    $("#fwi-picture-note").hidden = false;
  }
  draw();
}

function draw() {
  const canvas = $("#inversion");
  const f = state.frame;
  if (!state.grid || !f?.section) return clearCanvas(canvas);
  const { points, spacing } = state.grid;
  const extents = { x: (points[0] - 1) * spacing[0], y: (points[1] - 1) * spacing[1], z: (points[2] - 1) * spacing[2] };
  const [h, v] = { xz: ["x", "z"], xy: ["x", "y"], yz: ["y", "z"] }[state.view.plane];
  const axis = NORMAL_AXIS[state.view.plane];
  const coord = { x: 0, y: 1, z: 2 };
  const near = (p) => Math.abs(p[axis] / spacing[axis] - state.view.index) <= 1;
  const markers = hooks
    .markers()
    .filter((m) => near(m.position))
    .map((m) => ({ x: m.position[coord[h]], y: m.position[coord[v]], kind: m.kind }));
  fitHeight(canvas, extents[h] / extents[v]);
  plot(canvas, f.section.bitmap, {
    x: [0, extents[h], h === "z" ? "depth" : h],
    y: [0, extents[v], v === "z" ? "depth" : v],
    aspect: extents[h] / extents[v],
    smooth: true,
    markers,
    format: formatLength,
  });
  const { parameter, show } = state.view;
  const unit = show === "descent" ? "" : ` ${UNITS[parameter]}`;
  $("#fwi-legend-min").textContent = `${formatSigned(f.low)}${unit}`;
  $("#fwi-legend-max").textContent = `${show === "model" ? "" : "+"}${formatSigned(f.clip)}${unit}`;
  $("#fwi-legend-bar").dataset.scale = show === "model" ? "sequential" : "signed";
  const what = { model: "the model", update: "the update, m − m₀", "true-update": "the true update", descent: "−gradient, the way the misfit falls fastest" }[show];
  $("#inversion-title").textContent =
    `${LABEL[parameter]}: ${what}, on the ${state.view.plane} plane at ${["x", "y", "z"][axis]} = ${formatLength(state.view.index * spacing[axis])}`;
}

function drawMisfit() {
  const s = state.status;
  const canvas = $("#misfit");
  canvas.hidden = !(s?.starting_misfit != null && s.history.length);
  if (canvas.hidden) return;
  misfitChart(canvas, [s.starting_misfit, ...s.history.map((h) => h.misfit)]);
}

/// The model the solver sent, as a file.
export function savedModel({ parameter, bytes }) {
  download(new Blob([bytes], { type: "application/octet-stream" }), `${hooks.name()}-${parameter}.npy`);
}

/// Asks the solver what the inversion would need in memory; the latest answer is shown.
export function estimateFwi() {
  if ($('[data-panel="fwi"]').hidden) return;
  let start;
  try {
    start = pair(null).start;
    showError(null);
  } catch (e) {
    $("#fwi-memory").hidden = true;
    return showError(e.message ?? String(e));
  }
  state.estimate += 1;
  hooks.solver()?.postMessage({ type: "estimate", toml: start, id: `fwi-${state.estimate}`, kind: "inversion" });
}

export function showEstimate({ id, memory, error }) {
  if (id !== `fwi-${state.estimate}`) return;
  showMemory($("#fwi-memory"), memory, error, { what: "The inversion's memory" });
}

/// Shows why the inversion can't be set up, or clears the last such message.
function showError(message) {
  const box = $("#fwi-stage");
  if (message) {
    box.textContent = message;
    box.classList.add("warn");
  } else if (box.classList.contains("warn")) {
    box.textContent = "";
    box.classList.remove("warn");
  }
}
