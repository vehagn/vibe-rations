// The page: picks a configuration, hands it to the solver in worker.js, and draws what comes
// back. All the computing happens in the worker, on this machine.

const $ = (selector) => document.querySelector(selector);

const UNITS = { vx: "m/s", vy: "m/s", vz: "m/s", p: "Pa" };
// Room around a plot for the tick labels: the y label sits left of them, and the last x label
// reaches past the edge.
const MARGIN = { left: 74, right: 30, top: 10, bottom: 42 };
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
};

// ?threads=N picks the solver's thread count, e.g. 1 to compare.
const workerUrl = new URL("./worker.js", import.meta.url);
const askedThreads = new URLSearchParams(location.search).get("threads");
if (askedThreads) workerUrl.searchParams.set("threads", askedThreads);
const worker = new Worker(workerUrl, { type: "module" });
worker.onmessage = ({ data }) => on[data.type]?.(data);
worker.onerror = (e) => showError(`The solver failed to start: ${e.message}`);

// The example list; set at the end of this module, before any message can arrive.
let examplesLoaded;

const on = {
  async ready({ threads, isolated }) {
    await examplesLoaded;
    $("#stat-threads").textContent = threads > 1 ? `${threads}` : isolated ? "1" : "1 (page not isolated)";
    setStatus("ready", "Solver ready");
    loadExample($("#example").value);
  },
  loaded({ info }) {
    state.info = info;
    clearError();
    setupView(info);
    showInfo(info);
    state.step = 0;
    state.totalSteps = info.steps;
    updateProgress({ step: 0, totalSteps: info.steps, time: 0, speed: 0 });
    setStatus("ready", "Ready to run");
    sendView();
  },
  state({ running, done }) {
    state.running = running;
    $("#run").textContent = running ? "Pause" : state.step > 0 && !done ? "Resume" : "Run";
    $("#run").disabled = !state.info || (done ?? false) || state.step >= state.totalSteps;
    if (running) setStatus("running", "Running");
    else if (done || (state.totalSteps && state.step >= state.totalSteps)) setStatus("done", "Done");
    else if (state.info) setStatus("ready", state.step > 0 ? "Paused" : "Ready to run");
  },
  frame(f) {
    state.section?.bitmap.close();
    state.section = f.section;
    if (f.gather) {
      state.gather?.bitmap.close();
      state.gather = f.gather;
    }
    state.clip = f.clip;
    state.step = f.step;
    updateProgress(f);
    draw();
  },
  npy({ component, bytes }) {
    download(new Blob([bytes], { type: "application/octet-stream" }), `${component}.npy`);
  },
  error({ during, message }) {
    const what = { load: "The configuration has a problem", run: "The run stopped" }[during];
    showError(what ? `${what}: ${message}` : message);
    if (during === "load") setStatus("error", "Configuration error");
    if (during === "run") setStatus("error", "Stopped");
  },
};

// ---------------------------------------------------------------------------- scenarios

async function loadExamples() {
  const select = $("#example");
  const files = await fetch("examples/index.json").then((r) => r.json());
  for (const file of files) {
    const option = document.createElement("option");
    option.value = file;
    option.textContent = file.replace(/\.toml$/, "").replace(/[-_]/g, " ");
    select.append(option);
  }
  const wanted = new URLSearchParams(location.search).get("example");
  select.value = files.includes(`${wanted}.toml`) ? `${wanted}.toml` : files.includes("explosion.toml") ? "explosion.toml" : files[0];
}

async function loadExample(file) {
  const toml = await fetch(`examples/${file}`).then((r) => r.text());
  $("#editor").value = toml;
  $("#example-description").textContent = describe(toml);
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

function load(toml) {
  state.info = null;
  state.section?.bitmap.close();
  state.gather?.bitmap.close();
  state.section = state.gather = null;
  $("#run").disabled = true;
  setStatus("busy", "Setting up");
  worker.postMessage({ type: "load", toml });
}

// ---------------------------------------------------------------------------- view

function setupView(info) {
  const [nx, ny, nz] = info.points;
  const source = info.sources[0]?.position ?? [0, 0, 0];
  // Through the source along the plane's normal.
  state.view.index = Math.round(source[NORMAL_AXIS[state.view.plane]] / info.spacing[NORMAL_AXIS[state.view.plane]]);
  const gather = $("#gather-component");
  gather.replaceChildren(
    ...info.components.map((c) => Object.assign(document.createElement("option"), { value: c, textContent: c })),
  );
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
});

$("#gather-component").addEventListener("change", (e) => {
  state.view.gather = e.target.value;
  sendView();
});

// ---------------------------------------------------------------------------- controls

$("#example").addEventListener("change", (e) => loadExample(e.target.value));
$("#apply").addEventListener("click", () => {
  $("#example-description").textContent = describe($("#editor").value);
  load($("#editor").value);
});
$("#reset").addEventListener("click", () => load($("#editor").value));
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
  });
}

/// Makes a canvas as tall as an image of this aspect needs at its width, within limits, so
/// there's no empty band under a wide section.
function fitHeight(canvas, aspect) {
  const width = canvas.getBoundingClientRect().width;
  const fit = (width - MARGIN.left - MARGIN.right) / aspect + MARGIN.top + MARGIN.bottom;
  const max = Math.min(window.innerHeight * 0.65, 680);
  canvas.style.height = `${Math.round(Math.max(180, Math.min(max, fit)))}px`;
}

function clearCanvas(canvas) {
  const ctx = canvas.getContext("2d");
  ctx.clearRect(0, 0, canvas.width, canvas.height);
}

/// Draws `image` into `canvas` with axes. `x` and `y` are [first, last, label, format?]; with
/// `aspect` (width / height) the image keeps its proportions.
function plot(canvas, image, { x, y, aspect, smooth, markers }) {
  const dpr = window.devicePixelRatio || 1;
  const box = canvas.getBoundingClientRect();
  if (box.width === 0) return;
  canvas.width = Math.round(box.width * dpr);
  canvas.height = Math.round(box.height * dpr);
  const ctx = canvas.getContext("2d");
  ctx.setTransform(dpr, 0, 0, dpr, 0, 0);
  ctx.clearRect(0, 0, box.width, box.height);

  const margin = MARGIN;
  let w = box.width - margin.left - margin.right;
  let h = box.height - margin.top - margin.bottom;
  if (aspect) {
    if (w / h > aspect) w = h * aspect;
    else h = w / aspect;
  }
  const left = margin.left + (box.width - margin.left - margin.right - w) / 2;
  const top = margin.top;
  ctx.imageSmoothingEnabled = smooth;
  ctx.imageSmoothingQuality = "high";
  ctx.drawImage(image, left, top, w, h);

  const style = getComputedStyle(document.documentElement);
  const ink = style.getPropertyValue("--muted").trim();
  ctx.strokeStyle = style.getPropertyValue("--line").trim();
  ctx.lineWidth = 1;
  ctx.strokeRect(left + 0.5, top + 0.5, w - 1, h - 1);
  ctx.fillStyle = ink;
  ctx.font = `12px ${style.getPropertyValue("--font").trim()}`;

  const fx = x[3] ?? formatLength;
  const fy = y[3] ?? formatLength;
  const px = (v) => left + ((v - x[0]) / (x[1] - x[0] || 1)) * w;
  const py = (v) => top + ((v - y[0]) / (y[1] - y[0] || 1)) * h;
  ctx.textAlign = "center";
  ctx.textBaseline = "top";
  for (const t of ticks(x[0], x[1], Math.max(2, Math.floor(w / 85)))) {
    ctx.fillRect(px(t), top + h, 1, 4);
    ctx.fillText(fx(t), px(t), top + h + 6);
  }
  ctx.fillText(x[2], left + w / 2, top + h + 24);
  ctx.textAlign = "right";
  ctx.textBaseline = "middle";
  for (const t of ticks(y[0], y[1], Math.max(2, Math.floor(h / 45)))) {
    ctx.fillRect(left - 4, py(t), 4, 1);
    ctx.fillText(fy(t), left - 7, py(t));
  }
  ctx.save();
  ctx.translate(12, top + h / 2);
  ctx.rotate(-Math.PI / 2);
  ctx.textAlign = "center";
  ctx.fillText(y[2], 0, 0);
  ctx.restore();

  for (const m of markers) {
    const [cx, cy] = [px(m.x), py(m.y)];
    ctx.beginPath();
    if (m.kind === "source") {
      star(ctx, cx, cy, 7, 3);
      ctx.fillStyle = "#f5b400";
    } else {
      ctx.moveTo(cx - 4, cy - 5);
      ctx.lineTo(cx + 4, cy - 5);
      ctx.lineTo(cx, cy + 2);
      ctx.closePath();
      ctx.fillStyle = "#1f8a4c";
    }
    ctx.fill();
    ctx.strokeStyle = "#ffffff";
    ctx.stroke();
  }
}

function star(ctx, x, y, outer, inner) {
  for (let i = 0; i < 10; i++) {
    const r = i % 2 ? inner : outer;
    const a = (Math.PI / 5) * i - Math.PI / 2;
    ctx[i ? "lineTo" : "moveTo"](x + r * Math.cos(a), y + r * Math.sin(a));
  }
  ctx.closePath();
}

/// About `count` round numbers between `lo` and `hi`.
function ticks(lo, hi, count) {
  const span = hi - lo;
  if (!(span > 0)) return [lo];
  const raw = span / Math.min(count, 8);
  const mag = 10 ** Math.floor(Math.log10(raw));
  const step = [1, 2, 5, 10].map((m) => m * mag).find((s) => s >= raw);
  const out = [];
  for (let t = Math.ceil(lo / step) * step; t <= hi + step * 1e-9; t += step) out.push(+t.toPrecision(12));
  return out;
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
    ["Accuracy", `${info.points_per_wavelength.toFixed(1)} points per wavelength`],
    ["Boundary", `${info.pml_width}-node layers${info.free_surface ? ", free surface" : ""}`],
    ["Receivers", `${info.receivers} × ${info.samples} samples`],
    ["Memory", `${(info.memory_bytes / 2 ** 20).toFixed(0)} MiB`],
  ];
  $("#info").replaceChildren(
    ...rows.flatMap(([k, v]) => [
      Object.assign(document.createElement("dt"), { textContent: k }),
      Object.assign(document.createElement("dd"), { textContent: v }),
    ]),
  );
  $("#warnings").replaceChildren(
    ...info.warnings.map((w) => Object.assign(document.createElement("li"), { textContent: w })),
  );
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

// ---------------------------------------------------------------------------- formatting

function formatLength(m) {
  return Math.abs(m) >= 2000 ? `${+(m / 1000).toFixed(2)} km` : `${+m.toFixed(1)} m`;
}

function formatValue(v) {
  if (!(v > 0)) return "0";
  return v >= 0.01 && v < 1e4 ? v.toPrecision(2) : v.toExponential(1);
}

function formatGain(g) {
  return g >= 1 ? String(+g.toFixed(2)) : `1/${+(1 / g).toFixed(2)}`;
}

function formatDuration(s) {
  if (s < 60) return `${Math.ceil(s)} s`;
  return `${Math.floor(s / 60)} min ${Math.round(s % 60)} s`;
}

function download(blob, name) {
  const a = Object.assign(document.createElement("a"), { href: URL.createObjectURL(blob), download: name });
  a.click();
  setTimeout(() => URL.revokeObjectURL(a.href), 1000);
}

examplesLoaded = loadExamples();
await examplesLoaded;
