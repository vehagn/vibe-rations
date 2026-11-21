// The settings form: controls that show what the configuration sets and write changes back
// into its text (config.js), which stays the source of truth. Controls name what they edit:
// data-path (a dotted path, with data-index for an entry of an array), data-source (a key of
// the chosen source), data-inversion (a key of [inversion]) and data-range (an [inversion]
// range).

import { editConfig, readSettings } from "./config.js";
import { $, el, formatLength } from "./util.js";

let current = null;
let source = 0;
let hooks = { text: () => "", change: () => {}, info: () => null };

/// Binds the controls. `text()` is the configuration being edited, `change(text)` takes the
/// edited text, and `info()` is the applied run's info (for its time step), if any.
export function initForm(h) {
  hooks = h;
  for (const input of document.querySelectorAll("[data-path]")) {
    input.addEventListener("change", () => editPath(input));
  }
  for (const input of document.querySelectorAll("[data-source]")) {
    input.addEventListener("change", () => editSource(input));
  }
  for (const input of document.querySelectorAll("[data-inversion]")) {
    input.addEventListener("change", () => editInversion(input));
  }
  for (const input of document.querySelectorAll("[data-range]")) {
    input.addEventListener("change", () => editRange(input));
  }
  for (const input of document.querySelectorAll('input[name="dt-mode"]')) {
    input.addEventListener("change", () => switchTimeStep(input.value));
  }
  $("#courant").addEventListener("change", (e) =>
    apply([
      { path: "time.courant", value: number(e.target) },
      { path: "time.dt", value: null },
    ]),
  );
  $("#dt").addEventListener("change", (e) => {
    const dt = number(e.target);
    if (dt != null) apply([{ path: "time.dt", value: dt }, { path: "time.courant", value: null }]);
  });
  $("#source-pick").addEventListener("change", (e) => {
    source = Number(e.target.value);
    fill();
  });
  $("#source-copy").addEventListener("click", () => {
    const n = current.sources.length;
    if (apply([{ path: "source", copy: source }])) source = n;
    fill();
  });
  $("#source-remove").addEventListener("click", () => {
    if (current.sources.length < 2) return;
    apply([{ path: `source.${source}`, value: null }]);
    source = Math.max(0, source - 1);
    fill();
  });
  for (const box of document.querySelectorAll("#components input")) {
    box.addEventListener("change", () => {
      const picked = [...document.querySelectorAll("#components input:checked")].map((b) => b.value);
      if (!picked.length) return fill();
      apply([{ path: "receivers.components", value: picked }]);
    });
  }
  for (const box of document.querySelectorAll("#fwi-parameters input")) {
    box.addEventListener("change", () => {
      const picked = [...document.querySelectorAll("#fwi-parameters input:checked")].map((b) => b.value);
      if (!picked.length) return fill();
      apply([...inversionSection(), { path: "inversion.parameters", value: picked }]);
    });
  }
  $("#layer-add").addEventListener("click", addLayer);
  $("#sphere-add").addEventListener("click", addSphere);
}

/// Shows what `text` sets; a configuration that doesn't parse leaves the form as it was,
/// disabled, with the reason.
export function refreshForm(text) {
  let s;
  try {
    s = readSettings(text);
  } catch (e) {
    $("#settings-error").textContent = `The form follows the configuration once it reads: ${e.message ?? e}`;
    $("#settings-error").hidden = false;
    setDisabled(true);
    return null;
  }
  $("#settings-error").hidden = true;
  setDisabled(false);
  current = s;
  source = Math.min(source, s.sources.length - 1);
  fill();
  return s;
}

/// The settings last read.
export const settings = () => current;

function setDisabled(disabled) {
  for (const e of document.querySelectorAll("#settings input, #settings select, #settings button, #fwi .grid2 input, #fwi .check-line input, #fwi-parameters input")) {
    e.disabled = disabled;
  }
}

/// Applies `edits` to the configuration; false (and the reason shown) if they fail.
function apply(edits) {
  try {
    hooks.change(editConfig(hooks.text(), edits));
    return true;
  } catch (e) {
    $("#settings-error").textContent = String(e.message ?? e);
    $("#settings-error").hidden = false;
    return false;
  }
}

/// A number input's value, or null when it is empty.
function number(input) {
  return input.value.trim() === "" ? null : Number(input.value);
}

function get(object, path) {
  return path.split(".").reduce((o, k) => o?.[k], object);
}

function editPath(input) {
  const path = input.dataset.path;
  let value =
    input.type === "checkbox" ? input.checked : input.type === "number" || input.dataset.type === "int" ? number(input) : input.value;
  if (value === null) return fill();
  if (input.dataset.index !== undefined) {
    const array = [...get(current, path)];
    array[Number(input.dataset.index)] = value;
    // One spacing for all three axes is written as one number.
    value = path === "grid.spacing" && array.every((h) => h === array[0]) ? array[0] : array;
  }
  apply([{ path, value }]);
}

function editSource(input) {
  const key = input.dataset.source;
  const path = `source.${source}.${key}`;
  if (input.dataset.index !== undefined) {
    const position = [...current.sources[source].position];
    const v = number(input);
    if (v === null) return fill();
    position[Number(input.dataset.index)] = v;
    return apply([{ path, value: position }]);
  }
  const value = input.type === "number" ? number(input) : input.value;
  if (value === null) return fill();
  apply([{ path, value }]);
}

/// An [inversion] section to write into, if there is none yet.
function inversionSection() {
  return current.inversion ? [] : [{ path: "inversion.observed", value: "observed" }];
}

function editInversion(input) {
  const key = input.dataset.inversion;
  const value = input.type === "checkbox" ? input.checked : number(input);
  // An empty fixed_above means none; other numbers can't be empty.
  if (value === null && key !== "fixed_above") return fill();
  apply([...inversionSection(), { path: `inversion.${key}`, value }]);
}

function editRange(input) {
  const key = input.dataset.range;
  const pair = [...document.querySelectorAll(`[data-range="${key}"]`)].map(number);
  const value = pair.every((v) => v !== null) ? pair : null;
  if (value === null && pair.some((v) => v !== null)) return; // wait for the other end
  apply([...inversionSection(), { path: `inversion.${key}`, value }]);
}

function switchTimeStep(mode) {
  const t = current.time;
  if (mode === "dt") {
    // Start from the step the applied run takes.
    const dt = t.dt ?? hooks.info()?.dt;
    $("#dt-field").hidden = false;
    $("#courant-field").hidden = true;
    if (dt) apply([{ path: "time.dt", value: +dt.toPrecision(4) }, { path: "time.courant", value: null }]);
  } else {
    $("#dt-field").hidden = true;
    $("#courant-field").hidden = false;
    apply([{ path: "time.dt", value: null }]);
  }
}

function addLayer() {
  const layers = current.model.layers;
  const last = layers[layers.length - 1];
  const depth = (current.grid.points[2] - 1) * current.grid.spacing[2];
  const top = +(last ? (last.top + depth) / 2 : 0).toFixed(1);
  const props = last ?? { vp: 3000, vs: 1730, rho: 2400 };
  apply([
    {
      path: `model.layer.${layers.length}`,
      value: { top, vp: round(props.vp * 1.2), vs: round(props.vs * 1.2), rho: round(props.rho * 1.05) },
    },
  ]);
}

function addSphere() {
  const { points, spacing } = current.grid;
  const extent = points.map((n, a) => (n - 1) * spacing[a]);
  const center = extent.map((e) => round(e / 2));
  const radius = round(Math.min(...extent) / 8);
  const layers = current.model.layers;
  const at = [...layers].reverse().find((l) => l.top <= center[2]) ?? layers[0];
  apply([
    {
      path: `model.sphere.${current.model.spheres.length}`,
      value: { center, radius, vp: round(at.vp * 1.1), vs: round(at.vs * 1.1), rho: at.rho },
    },
  ]);
}

const round = (v) => Math.round(v * 10) / 10;

/// Puts the settings into the controls.
function fill() {
  const s = current;
  if (!s) return;
  for (const input of document.querySelectorAll("[data-path]")) {
    let v = get(s, input.dataset.path);
    if (input.dataset.index !== undefined) v = v?.[Number(input.dataset.index)];
    setValue(input, v);
  }
  const [nx, ny, nz] = s.grid.points;
  const ext = s.grid.points.map((n, a) => formatLength((n - 1) * s.grid.spacing[a]));
  $("#grid-extent").textContent = `${ext.join(" × ")}, ${(nx * ny * nz / 1e6).toFixed(2)} M nodes`;

  // The time step: a Courant fraction (empty for the default) or dt.
  const mode = s.time.dt != null ? "dt" : "courant";
  document.querySelector(`input[name="dt-mode"][value="${mode}"]`).checked = true;
  $("#dt-field").hidden = mode !== "dt";
  $("#courant-field").hidden = mode !== "courant";
  $("#courant").placeholder = `${s.time.default_courant} (default)`;
  $("#courant").value = s.time.courant ?? "";
  $("#dt").value = s.time.dt ?? "";

  // Sources.
  $("#source-pick").replaceChildren(
    ...s.sources.map((src, i) => el("option", { value: i, textContent: `Source ${i + 1} of ${s.sources.length}` })),
  );
  $("#source-pick").value = source;
  $("#source-remove").disabled = s.sources.length < 2;
  const src = s.sources[source];
  if (src) {
    $("#source-mechanism").textContent = `${src.mechanism[0].toUpperCase()}${src.mechanism.slice(1)}; its amplitude and mechanism are in the configuration.`;
    for (const input of document.querySelectorAll("[data-source]")) {
      const key = input.dataset.source;
      const v = key === "position" ? src.position[Number(input.dataset.index)] : key === "wavelet.kind" ? src.wavelet : src.frequency;
      setValue(input, v);
    }
  }

  // Receivers.
  for (const box of document.querySelectorAll("#components input")) box.checked = s.receivers.components.includes(box.value);
  const r = s.receivers;
  $("#receiver-count").textContent =
    `${r.count} receiver(s): ${r.lines} line(s), ${r.points} point(s), set in the configuration.`;

  // The model.
  $("#model-files").hidden = !s.model.files;
  for (const id of ["#layers", "#spheres", "#layer-add", "#sphere-add"]) $(id).hidden = s.model.files;
  fillRows(
    "#layers",
    s.model.layers,
    (l, i) => [
      ["top", l.top],
      ["vp", l.vp],
      ["vs", l.vs],
      ["rho", l.rho],
    ].map(([k, v]) => [`model.layer.${i}.${k}`, v]),
    (i) => (s.model.layers.length > 1 ? `model.layer.${i}` : null),
  );
  fillRows(
    "#spheres",
    s.model.spheres,
    (sp, i) => [
      ...sp.center.map((c, a) => [`model.sphere.${i}.center`, c, a]),
      [`model.sphere.${i}.radius`, sp.radius],
      [`model.sphere.${i}.vp`, sp.vp],
      [`model.sphere.${i}.vs`, sp.vs],
      [`model.sphere.${i}.rho`, sp.rho],
    ],
    (i) => `model.sphere.${i}`,
  );
  $("#spheres").hidden = s.model.files || !s.model.spheres.length;

  // [inversion], or what it would be.
  const inv = s.inversion ?? s.inversion_defaults;
  for (const box of document.querySelectorAll("#fwi-parameters input")) box.checked = inv.parameters.includes(box.value);
  for (const input of document.querySelectorAll("[data-inversion]")) setValue(input, inv[input.dataset.inversion]);
  for (const input of document.querySelectorAll("[data-range]")) setValue(input, inv[input.dataset.range]?.[Number(input.dataset.index)]);
}

/// Table rows of number inputs, one per entry, with a remove button where `removable(i)`
/// gives the path to remove.
function fillRows(table, entries, cells, removable) {
  const body = $(`${table} tbody`);
  body.replaceChildren(
    ...entries.map((entry, i) => {
      const inputs = cells(entry, i).map(([path, v, index]) => {
        const input = el("input", { type: "number", step: "any", value: v });
        input.addEventListener("change", () => {
          const value = number(input);
          if (value === null) return fill();
          if (index === undefined) return apply([{ path, value }]);
          // One coordinate of a sphere's centre.
          const center = [...entry.center];
          center[index] = value;
          apply([{ path, value: center }]);
        });
        return el("td", {}, [input]);
      });
      const path = removable(i);
      const remove = el("button", { textContent: "×", title: "Remove", disabled: !path });
      remove.addEventListener("click", () => path && apply([{ path, value: null }]));
      return el("tr", {}, [...inputs, el("td", {}, [remove])]);
    }),
  );
}

function setValue(input, v) {
  // Leave alone what the visitor is typing.
  if (input === document.activeElement && input.type !== "checkbox") return;
  if (input.type === "checkbox") input.checked = Boolean(v);
  else input.value = v ?? "";
}
