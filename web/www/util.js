// Small helpers the page's modules share: finding elements, making them, formatting numbers.

export const $ = (selector) => document.querySelector(selector);

/// An element with properties (and dataset entries under `data`) and children.
export function el(tag, props = {}, children = []) {
  const { data, ...rest } = props;
  const e = Object.assign(document.createElement(tag), rest);
  if (data) Object.assign(e.dataset, data);
  e.append(...children);
  return e;
}

/// Bytes for people, as the solver writes them: B, KiB, MiB or GiB.
export function formatBytes(bytes) {
  const units = ["B", "KiB", "MiB", "GiB"];
  let v = bytes;
  let u = 0;
  while (v >= 1024 && u < units.length - 1) {
    v /= 1024;
    u += 1;
  }
  return u === 0 ? `${bytes} B` : `${v < 10 ? v.toFixed(2) : v.toFixed(0)} ${units[u]}`;
}

export function formatLength(m) {
  return Math.abs(m) >= 2000 ? `${+(m / 1000).toFixed(2)} km` : `${+m.toFixed(1)} m`;
}

export function formatValue(v) {
  if (!(v > 0)) return "0";
  return v >= 0.01 && v < 1e4 ? v.toPrecision(2) : v.toExponential(1);
}

/// A signed or unsigned value with three significant digits, for colour scales.
export function formatSigned(v) {
  if (v === 0) return "0";
  const a = Math.abs(v);
  const s = a >= 0.01 && a < 1e5 ? String(+a.toPrecision(3)) : a.toExponential(2);
  return v < 0 ? `−${s}` : s;
}

export function formatGain(g) {
  return g >= 1 ? String(+g.toFixed(2)) : `1/${+(1 / g).toFixed(2)}`;
}

export function formatDuration(s) {
  if (s < 60) return `${Math.ceil(s)} s`;
  return `${Math.floor(s / 60)} min ${Math.round(s % 60)} s`;
}

export function download(blob, name) {
  const a = el("a", { href: URL.createObjectURL(blob), download: name });
  a.click();
  setTimeout(() => URL.revokeObjectURL(a.href), 1000);
}

/// Fills a <details> memory box from a memory estimate: the total in its summary, the parts
/// in its list. `limit` and `where` say what it is measured against.
export function showMemory(box, memory, error, { what = "Memory" } = {}) {
  box.hidden = Boolean(error);
  if (error) return;
  const server = memory.where === "server";
  const limit = memory.limit == null ? null : formatBytes(memory.limit);
  const over = memory.limit != null && memory.peak > memory.limit;
  box.dataset.over = over;
  box.querySelector("summary").textContent = over
    ? server
      ? `${what}: ${formatBytes(memory.peak)}, more than the server has available (${limit})`
      : `${what}: ${formatBytes(memory.peak)}, more than a browser can give (${limit})`
    : server
      ? `${what}: ${formatBytes(memory.peak)}` + (limit ? ` of the server's ${limit} available` : " on the server")
      : `${what}: ${formatBytes(memory.peak)} of the browser's ${limit}`;
  box.querySelector("dl").replaceChildren(
    ...memory.parts.flatMap(({ what, bytes, transient }) => [
      el("dt", { textContent: formatBytes(bytes) }),
      el("dd", { textContent: transient ? `${what} (only then)` : what }),
    ]),
  );
}

/// Copies text to the clipboard, saying so on `button` for a moment.
export async function copy(text, button) {
  try {
    await navigator.clipboard.writeText(text);
    flash(button, "Copied");
  } catch {
    // Without clipboard access (an insecure origin), select it for the visitor to copy.
    flash(button, "Select and copy");
  }
}

function flash(button, text) {
  const before = button.textContent;
  button.textContent = text;
  setTimeout(() => (button.textContent = before), 1200);
}
