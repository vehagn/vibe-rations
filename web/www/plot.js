// Drawing pictures with axes, and line charts, on canvases.

// Room around a plot for the tick labels: the y label sits left of them, and the last x label
// reaches past the edge.
export const MARGIN = { left: 74, right: 30, top: 10, bottom: 42 };

/// Makes a canvas as tall as an image of this aspect needs at its width, within limits, so
/// there's no empty band under a wide section.
export function fitHeight(canvas, aspect) {
  const width = canvas.getBoundingClientRect().width;
  const fit = (width - MARGIN.left - MARGIN.right) / aspect + MARGIN.top + MARGIN.bottom;
  const max = Math.min(window.innerHeight * 0.65, 680);
  canvas.style.height = `${Math.round(Math.max(180, Math.min(max, fit)))}px`;
}

export function clearCanvas(canvas) {
  const ctx = canvas.getContext("2d");
  ctx.clearRect(0, 0, canvas.width, canvas.height);
}

/// A canvas sized to its box at the screen's resolution, with its context and style.
function prepare(canvas) {
  const dpr = window.devicePixelRatio || 1;
  const box = canvas.getBoundingClientRect();
  if (box.width === 0) return null;
  canvas.width = Math.round(box.width * dpr);
  canvas.height = Math.round(box.height * dpr);
  const ctx = canvas.getContext("2d");
  ctx.setTransform(dpr, 0, 0, dpr, 0, 0);
  ctx.clearRect(0, 0, box.width, box.height);
  const style = getComputedStyle(document.documentElement);
  return { ctx, box, style, ink: style.getPropertyValue("--muted").trim() };
}

/// Draws `image` into `canvas` with axes. `x` and `y` are [first, last, label, format?]; with
/// `aspect` (width / height) the image keeps its proportions.
export function plot(canvas, image, { x, y, aspect, smooth, markers, format }) {
  const p = prepare(canvas);
  if (!p) return;
  const { ctx, box, style, ink } = p;
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

  ctx.strokeStyle = style.getPropertyValue("--line").trim();
  ctx.lineWidth = 1;
  ctx.strokeRect(left + 0.5, top + 0.5, w - 1, h - 1);
  ctx.fillStyle = ink;
  ctx.font = `12px ${style.getPropertyValue("--font").trim()}`;

  const fx = x[3] ?? format;
  const fy = y[3] ?? format;
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

/// A line of `values` against their index from 0, on a logarithmic axis of share of the first
/// value, with a dot per value: an inversion's misfit per iteration.
export function misfitChart(canvas, values) {
  const p = prepare(canvas);
  if (!p) return;
  const { ctx, box, style, ink } = p;
  const margin = { left: 74, right: 30, top: 14, bottom: 36 };
  const w = box.width - margin.left - margin.right;
  const h = box.height - margin.top - margin.bottom;
  const shares = values.map((v) => v / values[0]);
  const low = Math.min(...shares.filter((s) => s > 0), 1);
  const decades = Math.max(1, Math.ceil(-Math.log10(low) + 1e-9));
  const last = Math.max(1, values.length - 1);
  const px = (i) => margin.left + (i / last) * w;
  const py = (s) => margin.top + (-Math.log10(Math.max(s, 10 ** -decades)) / decades) * h;
  ctx.strokeStyle = style.getPropertyValue("--line").trim();
  ctx.strokeRect(margin.left + 0.5, margin.top + 0.5, w - 1, h - 1);
  ctx.fillStyle = ink;
  ctx.font = `12px ${style.getPropertyValue("--font").trim()}`;
  ctx.textAlign = "right";
  ctx.textBaseline = "middle";
  for (let d = 0; d <= decades; d++) {
    const s = 10 ** -d;
    ctx.fillRect(margin.left - 4, py(s), 4, 1);
    ctx.fillText(d === 0 ? "100%" : d <= 3 ? `${+(100 * s).toPrecision(1)}%` : `1e${2 - d}%`, margin.left - 7, py(s));
  }
  ctx.textAlign = "center";
  ctx.textBaseline = "top";
  for (let i = 0; i < values.length; i++) {
    if (values.length > 12 && i % Math.ceil(values.length / 12)) continue;
    ctx.fillRect(px(i), margin.top + h, 1, 4);
    ctx.fillText(String(i), px(i), margin.top + h + 6);
  }
  ctx.fillText("iteration", margin.left + w / 2, margin.top + h + 20);
  ctx.save();
  ctx.translate(12, margin.top + h / 2);
  ctx.rotate(-Math.PI / 2);
  ctx.fillText("misfit", 0, 0);
  ctx.restore();
  const accent = style.getPropertyValue("--accent").trim();
  ctx.strokeStyle = accent;
  ctx.fillStyle = accent;
  ctx.lineWidth = 2;
  ctx.beginPath();
  shares.forEach((s, i) => ctx[i ? "lineTo" : "moveTo"](px(i), py(s)));
  ctx.stroke();
  for (const [i, s] of shares.entries()) {
    ctx.beginPath();
    ctx.arc(px(i), py(s), 3.5, 0, 2 * Math.PI);
    ctx.fill();
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
export function ticks(lo, hi, count) {
  const span = hi - lo;
  if (!(span > 0)) return [lo];
  const raw = span / Math.min(count, 8);
  const mag = 10 ** Math.floor(Math.log10(raw));
  const step = [1, 2, 5, 10].map((m) => m * mag).find((s) => s >= raw);
  const out = [];
  for (let t = Math.ceil(lo / step) * step; t <= hi + step * 1e-9; t += step) out.push(+t.toPrecision(12));
  return out;
}
