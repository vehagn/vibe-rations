// Runs the solver off the main thread. The page sends it a configuration and view settings;
// it steps the simulation in chunks of about CHUNK_MS and sends back a picture of the section
// and of the gather after each chunk, so the page stays responsive and pausing is immediate.
// An inversion goes a unit at a time instead (a shot of observed traces, the start, an
// iteration), each sending its status and a picture.

const CHUNK_MS = 50;
// With live pictures off, longer chunks and progress only: more of the machine for the run.
const HEADLESS_CHUNK_MS = 500;

// Threads need SharedArrayBuffer, so a cross-origin-isolated page; otherwise one thread. The
// page can ask for a number with ?threads=N.
const asked = Number(new URL(self.location.href).searchParams.get("threads"));
const threads = self.crossOriginIsolated ? Math.max(1, asked || navigator.hardwareConcurrency || 1) : 1;
let Simulator;
let Inverter;

let sim = null;
let running = false;
let view = { field: "vz", plane: "xz", index: 0, gain: 1, gather: null, gatherRows: 600 };
let live = true;
// An inversion of synthetic data, stepped a unit at a time (a shot, the start, an iteration)
// while `inverting`; the forward run waits meanwhile.
let inv = null;
let inverting = false;
let inversionView = { parameter: "vs", show: "update", plane: "xz", index: 0, gain: 1 };
// Steps per chunk, adapted to the machine; and work done, for the speed figure.
let chunk = 1;
let work = { steps: 0, ms: 0 };

const send = (message, transfer = []) => postMessage(message, transfer);
const fail = (during, e) => send({ type: "error", during, message: String(e?.message ?? e) });

/// The run's progress, and its pictures if they are live.
async function frame() {
  const message = {
    type: "frame",
    step: sim.step,
    totalSteps: sim.totalSteps,
    time: sim.time,
    speed: work.ms > 0 ? (sim.cells * work.steps) / (work.ms * 1e3) : 0,
    section: null,
    gather: null,
  };
  if (!live) return send(message);
  const section = sim.section(view.field, view.plane, view.index, view.gain);
  message.section = await bitmap(section);
  message.clip = section.clip;
  section.free();
  if (view.gather) {
    const g = sim.gather(view.gather, view.gatherRows);
    message.gather = await bitmap(g);
    g.free();
  }
  const transfer = [message.section.bitmap];
  if (message.gather) transfer.push(message.gather.bitmap);
  send(message, transfer);
}

async function bitmap(f) {
  const { width, height } = f;
  const pixels = new Uint8ClampedArray(f.take_rgba().buffer);
  return { width, height, bitmap: await createImageBitmap(new ImageData(pixels, width, height)) };
}

async function loop() {
  while (running && sim.step < sim.totalSteps) {
    const t0 = performance.now();
    let taken;
    try {
      taken = sim.advance(chunk);
    } catch (e) {
      running = false;
      send({ type: "state", running });
      return fail("run", e);
    }
    const ms = performance.now() - t0;
    work.steps += taken;
    work.ms += ms;
    const target = live ? CHUNK_MS : HEADLESS_CHUNK_MS;
    chunk = Math.max(1, Math.min(5000, Math.round((chunk * target) / Math.max(ms, 1))));
    await frame();
    // Let messages (pause, view changes) in before the next chunk.
    await new Promise((resolve) => setTimeout(resolve, 0));
  }
  if (running) {
    running = false;
    send({ type: "state", running, done: true });
  }
}

function sendInversion() {
  let status = null;
  try {
    status = inv ? JSON.parse(inv.status()) : null;
  } catch {}
  send({ type: "inversion", status, running: inverting });
}

/// A picture of the inversion, or why there is none (no gradient yet, say).
async function inversionFrame() {
  if (!inv) return;
  const v = inversionView;
  try {
    const f = inv.picture(v.parameter, v.show, v.plane, v.index, v.gain);
    const message = { type: "inversionFrame", section: await bitmap(f), clip: f.clip, low: f.low };
    f.free();
    send(message, [message.section.bitmap]);
  } catch (e) {
    send({ type: "inversionFrame", section: null, note: String(e?.message ?? e) });
  }
}

async function invertLoop() {
  while (inverting) {
    let more;
    try {
      // Blocks until the unit is done: a forward or gradient run per shot.
      more = inv.advance();
    } catch (e) {
      inverting = false;
      // A panic leaves the inverter unusable (its borrow is never released): report, drop it.
      fail("invert", e);
      if (e instanceof WebAssembly.RuntimeError) inv = null;
      return sendInversion();
    }
    if (!more) inverting = false;
    sendInversion();
    await inversionFrame();
    // Let messages (stop, view changes) in before the next unit.
    await new Promise((resolve) => setTimeout(resolve, 0));
  }
}

const handlers = {
  async load({ toml }) {
    running = false;
    sim?.free();
    sim = null;
    try {
      sim = new Simulator(toml);
    } catch (e) {
      send({ type: "state", running });
      return fail("load", e);
    }
    chunk = 1;
    work = { steps: 0, ms: 0 };
    send({ type: "loaded", info: JSON.parse(sim.info()) });
    send({ type: "state", running });
  },
  async run() {
    if (!sim || running) return;
    inverting = false;
    running = sim.step < sim.totalSteps;
    send({ type: "state", running });
    loop();
  },
  async pause() {
    running = false;
    send({ type: "state", running });
  },
  async view(settings) {
    view = { ...view, ...settings.view };
    if (sim && !running) {
      try {
        await frame();
      } catch (e) {
        fail("view", e);
      }
    }
  },
  // The memory a configuration would need, worked out without setting it up: as a run, or as
  // the start of an inversion.
  async estimate({ toml, id, kind }) {
    try {
      const memory = kind === "inversion" ? Inverter.estimateMemory(toml) : Simulator.estimateMemory(toml);
      send({ type: "estimate", id, memory: JSON.parse(memory) });
    } catch (e) {
      send({ type: "estimate", id, error: String(e?.message ?? e) });
    }
  },
  async live({ on }) {
    live = on;
    if (sim && !running && on) {
      try {
        await frame();
      } catch (e) {
        fail("view", e);
      }
    }
  },
  async invert({ start, truth }) {
    running = false;
    send({ type: "state", running });
    inverting = false;
    inv?.free();
    inv = null;
    try {
      inv = new Inverter(start, truth);
    } catch (e) {
      sendInversion();
      return fail("invert", e);
    }
    inverting = true;
    sendInversion();
    invertLoop();
  },
  async stopInversion() {
    // Takes effect after the unit in progress.
    inverting = false;
    sendInversion();
  },
  async resumeInversion() {
    if (!inv || inverting || JSON.parse(inv.status()).stage === "done") return;
    inverting = true;
    sendInversion();
    invertLoop();
  },
  async inversionView({ view }) {
    inversionView = { ...inversionView, ...view };
    if (!inverting) await inversionFrame();
  },
  async modelNpy({ parameter }) {
    try {
      const bytes = inv.modelNpy(parameter);
      send({ type: "modelNpy", parameter, bytes }, [bytes.buffer]);
    } catch (e) {
      fail("download", e);
    }
  },
  async npy({ component }) {
    try {
      const bytes = sim.tracesNpy(component);
      send({ type: "npy", component, bytes }, [bytes.buffer]);
    } catch (e) {
      fail("download", e);
    }
  },
};

onmessage = ({ data }) => handlers[data.type]?.(data);

try {
  const pkg = await import(threads > 1 ? "./pkg-threads/vibe_rations_web.js" : "./pkg/vibe_rations_web.js");
  await pkg.default();
  if (threads > 1) await pkg.initThreadPool(threads);
  Simulator = pkg.Simulator;
  Inverter = pkg.Inverter;
  send({ type: "ready", threads, isolated: self.crossOriginIsolated });
} catch (e) {
  fail("init", e);
}
