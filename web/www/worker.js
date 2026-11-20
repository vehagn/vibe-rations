// Runs the solver off the main thread. The page sends it a configuration and view settings;
// it steps the simulation in chunks of about CHUNK_MS and sends back a picture of the section
// and of the gather after each chunk, so the page stays responsive and pausing is immediate.

const CHUNK_MS = 50;

// Threads need SharedArrayBuffer, so a cross-origin-isolated page; otherwise one thread. The
// page can ask for a number with ?threads=N.
const asked = Number(new URL(self.location.href).searchParams.get("threads"));
const threads = self.crossOriginIsolated ? Math.max(1, asked || navigator.hardwareConcurrency || 1) : 1;
let Simulator;

let sim = null;
let running = false;
let view = { field: "vz", plane: "xz", index: 0, gain: 1, gather: null, gatherRows: 600 };
// Steps per chunk, adapted to the machine; and work done, for the speed figure.
let chunk = 1;
let work = { steps: 0, ms: 0 };

const send = (message, transfer = []) => postMessage(message, transfer);
const fail = (during, e) => send({ type: "error", during, message: String(e?.message ?? e) });

async function frame() {
  const section = sim.section(view.field, view.plane, view.index, view.gain);
  const message = {
    type: "frame",
    step: sim.step,
    totalSteps: sim.totalSteps,
    time: sim.time,
    speed: work.ms > 0 ? (sim.cells * work.steps) / (work.ms * 1e3) : 0,
    section: await bitmap(section),
    clip: section.clip,
    gather: null,
  };
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
    chunk = Math.max(1, Math.min(500, Math.round((chunk * CHUNK_MS) / Math.max(ms, 1))));
    await frame();
    // Let messages (pause, view changes) in before the next chunk.
    await new Promise((resolve) => setTimeout(resolve, 0));
  }
  if (running) {
    running = false;
    send({ type: "state", running, done: true });
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
  send({ type: "ready", threads, isolated: self.crossOriginIsolated });
} catch (e) {
  fail("init", e);
}
