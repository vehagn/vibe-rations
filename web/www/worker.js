// Runs the solver off the main thread. The page sends it a configuration and view settings;
// it steps the simulation in chunks of about CHUNK_MS and sends back a picture of the section
// and of the gather after each chunk, so the page stays responsive and pausing is immediate.

const CHUNK_MS = 50;
// With live pictures off, longer chunks and progress only: more of the machine for the run.
const HEADLESS_CHUNK_MS = 500;

// Threads need SharedArrayBuffer, so a cross-origin-isolated page; otherwise one thread. The
// page can ask for a number with ?threads=N.
const asked = Number(new URL(self.location.href).searchParams.get("threads"));
const threads = self.crossOriginIsolated ? Math.max(1, asked || navigator.hardwareConcurrency || 1) : 1;
let Simulator, ModelFiles;

let sim = null;
let running = false;
let view = { field: "vz", plane: "xz", index: 0, gain: 1, gather: null, gatherRows: 600 };
let live = true;
// Steps per chunk, adapted to the machine; and work done, for the speed figure.
let chunk = 1;
let work = { steps: 0, ms: 0 };

const send = (message, transfer = []) => postMessage(message, transfer);
const fail = (during, e) => send({ type: "error", during, message: String(e?.message ?? e) });

// The model files configurations name, relative to the examples as the page's are; fetched
// once each, by URL.
const examples = new URL("examples/", self.location);
const fetched = new Map();

/// The model files `toml` names, fetched, for setting it up; free it afterwards.
async function modelFiles(toml) {
  const files = new ModelFiles();
  try {
    for (const name of Simulator.modelFiles(toml)) {
      const url = new URL(name, examples).href;
      if (!fetched.has(url)) {
        const bytes = fetch(url).then(async (r) => {
          if (!r.ok) throw new Error(`fetching the model file ${name}: ${r.status} ${r.statusText}`);
          return new Uint8Array(await r.arrayBuffer());
        });
        fetched.set(url, bytes);
        bytes.catch(() => fetched.delete(url));
      }
      files.add(name, await fetched.get(url));
    }
    return files;
  } catch (e) {
    files.free();
    throw e;
  }
}

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

// Loads so far: fetching model files lets a newer one start before an older one is done,
// and only the newest may set up.
let loads = 0;

const handlers = {
  async load({ toml }) {
    running = false;
    sim?.free();
    sim = null;
    const ticket = ++loads;
    let made;
    try {
      const files = await modelFiles(toml);
      try {
        if (ticket === loads) made = new Simulator(toml, files);
      } finally {
        files.free();
      }
    } catch (e) {
      if (ticket !== loads) return;
      send({ type: "state", running });
      return fail("load", e);
    }
    if (ticket !== loads) return made?.free();
    sim = made;
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
  // The memory a configuration would need, worked out without setting it up.
  async estimate({ toml, id }) {
    try {
      const files = await modelFiles(toml);
      try {
        send({ type: "estimate", id, memory: JSON.parse(Simulator.estimateMemory(toml, files)) });
      } finally {
        files.free();
      }
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
  ({ Simulator, ModelFiles } = pkg);
  send({ type: "ready", threads, isolated: self.crossOriginIsolated });
} catch (e) {
  fail("init", e);
}
