// The solver on the server that served this page (`vibe-rations serve`), behind the same
// messages as worker.js, so the page drives either the same way. Runs live on the server: they
// carry on when the page is closed, and the page can attach to them again. So do inversions,
// which are runs of their own (kind=inversion) that the server steps a unit at a time.

const POLL_MS = 200;
// With live pictures off: progress only, and seldom. The server sees nobody asking for
// pictures and runs the run headless, in longer chunks.
const HEADLESS_POLL_MS = 2000;

async function call(method, path, body) {
  const response = await fetch(path, { method, body });
  if (!response.ok) {
    const text = await response.text();
    let message = text;
    try {
      message = JSON.parse(text).error ?? text;
    } catch {}
    throw new Error(message || `${response.status} ${response.statusText}`);
  }
  return response;
}

const json = async (method, path, body) => (await call(method, path, body)).json();
const run = (id, rest = "") => `/api/runs/${encodeURIComponent(id)}${rest}`;

export class RemoteSolver {
  constructor() {
    this.onmessage = null;
    this.onerror = null;
    this.id = null;
    // Whether this page created the run and never ran it: such a run is deleted when replaced.
    this.unstarted = false;
    this.polling = false;
    this.stopped = false;
    this.view = { field: "vz", plane: "xz", index: 0, gain: 1, gather: null, gatherRows: 600 };
    this.live = true;
    // The inversion this page follows, and how it looks at it.
    this.inversion = null;
    this.following = false;
    this.inversionView = { parameter: "vs", show: "update", plane: "xz", index: 0, gain: 1 };
    this.ready();
  }

  send(data) {
    if (!this.stopped) this.onmessage?.({ data });
  }

  fail(during, e) {
    this.send({ type: "error", during, message: String(e?.message ?? e) });
  }

  postMessage(message) {
    const handler = {
      load: this.load,
      attach: this.attach,
      run: this.resume,
      pause: this.pause,
      view: this.setView,
      live: this.setLive,
      npy: this.npy,
      estimate: this.estimate,
      invert: this.invert,
      stopInversion: this.stopInversion,
      resumeInversion: this.resumeInversion,
      inversionView: this.setInversionView,
      modelNpy: this.modelNpy,
    }[message.type];
    handler?.call(this, message);
  }

  terminate() {
    this.stopped = true;
    this.forget();
  }

  async ready() {
    try {
      const info = await json("GET", "/api/info");
      this.send({ type: "ready", threads: info.threads, isolated: true, server: info });
    } catch (e) {
      this.fail("init", e);
    }
  }

  /// Deletes the run this page set up if it never ran, so trying configurations leaves nothing.
  forget() {
    if (this.id && this.unstarted) call("DELETE", run(this.id)).catch(() => {});
    this.id = null;
    this.unstarted = false;
  }

  async load({ toml, name }) {
    this.forget();
    // A newer load replaces this one: it deletes the run this one set up, or this one does,
    // and nothing is reported for it.
    const mine = (this.loads = (this.loads ?? 0) + 1);
    const replaced = () => mine !== this.loads;
    try {
      const created = await json("POST", `/api/runs?name=${encodeURIComponent(name ?? "page.toml")}`, toml);
      if (replaced()) return call("DELETE", run(created.id)).catch(() => {});
      this.id = created.id;
      this.unstarted = true;
      if (created.needs.length) {
        throw new Error(`model files (${created.needs.join(", ")}) can only be sent from the command line`);
      }
      const s = await json("POST", run(created.id, "/start?paused=true"));
      if (replaced()) return;
      this.send({ type: "loaded", info: s.info, id: s.id });
      this.send({ type: "state", running: false });
    } catch (e) {
      if (replaced()) return;
      this.send({ type: "state", running: false });
      this.fail("load", e);
    }
  }

  /// Picks up a run already on the server.
  async attach({ id }) {
    try {
      const s = await json("GET", run(id));
      if (!s.info) throw new Error("that run hasn't been set up yet");
      if (s.kind === "inversion") return await this.attachInversion(s);
      this.forget();
      this.id = id;
      const toml = await (await call("GET", run(id, "/config"))).text();
      this.send({ type: "loaded", info: s.info, id, toml, name: s.name });
      await this.report(s);
      if (s.phase === "running") this.poll();
    } catch (e) {
      this.fail("load", e);
    }
  }

  async resume() {
    if (!this.id) return;
    try {
      await json("POST", run(this.id, "/resume"));
      this.unstarted = false;
      this.send({ type: "state", running: true });
      this.poll();
    } catch (e) {
      this.fail("run", e);
    }
  }

  async pause() {
    if (!this.id) return;
    try {
      await json("POST", run(this.id, "/pause"));
    } catch {}
  }

  async setView({ view }) {
    this.view = { ...this.view, ...view };
    if (this.id && !this.polling) {
      try {
        await this.frame(await json("GET", run(this.id)));
      } catch (e) {
        this.fail("view", e);
      }
    }
  }

  async setLive({ on }) {
    this.live = on;
    if (on && this.id && !this.polling) {
      try {
        await this.frame(await json("GET", run(this.id)));
      } catch (e) {
        this.fail("view", e);
      }
    }
  }

  async npy({ component }) {
    try {
      const response = await call("GET", run(this.id, `/traces/${component}.npy`));
      this.send({ type: "npy", component, bytes: new Uint8Array(await response.arrayBuffer()) });
    } catch (e) {
      this.fail("download", e);
    }
  }

  async estimate({ toml, id, kind }) {
    try {
      const memory = await json("POST", `/api/estimate${kind ? `?kind=${kind}` : ""}`, toml);
      this.send({ type: "estimate", id, memory: { ...memory, where: "server" } });
    } catch (e) {
      this.send({ type: "estimate", id, error: String(e?.message ?? e) });
    }
  }

  /// Starts an inversion on the server: the start's configuration, then its true model.
  async invert({ start, truth, name }) {
    try {
      const created = await json("POST", `/api/runs?name=${encodeURIComponent(name ?? "inversion.toml")}&kind=inversion`, start);
      await call("PUT", run(created.id, "/files/truth.toml"), truth);
      const s = await json("POST", run(created.id, "/start"));
      this.inversion = s.id;
      this.reportInversion(s);
      this.followInversion();
    } catch (e) {
      this.send({ type: "inversion", status: null, running: false });
      this.fail("invert", e);
    }
  }

  /// Picks up an inversion already on the server: its start, its truth, where it is.
  async attachInversion(s) {
    const start = await (await call("GET", run(s.id, "/config"))).text();
    const truth = await (await call("GET", run(s.id, "/files/truth.toml"))).text();
    this.inversion = s.id;
    this.send({ type: "inversionAttached", id: s.id, start, truth, name: s.name });
    await this.reportInversion(s);
    if (s.phase === "running") this.followInversion();
  }

  async stopInversion() {
    if (!this.inversion) return;
    try {
      await json("POST", run(this.inversion, "/pause"));
    } catch {}
  }

  async resumeInversion() {
    if (!this.inversion) return;
    try {
      await json("POST", run(this.inversion, "/resume"));
      this.followInversion();
    } catch (e) {
      this.fail("invert", e);
    }
  }

  async setInversionView({ view }) {
    this.inversionView = { ...this.inversionView, ...view };
    if (this.inversion && !this.following) await this.inversionFrame();
  }

  async modelNpy({ parameter }) {
    try {
      const response = await call("GET", run(this.inversion, `/outputs/${parameter}.npy`));
      this.send({ type: "modelNpy", parameter, bytes: new Uint8Array(await response.arrayBuffer()) });
    } catch (e) {
      this.fail("download", e);
    }
  }

  /// Follows the inversion while it runs: its status and a picture after each poll.
  async followInversion() {
    if (this.following) return;
    this.following = true;
    const id = this.inversion;
    try {
      while (!this.stopped && this.inversion === id) {
        const s = await json("GET", run(id));
        const ongoing = await this.reportInversion(s);
        if (!ongoing) break;
        await new Promise((resolve) => setTimeout(resolve, 1000));
      }
    } catch (e) {
      this.fail("invert", e);
    } finally {
      this.following = false;
    }
  }

  /// Sends an inversion's status and picture; whether it is still running.
  async reportInversion(s) {
    const running = s.phase === "running";
    this.send({ type: "inversion", id: s.id, status: s.inversion, running });
    await this.inversionFrame();
    if (s.phase === "failed" || s.phase === "interrupted") this.fail("invert", s.error ?? s.phase);
    return running;
  }

  async inversionFrame() {
    const v = this.inversionView;
    const query = new URLSearchParams({ parameter: v.parameter, show: v.show, plane: v.plane, index: v.index, gain: v.gain });
    try {
      const response = await call("GET", run(this.inversion, `/model?${query}`));
      const clip = Number(response.headers.get("X-Clip"));
      const low = Number(response.headers.get("X-Low"));
      this.send({ type: "inversionFrame", section: await bitmap(response), clip, low });
    } catch (e) {
      this.send({ type: "inversionFrame", section: null, note: String(e?.message ?? e) });
    }
  }

  /// Follows the run while it runs, a picture at a time.
  async poll() {
    if (this.polling) return;
    this.polling = true;
    const id = this.id;
    try {
      while (!this.stopped && this.id === id) {
        const s = await json("GET", run(id));
        const ongoing = await this.report(s);
        if (!ongoing) break;
        await new Promise((resolve) => setTimeout(resolve, this.live ? POLL_MS : HEADLESS_POLL_MS));
      }
    } catch (e) {
      this.fail("run", e);
    } finally {
      this.polling = false;
    }
  }

  /// Sends the progress (and pictures) and the state for status `s`; whether it is still going.
  async report(s) {
    if (s.live || !this.live) await this.frame(s);
    const running = s.phase === "running";
    this.send({ type: "state", running, done: s.phase === "done" });
    if (s.phase === "failed" || s.phase === "interrupted") this.fail("run", s.error ?? s.phase);
    return running;
  }

  async frame(s) {
    const progress = { type: "frame", step: s.step, totalSteps: s.total_steps, time: s.time, speed: s.speed };
    if (!this.live) return this.send({ ...progress, section: null, gather: null });
    const v = this.view;
    const query = new URLSearchParams({ field: v.field, plane: v.plane, index: v.index, gain: v.gain });
    const section = await call("GET", run(this.id, `/section?${query}`));
    const clip = Number(section.headers.get("X-Clip"));
    const message = {
      type: "frame",
      step: s.step,
      totalSteps: s.total_steps,
      time: s.time,
      speed: s.speed,
      section: await bitmap(section),
      clip,
      gather: null,
    };
    if (v.gather) {
      const g = await call("GET", run(this.id, `/gather?component=${v.gather}&rows=${v.gatherRows}`));
      message.gather = await bitmap(g);
    }
    this.send(message);
  }
}

async function bitmap(response) {
  const image = await createImageBitmap(await response.blob());
  return { width: image.width, height: image.height, bitmap: image };
}

/// The server that served this page, if it is `vibe-rations serve`: { threads, available, … }.
export async function findServer() {
  try {
    const response = await fetch("/api/info");
    return response.ok ? await response.json() : null;
  } catch {
    return null;
  }
}

/// The runs on that server, newest first.
export const listRuns = () => json("GET", "/api/runs");
