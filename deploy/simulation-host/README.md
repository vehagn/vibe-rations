# A machine for simulations

`setup.sh` makes a Debian 13 machine into a place to run vibe-rations: the solver server for
`run --remote` and the web page's server mode, and the dev container for the studies. This
page records the machine it was written for, gauss, so that it can be set up again.

## gauss

Provisioned 2026-11-22 from Debian 13.8's netinst image, reached as root over SSH with a key
(`ssh gauss`), and dedicated to simulations.

| | |
|---|---|
| CPU | Intel Xeon E5-1650 v3, 6 cores, 12 threads, 3.5 GHz (Haswell: AVX2, no AVX-512) |
| Memory | 62 GiB |
| System disk | 512 GB SSD (ADATA SX900), `/` |
| Data disk | 2 TB spinning disk (WD Red WD20EFRX), ext4, `/srv/data` |
| Packages | Debian's `main contrib non-free non-free-firmware`, with security and updates |

### What was installed and set up (2026-11-22)

All by `setup.sh` (it was first done by hand, the same steps):

- Debian's packages: `docker.io` 26.1.5+dfsg1-9+deb13u1+b1 (with `containerd`
  1.7.24~ds1-6+deb13u1 and `runc` 1.1.15+ds1-2+b4), `docker-compose` 2.26.1-4 (`docker
  compose`), `git` 2.47.3, `just` 1.40.0. Docker starts at boot.
- The solver server: the container `vibe-rations` from `ghcr.io/vehagn/vibe-rations:main`
  (first `sha256:23ac8146bd44…`, built from `a2cffdf`), restarted unless stopped, published
  on `127.0.0.1:8088` only, with `/srv/data/runs` (owned by uid 65534, the image's user) as its
  `/data`. `/etc/vibe-rations.conf` holds `data=/srv/data/runs`, which `just remote-up` and
  `just remote-update` read so that an update keeps the runs where they are.
- A clone of the repository in `/srv/data/vibe-rations`, owned by uid 1000 (the dev
  container's user; the machine has no such user), marked safe for root's git
  (`git config --global --add safe.directory`), and its dev container built and started
  (`docker compose -f .devcontainer/compose.yaml up --build`). Its cargo caches are Docker
  volumes on the SSD; its `out/` and `data/` are on the data disk.
- Pulled while measuring: `docker.io/library/rust:1.99.1-trixie` (the dev container's base
  too).
- For the metered line: the `deb-src` lines of `/etc/apt/sources.list` commented out (the
  original is `/etc/apt/sources.list.orig`); Marmousi2's converted model files (1.4 MB
  compressed) sent to the clone's `examples/marmousi2/` from the Mac; and `just in unit` run
  once in the clone before the line was metered, so the workspace's crates are in the dev
  container's cache (they come again only when `Cargo.lock` changes).

### Downloaded while unmetered (2026-11-22)

Before the metered line, `prefetch.sh` fetched what later work may need:

- `rsync` 3.5.0+ds1-0+deb13u2, `tmux` 3.5a-3, `htop` 3.4.1-5 and `time` 1.9-0.2 (no updates
  were pending; no reboot needed).
- Every crate in `Cargo.lock`, in the dev container's cache.
- Marmousi2's SEG-Y (446 MB, `data/marmousi2`) and every model file `convert.sh` makes, plus
  `examples/marmousi2/25m-smooth`, the whole section smoothed over 250 m (`whole.toml`'s
  start).
- The command line built for gauss's CPU (`-C target-cpu=native`: AVX2) as
  `/srv/data/bin/vibe-rations-native`, 15% faster than the image's. Rebuilt on 2026-11-25
  from the branch that added the run-time AVX2 and the bands (`perf/kernels`), the old one
  kept as `vibe-rations-native.before-perf`; the default build is now as fast.
- numpy 2.5.5 and matplotlib 3.11.2 in `out/venv` in the clone, for the dev container's
  Python (`out/venv/bin/python`).
- The solver image built from the clone once (`vibe-rations-solver:local`), which caches its
  toolchains, the WebAssembly tools and the Debian base for later builds.

Left alone: users (everything runs as root or in containers), SSH, the firewall (none; the
server listens on localhost only). Nothing updates by itself: `unattended-upgrades` isn't
installed and no `APT::Periodic` setting is, so `apt-daily.timer` does nothing; update by hand
(`apt-get update && apt-get upgrade`).

### Using it

```bash
just tunnel gauss                                        # the web page's server mode, http://localhost:8089
vibe-rations run examples/layered.toml --remote gauss    # a run there, its results here
just remote-update gauss                                 # the server to the latest image
ssh gauss 'cd /srv/data/vibe-rations && git pull && just in marmousi2'   # dev-container recipes there
```

`run --remote gauss` opens its own SSH tunnel; it was first used against gauss. Its traces
matched a run on the Mac to $3 \times 10^{-12}$ of the peak, not bit for bit (each machine's
maths library rounds its last bit its own way).

### On a metered line

gauss's line is metered, so what crosses it is kept small, within reason:

- **Marmousi2.** The SEG-Y and its model files are there (above). On a machine without
  them, don't run `just marmousi2` over a metered line (467 MB): send the converted model
  files from a checkout that has them, 1.4 MB:

  ```bash
  cd examples/marmousi2 && tar czf - 25m 25m-* 50m | ssh gauss 'tar xzf - -C /srv/data/vibe-rations/examples/marmousi2 && chown -R 1000:1000 /srv/data/vibe-rations/examples/marmousi2'
  ```

- **Updates.** `git pull` is a few kB (the whole history is 1.2 MB). `just remote-update gauss`
  fetches only the image's changed layers: the binary (4.4 MB) and the page (2.6 MB), and
  the Debian base (79 MB unpacked) only when CI's base image has changed. Update when there
  is something to have, not by habit; don't prune Docker's images or build cache, which would
  all come again.
- **The page.** The first load through the tunnel downloads it, 2.6 MB; later loads ask
  whether each file has changed and get an empty 304 if not (ETags in the server's
  page route), until a new image changes them. Watching a run keeps fetching its pictures for
  as long as the page is open. Follow long runs with `vibe-rations remote gauss list` (a line of JSON) and leave
  the page closed; a run nobody watches also runs about 12% faster.
- **Results.** `run --remote` fetches all of a run's outputs, snapshots included. For long
  runs leave snapshots out or few, or run them in the clone there and fetch what you need.

### Speed

`examples/layered.toml` (129 × 129 × 81 nodes, 1666 steps), Mcells/s:

| | 12 threads | 6 threads |
|---|---|---|
| gauss, the published image (baseline x86-64) | 90 | 74 |
| gauss, built for its CPU (`RUSTFLAGS="-C target-cpu=native"`: AVX2) | 104 | 93 |
| the M1 Max natively (`just mac`), 10 threads | 193 | |
| the M1 Max's dev container, 4 cores | 156 | |
| gauss, the default build since 2026-11-25 (AVX2 at run time, bands of rows) | 111 | 106 |
| the same, `VIBE_RATIONS_SIMD=off` (SSE2) | 99 | |
| the same, built for its CPU | 112 | 106 |

The first four rows are from before the changes of 2026-11-25 (the research log's "Faster
updates"); the Mac hasn't been measured since. Measured again that day with main's binaries,
the first two rows were 91 and 74, 103 and 95. Since then the default build, which
the image is, picks the AVX2 loops itself and runs as fast as one built for the CPU; on the
whole Marmousi2 section's grid it went from 85 to 114 Mcell/s (108 built for the CPU before).
The solver is still limited by memory bandwidth, and a 2014 Xeon's is about half an M1
Max's; gauss is for long and large runs, with 62 GiB, and frees the Mac.
