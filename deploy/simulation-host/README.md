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

- **Marmousi2.** Don't run `just marmousi2` there (467 MB of SEG-Y). Send the converted
  model files from a checkout that has them, 1.4 MB:

  ```bash
  cd examples/marmousi2 && tar czf - 25m 25m-* 50m | ssh gauss 'tar xzf - -C /srv/data/vibe-rations/examples/marmousi2 && chown -R 1000:1000 /srv/data/vibe-rations/examples/marmousi2'
  ```

- **Updates.** `git pull` is a few kB (the whole history is 1.2 MB). `just remote-update gauss`
  fetches only the image's changed layers: the binary (4.4 MB) and the page (2.6 MB), and
  the Debian base (79 MB unpacked) only when CI's base image has changed. Update when there
  is something to have, not by habit; don't prune Docker's images or build cache, which would
  all come again.
- **The page.** Each load through the tunnel downloads the page again, 2.6 MB (the server
  sends `no-cache` without a validator), and watching a run keeps fetching its pictures for
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

The solver is limited by memory bandwidth, and a 2014 Xeon's is about half an M1 Max's; gauss
is for long and large runs, with 62 GiB, and frees the Mac. Building the image for x86-64-v3
(AVX2, every x86 server CPU since 2013) would give the server the 15–25%.
