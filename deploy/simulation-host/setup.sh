#!/usr/bin/env bash
# Sets up a Debian 13 machine to run vibe-rations simulations, as gauss was set up
# (README.md beside this script records what that machine got and why):
#
#   - Debian's Docker, Compose, git and just;
#   - the solver server (`vibe-rations serve`, the published image) on localhost:8088, its
#     runs and results in $DATA/runs, where /etc/vibe-rations.conf tells `just remote-up` and
#     `just remote-update` to keep them;
#   - a clone of the repository in $DATA/vibe-rations with its dev container built, for
#     `just in <recipe>` there: the studies, `just marmousi2`, a build for the machine's CPU.
#
# Run it as root on the machine, or from a checkout:
#
#   ssh root@host 'bash -s' < deploy/simulation-host/setup.sh
#   ssh root@host 'DATA=/srv/data bash -s' < deploy/simulation-host/setup.sh
#
# Running it again changes nothing that is in place: it doesn't replace a running server
# (`just remote-update host` does) or touch a clone's work. It leaves users, SSH and the
# firewall alone; the server listens on localhost only, so SSH is the way in.
#
# For a metered line it fetches nothing it needn't: no source-package indices, no automatic
# updates. Marmousi2's model files go from a checkout that has them (README.md, "On a
# metered line"), not `just marmousi2` there, which downloads 467 MB.
set -euo pipefail

# The big disk: runs, results, the clone and its outputs. Docker's images and build caches
# stay on the system disk (/var/lib/docker), where speed helps compiling.
DATA=${DATA:-/srv/data}
IMAGE=${IMAGE:-ghcr.io/vehagn/vibe-rations:main}
REPO=${REPO:-https://github.com/vehagn/vibe-rations.git}
# The solver image runs as nobody (65534); the dev container as dev (1000).
SERVER_UID=65534
DEV_UID=1000

export DEBIAN_FRONTEND=noninteractive
# Source packages aren't needed, and their indices cost tens of MB an update.
[ -e /etc/apt/sources.list.orig ] || cp /etc/apt/sources.list /etc/apt/sources.list.orig
sed -i 's/^deb-src /# deb-src /' /etc/apt/sources.list
apt-get update -qq
apt-get install -y -qq docker.io docker-compose git just
systemctl enable --now docker

# The solver server.
install -d -o "$SERVER_UID" -g "$SERVER_UID" -m 0755 "$DATA/runs"
echo "data=$DATA/runs" >/etc/vibe-rations.conf
docker pull -q "$IMAGE"
docker start vibe-rations 2>/dev/null ||
    docker run -d --name vibe-rations --restart unless-stopped \
        -p 127.0.0.1:8088:8088 -v "$DATA/runs":/data "$IMAGE"

# The repository and its dev container. Its files belong to the container's user, who
# writes out/ and data/ in it; root may still run git there.
[ -d "$DATA/vibe-rations" ] || git clone -q "$REPO" "$DATA/vibe-rations"
chown -R "$DEV_UID:$DEV_UID" "$DATA/vibe-rations"
git config --global --get-all safe.directory 2>/dev/null | grep -qx "$DATA/vibe-rations" ||
    git config --global --add safe.directory "$DATA/vibe-rations"
docker compose -f "$DATA/vibe-rations/.devcontainer/compose.yaml" up --detach --build --wait

docker ps --format '{{.Names}}: {{.Status}}'
