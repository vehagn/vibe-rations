# vibe-rations serve: the solver for remote runs, with the web page that drives it.
#
#   docker build --file Containerfile --tag vibe-rations .
#   docker run -d --name vibe-rations --restart unless-stopped \
#     -p 127.0.0.1:8088:8088 -v vibe-rations:/data vibe-rations
#
# Reach it through SSH (README, "Remote runs"). It listens on every interface of the
# container, so publish its port on the host's localhost only, as above: SSH is what
# authenticates and encrypts.
#
# Both build stages run on the build machine's platform: the WebAssembly is the same
# everywhere, and the solver is cross-compiled, so arm64 needs no emulation.

# The web page, as in web/Containerfile.
FROM --platform=$BUILDPLATFORM docker.io/library/rust:1.99.1-trixie AS page

# Must match the wasm-bindgen crate version in web/Cargo.toml exactly.
ARG WASM_BINDGEN_VERSION=0.2.130
# WASM_THREADS_TOOLCHAIN in web/build.sh must name the same one.
ARG WASM_THREADS_TOOLCHAIN=nightly-2026-11-18

SHELL ["/bin/bash", "-euo", "pipefail", "-c"]
RUN apt-get update; \
    apt-get install -y --no-install-recommends binaryen; \
    rm -rf /var/lib/apt/lists/*; \
    rustup target add wasm32-unknown-unknown; \
    rustup toolchain install "${WASM_THREADS_TOOLCHAIN}" --profile minimal \
      --component rust-src --target wasm32-unknown-unknown; \
    cargo install --locked "wasm-bindgen-cli@${WASM_BINDGEN_VERSION}"

WORKDIR /src
COPY . .
RUN WASM_THREADS_TOOLCHAIN="${WASM_THREADS_TOOLCHAIN}" web/build.sh

# The solver, for the platform the image is for.
FROM --platform=$BUILDPLATFORM docker.io/library/rust:1.99.1-trixie AS solver

ARG BUILDARCH
ARG TARGETARCH

SHELL ["/bin/bash", "-euo", "pipefail", "-c"]
RUN case "${TARGETARCH}" in \
      amd64) target=x86_64-unknown-linux-gnu; gcc=gcc-x86-64-linux-gnu ;; \
      arm64) target=aarch64-unknown-linux-gnu; gcc=gcc-aarch64-linux-gnu ;; \
      *) echo "no build for ${TARGETARCH}" >&2; exit 1 ;; \
    esac; \
    echo "${target}" > /target; \
    if [[ "${TARGETARCH}" != "${BUILDARCH}" ]]; then \
      apt-get update; \
      apt-get install -y --no-install-recommends "${gcc}" "libc6-dev-${TARGETARCH}-cross"; \
      rm -rf /var/lib/apt/lists/*; \
    fi; \
    rustup target add "${target}"

WORKDIR /src
COPY . .
RUN target="$(cat /target)"; \
    if [[ "${TARGETARCH}" != "${BUILDARCH}" ]]; then \
      export "CARGO_TARGET_$(tr a-z- A-Z_ <<< "${target}")_LINKER=${target%%-unknown*}-linux-gnu-gcc"; \
    fi; \
    cargo build --release --locked --target "${target}" --bin vibe-rations; \
    mkdir -p /out/data; \
    cp "target/${target}/release/vibe-rations" /out/

FROM docker.io/library/debian:trixie-slim

COPY --from=solver /out/vibe-rations /usr/local/bin/vibe-rations
COPY --from=page /src/web/dist /srv/page
# A volume mounted here starts out with this ownership.
COPY --from=solver --chown=65534:65534 /out/data /data
USER 65534:65534
VOLUME /data
EXPOSE 8088
ENTRYPOINT ["vibe-rations"]
CMD ["serve", "--listen", "0.0.0.0:8088", "--data", "/data", "--page", "/srv/page"]
