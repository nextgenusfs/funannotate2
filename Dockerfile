# syntax=docker/dockerfile:1.7
# funannotate2 Docker image — multi-stage build assembled with pixi
# Target: linux/amd64 (x86_64 Linux / HPC)

ARG PIXI_VERSION=0.67.0
ARG UBUNTU_VERSION=22.04
ARG PYTANTAN_VERSION=0.1.4

# ---------------------------------------------------------------------------
# Stage 1: build — resolve + install the pixi environment from pixi.lock
# ---------------------------------------------------------------------------
FROM --platform=linux/amd64 ghcr.io/prefix-dev/pixi:${PIXI_VERSION} AS build

WORKDIR /app

# Copy only what's needed for the pixi install + the local funannotate2 source
# referenced by `funannotate2 = { path = "." }` in pixi.toml.
COPY pixi.toml pixi.lock pyproject.toml setup.py README.md LICENSE ./
COPY funannotate2 ./funannotate2

# Install from the lockfile; no-op if the lockfile already matches.
RUN pixi install --locked

# Rebuild pytantan from source with AVX2 disabled (SSE4-only baseline).
# The PyPI/bioconda wheel ships with AVX2 SIMD which SIGILLs on CPUs that lack
# it — notably Rosetta 2 on Apple Silicon — and AVX2 gives no meaningful win
# for this pipeline. `pip` uses build isolation, so scikit-build-core/cython/
# scoring-matrices are pulled in transparently; we only need the C/C++ toolchain.
ARG PYTANTAN_VERSION
RUN apt-get update && \
    apt-get install -y --no-install-recommends \
        git build-essential cmake zlib1g-dev ca-certificates && \
    rm -rf /var/lib/apt/lists/*
RUN SKBUILD_CMAKE_ARGS="-DHAVE_AVX2:BOOL=OFF;-DAVX2_C_FLAGS:STRING=" \
    /app/.pixi/envs/default/bin/pip install --no-deps --force-reinstall \
        "pytantan @ git+https://github.com/althonos/pytantan.git@v${PYTANTAN_VERSION}" && \
    rm -rf /root/.cache/pip

# Pre-generate an activation script so the final image doesn't need pixi.
RUN mkdir -p /app/bin && \
    { echo '#!/bin/bash'; \
      echo 'set -e'; \
      pixi shell-hook --shell bash; \
      echo 'exec "$@"'; \
    } > /app/bin/entrypoint.sh && \
    chmod +x /app/bin/entrypoint.sh

# ---------------------------------------------------------------------------
# Stage 2: dbs — download/build the funannotate2 databases (minus BUSCO)
# ---------------------------------------------------------------------------
FROM --platform=linux/amd64 ubuntu:${UBUNTU_VERSION} AS dbs

ENV DEBIAN_FRONTEND=noninteractive \
    LANG=C.UTF-8 \
    LC_ALL=C.UTF-8

RUN apt-get update && \
    apt-get install -y --no-install-recommends ca-certificates wget gzip && \
    rm -rf /var/lib/apt/lists/*

# Pull in the pixi env built above so we can run funannotate2 install.
COPY --from=build /app/.pixi/envs/default /app/.pixi/envs/default
COPY --from=build /app/bin/entrypoint.sh /app/bin/entrypoint.sh

ENV FUNANNOTATE2_DB=/opt/funannotate2_db \
    PATH=/app/.pixi/envs/default/bin:/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin

RUN mkdir -p "${FUNANNOTATE2_DB}" && \
    funannotate2 install -d all && \
    rm -rf /root/.cache /tmp/* /var/tmp/*

# ---------------------------------------------------------------------------
# Stage 3: final — minimal runtime with env + pre-built databases
# ---------------------------------------------------------------------------
FROM --platform=linux/amd64 ubuntu:${UBUNTU_VERSION} AS final

ENV DEBIAN_FRONTEND=noninteractive \
    LANG=C.UTF-8 \
    LC_ALL=C.UTF-8

RUN apt-get update && \
    apt-get install -y --no-install-recommends \
        ca-certificates \
        tini \
        procps \
        wget && \
    rm -rf /var/lib/apt/lists/*

# Pixi env (python + bioconda tooling + funannotate2 + addons + helixerlite)
COPY --from=build /app/.pixi/envs/default /app/.pixi/envs/default
COPY --from=build /app/bin/entrypoint.sh /app/bin/entrypoint.sh

# Pre-built databases (~3 GB; BUSCO lineages download at runtime)
COPY --from=dbs /opt/funannotate2_db /opt/funannotate2_db

ENV FUNANNOTATE2_DB=/opt/funannotate2_db \
    AUGUSTUS_CONFIG_PATH=/app/.pixi/envs/default/config \
    PATH=/app/.pixi/envs/default/bin:/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin

WORKDIR /data

LABEL org.opencontainers.image.source="https://github.com/nextgenusfs/funannotate2" \
      org.opencontainers.image.description="funannotate2 genome annotation pipeline (with funannotate2-addons + helixerlite) and pre-built databases" \
      org.opencontainers.image.licenses="BSD-2-Clause"

ENTRYPOINT ["/usr/bin/tini", "--", "/app/bin/entrypoint.sh"]
CMD ["funannotate2", "--help"]

