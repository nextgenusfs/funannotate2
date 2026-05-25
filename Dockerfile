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

# Rebuild pytantan from source with all SIMD backends disabled.
# Background: the PyPI/bioconda wheel for 0.1.4 ships with AVX2 SIMD baked
# into pytantan/platform/avx2.*.so, and pytantan's lib module imports that
# platform module unconditionally at `import pytantan` time when the
# extension was built with AVX2 support. That SIGILLs on x86_64 hosts that
# lack AVX2 -- notably Rosetta 2 on Apple Silicon, which only emulates up
# to SSE4.2. Building without AVX2/SSE4/NEON drops to pytantan's pure C++
# `generic` backend on all hosts, which is the only safe option for an
# image that runs under Rosetta 2.
#
# Notes:
#   - `--config-settings=cmake.define.HAVE_*=OFF` is the supported
#     scikit-build-core channel for forwarding cmake -D options. The
#     `CMAKE_ARGS` env-var is filtered inconsistently and was the reason
#     PRs #64/#65 published images that still SIGILLed.
#   - pip uses build isolation, so cython/scikit-build-core are pulled in
#     transparently; we only need the C/C++ toolchain and binutils
#     (for the verification step below).
ARG PYTANTAN_VERSION
RUN apt-get update && \
    apt-get install -y --no-install-recommends \
        git build-essential cmake zlib1g-dev binutils ca-certificates && \
    rm -rf /var/lib/apt/lists/*
RUN /app/.pixi/envs/default/bin/pip install \
        --no-cache-dir --no-deps --force-reinstall -v \
        --config-settings=cmake.define.HAVE_AVX2=OFF \
        --config-settings=cmake.define.HAVE_SSE4=OFF \
        --config-settings=cmake.define.HAVE_NEON=OFF \
        --config-settings=cmake.define.AVX2_C_FLAGS= \
        --config-settings=cmake.define.SSE4_C_FLAGS= \
        --config-settings=cmake.define.NEON_C_FLAGS= \
        "pytantan @ git+https://github.com/althonos/pytantan.git@v${PYTANTAN_VERSION}"

# Fail the build loudly if any SIMD-tainted code survived the rebuild.
# This protects against silent regressions (e.g. scikit-build-core changing
# how it parses --config-settings, or a future pytantan adding a new
# detection path).
RUN PT_DIR="$(/app/.pixi/envs/default/bin/python -c 'import os, pytantan; print(os.path.dirname(pytantan.__file__))')" && \
    echo "pytantan install dir: $PT_DIR" && \
    ls -la "$PT_DIR/platform/" && \
    if ls "$PT_DIR"/platform/avx2*.so "$PT_DIR"/platform/sse4*.so "$PT_DIR"/platform/neon*.so 2>/dev/null | grep -q .; then \
        echo "FATAL: pytantan/platform/ still contains SIMD-specific extension modules" >&2; \
        exit 1; \
    fi && \
    for so in $(find "$PT_DIR" -name '*.so'); do \
        echo "objdump scan: $so"; \
        if objdump -d -M intel "$so" 2>/dev/null \
             | grep -Eo '\b(ymm[0-9]+|zmm[0-9]+|vpbroadcast[a-z]*|vextracti128|vinserti128)\b' \
             | sort -u | head; then \
            echo "FATAL: AVX2/AVX-512 mnemonics found in $so" >&2; \
            exit 1; \
        fi; \
    done && \
    /app/.pixi/envs/default/bin/python -c \
        "from pytantan import Alphabet, RepeatFinder, default_scoring_matrix; print('pytantan import OK')" && \
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

