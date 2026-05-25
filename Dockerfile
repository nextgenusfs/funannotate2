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
#
# Background: pytantan's wheels ship with AVX2 enabled, and prior to v0.1.4
# the generic path was also polluted with AVX2 flags via `add_compile_options`
# in the project-level CMakeLists.txt. Either way, on x86_64 CPUs without AVX2
# (notably Rosetta 2 on Apple Silicon) any code path that executes an AVX2
# instruction raises SIGILL. The runtime-dispatch added in v0.1.4 does not
# help us here because we want a build that is portable to *any* x86_64.
#
# pytantan's CMakeLists.txt uses FindAVX2/FindSSE4/FindNEON which auto-detect
# on the build host (GitHub Actions runners have AVX2), so we have to force
# the HAVE_* flags OFF. Pre-defining them in FindAVX2.cmake's `if((DEFINED ...))`
# guard short-circuits detection entirely.
#
# We pass these via pip's `--config-settings=cmake.define.*` rather than
# CMAKE_ARGS/SKBUILD_CMAKE_ARGS env vars, since that route goes directly into
# scikit-build-core's parser and is not subject to env-var filtering.
ARG PYTANTAN_VERSION
RUN apt-get update && \
    apt-get install -y --no-install-recommends \
        git build-essential cmake zlib1g-dev binutils ca-certificates && \
    rm -rf /var/lib/apt/lists/*
RUN /app/.pixi/envs/default/bin/pip install \
        --no-deps --no-cache-dir --force-reinstall -v \
        --config-settings=cmake.define.HAVE_AVX2=OFF \
        --config-settings=cmake.define.HAVE_SSE4=OFF \
        --config-settings=cmake.define.HAVE_NEON=OFF \
        --config-settings=cmake.define.AVX2_C_FLAGS= \
        --config-settings=cmake.define.SSE4_C_FLAGS= \
        --config-settings=cmake.define.NEON_C_FLAGS= \
        "pytantan @ git+https://github.com/althonos/pytantan.git@v${PYTANTAN_VERSION}" && \
    rm -rf /root/.cache/pip

# Verify the rebuild was effective: only `generic` should be present in the
# platform/ directory and neither lib.*.so nor generic.*.so should contain
# any AVX/AVX2 instructions (no `ymm`/`zmm` register references, no `vex`-
# encoded vector ops). Fail the build loudly otherwise so we never ship an
# image that SIGILLs at import time.
RUN set -eux; \
    PY=/app/.pixi/envs/default/bin/python; \
    PT_DIR=$("$PY" -c 'import pytantan, os; print(os.path.dirname(pytantan.__file__))'); \
    echo "pytantan installed at: $PT_DIR"; \
    ls -la "$PT_DIR/platform/"; \
    SIMD_MODS=$(ls "$PT_DIR/platform/" | grep -E '^(avx2|sse4|neon)\.' || true); \
    if [ -n "$SIMD_MODS" ]; then \
        echo "ERROR: SIMD platform modules were built despite HAVE_*=OFF: $SIMD_MODS"; \
        exit 1; \
    fi; \
    BAD=""; \
    for so in $(find "$PT_DIR" -maxdepth 2 -name '*.so'); do \
        if objdump -d -M intel --no-show-raw-insn "$so" 2>/dev/null \
            | grep -Eq '\b(ymm[0-9]+|zmm[0-9]+|vpbroadcast|vextracti128|vinserti128)\b'; then \
            echo "ERROR: $so contains AVX/AVX2 instructions"; \
            BAD="${BAD} ${so}"; \
        fi; \
    done; \
    if [ -n "$BAD" ]; then exit 1; fi; \
    "$PY" -c "import pytantan; from pytantan import Alphabet, RepeatFinder, default_scoring_matrix; print('pytantan smoke test OK', pytantan.__version__)"

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

