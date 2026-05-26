# syntax=docker/dockerfile:1.7
# funannotate2 Docker image — multi-stage build assembled with pixi
# Target: linux/amd64 (x86_64 Linux / HPC)

ARG PIXI_VERSION=0.67.0
ARG UBUNTU_VERSION=22.04

# ---------------------------------------------------------------------------
# Stage 1: build — resolve + install the pixi environment from pixi.lock
# ---------------------------------------------------------------------------
FROM --platform=linux/amd64 ghcr.io/prefix-dev/pixi:${PIXI_VERSION} AS build

WORKDIR /app

# Build tools needed by uv when pixi compiles pytantan from the git source
# (see pixi.toml). Installed before `pixi install` so the source build can run.
RUN apt-get update && \
    apt-get install -y --no-install-recommends \
        git build-essential cmake zlib1g-dev binutils ca-certificates && \
    rm -rf /var/lib/apt/lists/*

# Copy only what's needed for the pixi install + the local funannotate2 source
# referenced by `funannotate2 = { path = "." }` in pixi.toml.
COPY pixi.toml pixi.lock pyproject.toml setup.py README.md LICENSE ./
COPY funannotate2 ./funannotate2

# Force pytantan's generic-only build (no SSE4/AVX2/NEON). The fork's
# CMakeLists.txt reads PYTANTAN_DISABLE_SIMD and skips SIMD detection, so the
# resulting wheel contains only platform/generic.so. This keeps the image
# safe on Rosetta 2 (Apple Silicon) and pre-Haswell x86_64 servers.
ENV CMAKE_ARGS="-DPYTANTAN_DISABLE_SIMD=ON"

RUN pixi install --locked

# Verify pytantan was built without SIMD backends. Two independent checks:
#   (a) No avx*/sse4*/neon* platform module .so files should exist.
#   (b) No .so under pytantan may contain AVX2/AVX512 instructions (ymm/zmm
#       register refs or VEX-encoded vector mnemonics). Plain SSE2 (xmm only,
#       non-VEX) is the x86_64 baseline and is fine.
# Fail the build loudly otherwise so we never ship an image that SIGILLs at
# import time.
RUN set -eux; \
    PY=/app/.pixi/envs/default/bin/python; \
    PT_DIR=$("$PY" -c 'import pytantan, os; print(os.path.dirname(pytantan.__file__))'); \
    echo "pytantan installed at: $PT_DIR"; \
    "$PY" -c "import pytantan; print('pytantan version:', pytantan.__version__)"; \
    if [ -d "$PT_DIR/platform" ]; then ls -la "$PT_DIR/platform/"; fi; \
    if [ -d "$PT_DIR/platform" ]; then \
        SIMD_MODS=$(ls "$PT_DIR/platform/" | grep -Ei '^(avx|sse|neon)' || true); \
        if [ -n "$SIMD_MODS" ]; then \
            echo "ERROR: SIMD platform modules present despite PYTANTAN_DISABLE_SIMD=ON:"; \
            echo "$SIMD_MODS"; \
            exit 1; \
        fi; \
    fi; \
    BAD=""; \
    for so in $(find "$PT_DIR" -name '*.so'); do \
        hits=$(objdump -d -M intel --no-show-raw-insn "$so" 2>/dev/null \
            | grep -Ec '\b(ymm[0-9]+|zmm[0-9]+|vpbroadcast|vextracti128|vinserti128|vfmadd|vpermd|vpgatherdd)\b' || true); \
        echo "$so -> $hits AVX/AVX2/AVX512 hits"; \
        if [ "$hits" -gt 0 ]; then BAD="${BAD} ${so}"; fi; \
    done; \
    if [ -n "$BAD" ]; then \
        echo "ERROR: AVX-bearing .so files in pytantan:$BAD"; \
        exit 1; \
    fi; \
    "$PY" -c "import pytantan; from pytantan.lib import RepeatFinder, default_scoring_matrix; print('pytantan smoke test OK', pytantan.__version__)"

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

