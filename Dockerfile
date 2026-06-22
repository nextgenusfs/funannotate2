# syntax=docker/dockerfile:1.7
# funannotate2 Docker image — multi-stage build assembled with pixi
# Target: linux/amd64 (x86_64 Linux / HPC)

ARG PIXI_VERSION=0.67.0
# Ubuntu 24.04 (noble) is required for augustus 3.5.0 from apt
# (22.04/jammy only has 3.4.0). augustus is intentionally installed via apt
# rather than bioconda — see pixi.toml and the final stage below for details.
ARG UBUNTU_VERSION=24.04

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

RUN pixi lock && pixi install --locked

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
# The post-hook PATH re-prepend keeps /opt/glimmerhmm/bin ahead of
# /app/.pixi/envs/default/bin even after `pixi shell-hook` activation
# prepends the pixi env — so our generic-baseline glimmerhmm shadows the
# bioconda v3 build at runtime. See the glimmerhmm-build stage for why.
RUN mkdir -p /app/bin && \
    { echo '#!/bin/bash'; \
      echo 'set -e'; \
      pixi shell-hook --shell bash; \
      echo 'export PATH="/opt/glimmerhmm/bin:$PATH"'; \
      echo 'exec "$@"'; \
    } > /app/bin/entrypoint.sh && \
    chmod +x /app/bin/entrypoint.sh

# ---------------------------------------------------------------------------
# Stage 1b: glimmerhmm — compile glimmerhmm from upstream source with a
# generic x86_64 baseline so it runs on every host (incl. Rosetta 2).
# The bioconda glimmerhmm 3.0.4 build sets -march=x86-64-v3 (AVX2/BMI2)
# and SIGILLs under Rosetta 2 on Apple Silicon / pre-Haswell x86_64
# (the same root cause as augustus and pytantan). glimmerhmm is NOT in
# Ubuntu apt, so we build it ourselves here.
# Install layout matches bioconda's: bin/{glimmerhmm,glimmhmm.pl,trainGlimmerHMM}
# + share/glimmerhmm/{train/*,trained_dir/*}. The trainGlimmerHMM perl
# script locates its support binaries via $RealBin/../share/glimmerhmm/train
# (the same upstream patch applied by both bioconda and this recipe).
# The bioconda glimmerhmm in the pixi env is kept as a shadowed fallback
# (and is what runs on osx-arm64 native dev) — see the PATH order in the
# final stage and the comment in pixi.toml.
# ---------------------------------------------------------------------------
FROM --platform=linux/amd64 ubuntu:${UBUNTU_VERSION} AS glimmerhmm-build

ENV DEBIAN_FRONTEND=noninteractive \
    LANG=C.UTF-8 \
    LC_ALL=C.UTF-8

RUN apt-get update && \
    apt-get install -y --no-install-recommends \
        build-essential \
        ca-certificates \
        wget \
        binutils && \
    rm -rf /var/lib/apt/lists/*

ARG GLIMMERHMM_VERSION=3.0.4
ARG GLIMMERHMM_SHA256=43e321792b9f49a3d78154cbe8ddd1fb747774dccb9e5c62fbcc37c6d0650727

WORKDIR /build
RUN wget -q https://ccb.jhu.edu/software/glimmerhmm/dl/GlimmerHMM-${GLIMMERHMM_VERSION}.tar.gz && \
    echo "${GLIMMERHMM_SHA256}  GlimmerHMM-${GLIMMERHMM_VERSION}.tar.gz" | sha256sum -c - && \
    tar -xzf GlimmerHMM-${GLIMMERHMM_VERSION}.tar.gz

WORKDIR /build/GlimmerHMM

# Same upstream fixes the bioconda recipe applies (makefile typos +
# self-locating perl entry points). Without these the `all` target in
# train/makefile names "escoreSTOP2" / "rfapp" and clean targets a
# nonexistent "trainGlimmerHMM", and the perl entry points don't find
# their support binaries when invoked via PATH.
RUN sed -i 's|^escoreSTOP2:|scoreSTOP2:|g' train/makefile && \
    sed -i 's|^rfapp:|erfapp:|g' train/makefile && \
    sed -i 's| trainGlimmerHMM||g' train/makefile && \
    sed -i 's|all:    build-icm|all:    misc.o build-icm.o build-icm-noframe.o build-icm|g' train/makefile && \
    sed -i '1 s|^.*$|#!/usr/bin/env perl|g' train/trainGlimmerHMM && \
    sed -i 's|FindBin;|FindBin qw($RealBin);|g' train/trainGlimmerHMM && \
    sed -i 's|$FindBin::Bin;|"$RealBin/../share/glimmerhmm/train";|g' train/trainGlimmerHMM && \
    sed -i '1 s|^.*$|#!/usr/bin/env perl|g' bin/glimmhmm.pl

# -O3 matches the bioconda recipe; -march=x86-64 -mtune=generic gives the
# v1 baseline that every Rosetta-2-emulated CPU can execute.
ENV CFLAGS="-O3 -march=x86-64 -mtune=generic -Wno-format -Wno-deprecated-declarations -Wno-unused-variable -Wno-unused-but-set-variable -Wno-comment" \
    CXXFLAGS="-O3 -march=x86-64 -mtune=generic -Wno-format -Wno-deprecated-declarations -Wno-unused-variable -Wno-unused-but-set-variable -Wno-comment"

RUN make -C sources CC=g++ CFLAGS="${CXXFLAGS}" -j"$(nproc)" && \
    make -C train clean && \
    make -C train all C=gcc CC=g++ CFLAGS="${CXXFLAGS}" -j"$(nproc)"

# Install into /opt/glimmerhmm with the same bin/ + share/ layout
# bioconda uses, so trainGlimmerHMM's $RealBin/../share/... lookup
# still resolves once the tree is copied into the final stage.
RUN mkdir -p /opt/glimmerhmm/bin /opt/glimmerhmm/share/glimmerhmm/train && \
    install -m 0755 bin/glimmhmm.pl sources/glimmerhmm train/trainGlimmerHMM /opt/glimmerhmm/bin/ && \
    install -m 0755 train/build-icm train/build-icm-noframe train/build1 train/build2 train/erfapp /opt/glimmerhmm/share/glimmerhmm/train/ && \
    install -m 0755 train/falsecomp train/findsites train/karlin train/score train/score2 /opt/glimmerhmm/share/glimmerhmm/train/ && \
    install -m 0755 train/scoreATG train/scoreATG2 train/scoreSTOP train/scoreSTOP2 train/splicescore /opt/glimmerhmm/share/glimmerhmm/train/ && \
    cp -f train/*.pm /opt/glimmerhmm/share/glimmerhmm/train/ && \
    cp -Rf trained_dir /opt/glimmerhmm/share/glimmerhmm/

# Verify no AVX/AVX2/AVX512 instructions slipped into the compiled
# binaries (defense-in-depth — same check pytantan uses). Plain SSE2
# (xmm only, non-VEX) is the x86_64 baseline and is fine.
RUN set -eux; \
    BAD=""; \
    for bin in /opt/glimmerhmm/bin/glimmerhmm \
               /opt/glimmerhmm/share/glimmerhmm/train/build1 \
               /opt/glimmerhmm/share/glimmerhmm/train/build2 \
               /opt/glimmerhmm/share/glimmerhmm/train/build-icm \
               /opt/glimmerhmm/share/glimmerhmm/train/build-icm-noframe; do \
        hits=$(objdump -d -M intel --no-show-raw-insn "$bin" 2>/dev/null \
            | grep -Ec '\b(ymm[0-9]+|zmm[0-9]+|vpbroadcast|vextracti128|vinserti128|vfmadd|vpermd|vpgatherdd)\b' || true); \
        echo "$bin -> $hits AVX/AVX2/AVX512 hits"; \
        if [ "$hits" -gt 0 ]; then BAD="${BAD} ${bin}"; fi; \
    done; \
    if [ -n "$BAD" ]; then \
        echo "ERROR: AVX-bearing binaries in glimmerhmm:$BAD"; \
        exit 1; \
    fi

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

# augustus + augustus-data are installed from Ubuntu apt rather than bioconda.
# The bioconda augustus 3.5.0 binary is built with a modern x86_64 microarch
# baseline (AVX2/BMI2) and SIGILLs under Rosetta 2 on Apple Silicon and on
# pre-Haswell x86_64. Ubuntu noble's augustus 3.5.0+dfsg targets a generic
# x86_64 baseline and runs on every host where this image lands.
RUN apt-get update && \
    apt-get install -y --no-install-recommends \
        ca-certificates \
        tini \
        procps \
        wget \
        augustus \
        augustus-data && \
    rm -rf /var/lib/apt/lists/*

# Pixi env (python + bioconda tooling + funannotate2 + addons + helixerlite)
COPY --from=build /app/.pixi/envs/default /app/.pixi/envs/default
COPY --from=build /app/bin/entrypoint.sh /app/bin/entrypoint.sh

# Self-compiled glimmerhmm (generic x86_64 baseline). Layered onto PATH
# ahead of /app/.pixi/envs/default/bin so it shadows the bioconda binary
# at runtime — the bioconda glimmerhmm in the pixi env is built with
# -march=x86-64-v3 (AVX2/BMI2) and SIGILLs under Rosetta 2. See the
# glimmerhmm-build stage above for the rationale; the pixi-env glimmerhmm
# remains in the image as a shadowed fallback (and is what runs natively
# on osx-arm64 outside docker).
COPY --from=glimmerhmm-build /opt/glimmerhmm /opt/glimmerhmm

# Pre-built databases (~3 GB; BUSCO lineages download at runtime)
COPY --from=dbs /opt/funannotate2_db /opt/funannotate2_db

# AUGUSTUS_CONFIG_PATH points at the apt-shipped config tree. funannotate2's
# config.py derives AUGUSTUS_BASE from dirname(AUGUSTUS_CONFIG_PATH)+"/scripts",
# which resolves to /usr/share/augustus/scripts where the apt augustus package
# installs new_species.pl, optimize_augustus.pl, etc.
ENV FUNANNOTATE2_DB=/opt/funannotate2_db \
    AUGUSTUS_CONFIG_PATH=/usr/share/augustus/config \
    PATH=/opt/glimmerhmm/bin:/app/.pixi/envs/default/bin:/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin

WORKDIR /data

LABEL org.opencontainers.image.source="https://github.com/nextgenusfs/funannotate2" \
      org.opencontainers.image.description="funannotate2 genome annotation pipeline (with funannotate2-addons + helixerlite) and pre-built databases" \
      org.opencontainers.image.licenses="BSD-2-Clause"

ENTRYPOINT ["/usr/bin/tini", "--", "/app/bin/entrypoint.sh"]
CMD ["funannotate2", "--help"]

