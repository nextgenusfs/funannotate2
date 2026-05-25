# syntax=docker/dockerfile:1.7
# funannotate2 Docker image — multi-stage build assembled with pixi
# Target: linux/amd64 (x86_64 Linux / HPC)

ARG PIXI_VERSION=0.67.0
ARG UBUNTU_VERSION=22.04
ARG PYTANTAN_VERSION=0.1.3

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
# Background: pytantan's PyPI wheels ship with AVX2 SIMD enabled, which SIGILLs
# on x86_64 hosts that don't translate every AVX2 opcode (Rosetta 2 on Apple
# Silicon has known gaps, and some pre-Haswell servers lack AVX2 outright).
#
# We pin to v0.1.3 and rebuild from a patched checkout. Two independent gates
# guarantee the resulting .so files are free of AVX2/AVX512:
#   (1) CMakeLists.txt is patched to set HAVE_SSE4/HAVE_AVX2/HAVE_NEON to OFF
#       and drop the include() lines for the corresponding Find*.cmake modules,
#       so SIMD detection never runs and add_compile_options(-mavx2) cannot
#       fire even if the env-var plumbing is wrong.
#   (2) The build is invoked with scikit-build-core's real config-settings
#       interface (cmake.define.*) and CMAKE_ARGS, both of which are honored
#       (the previously-used SKBUILD_CMAKE_ARGS env var is NOT a real knob).
ARG PYTANTAN_VERSION
RUN apt-get update && \
    apt-get install -y --no-install-recommends \
        git build-essential cmake zlib1g-dev binutils ca-certificates && \
    rm -rf /var/lib/apt/lists/*

# CACHEBUST forces the pytantan rebuild RUN below to re-execute when bumped,
# bypassing BuildKit's GHA layer cache. Bump on any pytantan-related change.
ARG PYTANTAN_CACHEBUST=2

# Stage-local ENVs so the values are available inside the quoted heredoc below
# without subjecting the Python source to Dockerfile-level ${...} expansion.
ENV _PT_VERSION=${PYTANTAN_VERSION} \
    _PT_CACHEBUST=${PYTANTAN_CACHEBUST}

RUN --mount=type=tmpfs,target=/tmp/build <<'BASH'
set -eux
echo "pytantan rebuild (cachebust=${_PT_CACHEBUST}, version=${_PT_VERSION})"
git clone --depth 1 --branch "v${_PT_VERSION}" \
    https://github.com/althonos/pytantan.git /tmp/build/pytantan
cd /tmp/build/pytantan

python3 <<'PY'
import pathlib, re

def strip_if_block(text, var, replacement):
    """Remove a top-level if(<var>) ... endif() block, handling nested if()."""
    lines = text.splitlines(keepends=True)
    out = []
    i = 0
    open_re = re.compile(r'^\s*if\s*\(')
    close_re = re.compile(r'^\s*endif\s*\(')
    target_re = re.compile(rf'^\s*if\s*\(\s*{re.escape(var)}\s*\)')
    while i < len(lines):
        if target_re.match(lines[i]):
            depth = 1
            i += 1
            while i < len(lines) and depth > 0:
                if open_re.match(lines[i]):
                    depth += 1
                elif close_re.match(lines[i]):
                    depth -= 1
                i += 1
            if replacement:
                out.append(replacement)
        else:
            out.append(lines[i])
            i += 1
    return "".join(out)

for path, repl_factory in [
    ("CMakeLists.txt", lambda v: f"set({v} OFF)\n"),
    ("src/pytantan/platform/CMakeLists.txt", lambda v: ""),
]:
    p = pathlib.Path(path)
    src = p.read_text()
    if path == "CMakeLists.txt":
        src = re.sub(
            r'include\("src/scripts/cmake/Find(SSE4|AVX2|NEON)\.cmake"\)\n',
            "", src,
        )
    for var in ("HAVE_SSE4", "HAVE_AVX2", "HAVE_NEON"):
        src = strip_if_block(src, var, repl_factory(var))
    p.write_text(src)
    print(f"=== patched {path} ===")
    print(src)
PY

CMAKE_ARGS="-DHAVE_SSE4:BOOL=OFF -DHAVE_AVX2:BOOL=OFF -DHAVE_NEON:BOOL=OFF" \
/app/.pixi/envs/default/bin/pip install -v \
    --no-deps --no-cache-dir --force-reinstall \
    --config-settings=cmake.define.HAVE_SSE4=OFF \
    --config-settings=cmake.define.HAVE_AVX2=OFF \
    --config-settings=cmake.define.HAVE_NEON=OFF \
    /tmp/build/pytantan
rm -rf /root/.cache/pip
BASH

# Verify the rebuild was effective. Two independent checks:
#   (a) No avx*/sse4*/neon* platform module .so files should exist.
#   (b) No .so under the pixi env may contain AVX2/AVX512 instructions
#       (ymm/zmm register refs or VEX-encoded vector mnemonics). Anything
#       SSE2 (xmm only, non-VEX) is safe — that's the x86_64 baseline.
# Fail the build loudly otherwise so we can never ship an image that SIGILLs
# at import time.
RUN set -eux; \
    PY=/app/.pixi/envs/default/bin/python; \
    SITE=$("$PY" -c 'import sysconfig; print(sysconfig.get_paths()["purelib"])'); \
    PT_DIR=$("$PY" -c 'import pytantan, os; print(os.path.dirname(pytantan.__file__))'); \
    echo "pytantan installed at: $PT_DIR"; \
    "$PY" -c "import pytantan; print('pytantan version:', pytantan.__version__)"; \
    if [ -d "$PT_DIR/platform" ]; then ls -la "$PT_DIR/platform/"; fi; \
    if [ -d "$PT_DIR/platform" ]; then \
        SIMD_MODS=$(ls "$PT_DIR/platform/" | grep -Ei '^(avx|sse|neon)' || true); \
        if [ -n "$SIMD_MODS" ]; then \
            echo "ERROR: SIMD platform modules present despite HAVE_*=OFF:"; \
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

