#!/usr/bin/env bash
# Install PASA (Rust-enabled) into a conda/pixi environment.
#
# Usage:
#   scripts/install.sh [--install-prefix PATH]
#
# If --install-prefix is not provided, uses $CONDA_PREFIX if set.
# Otherwise defaults to /opt/pasa-rust-3.0
#
# This script:
# - Builds all PASA components (Rust, C++, plugins)
# - Installs binaries and scripts to $INSTALL_PREFIX/bin
# - Installs Perl library and pipeline scripts to $INSTALL_PREFIX/src
# - Is idempotent: safe to run multiple times

set -euo pipefail

# Determine install prefix
INSTALL_PREFIX="${CONDA_PREFIX:-/opt/pasa-rust-3.0}"

# Parse arguments
while [[ $# -gt 0 ]]; do
    case "$1" in
        --install-prefix)
            INSTALL_PREFIX="$2"
            shift 2
            ;;
        *)
            echo "Error: unknown option $1" >&2
            exit 1
            ;;
    esac
done

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PASA_ROOT="$(dirname "${SCRIPT_DIR}")"
BIN_DIR="${INSTALL_PREFIX}/bin"
SRC_DIR="${INSTALL_PREFIX}/src"

# Check if build is needed (exit early if all binaries already exist)
if [ -x "${BIN_DIR}/pasa_rust" ] && [ -x "${BIN_DIR}/slclust_rust" ] && [ -x "${SRC_DIR}/Launch_PASA_pipeline.pl" ]; then
    echo "[PASA install] Already built at ${INSTALL_PREFIX}"
    exit 0
fi

echo "[PASA install] Installing to ${INSTALL_PREFIX}"
mkdir -p "${BIN_DIR}" "${SRC_DIR}"

# Track build success
BUILD_ERRORS=0

# Build Rust components using cargo
echo "[PASA install] Building Rust components..."
if ! (cd "${PASA_ROOT}/pasa_rust" && cargo build --release); then
    echo "[PASA install] ERROR: Rust build failed" >&2
    BUILD_ERRORS=$((BUILD_ERRORS + 1))
fi

# Copy Rust binaries
echo "[PASA install] Installing Rust binaries..."
RUST_BINS=(
    pasa
    slclust
    cdbyank_rust
    faidx_rust
)
for bin in "${RUST_BINS[@]}"; do
    SRC_BIN="${PASA_ROOT}/pasa_rust/target/release/${bin}"
    if [ -x "${SRC_BIN}" ]; then
        cp "${SRC_BIN}" "${BIN_DIR}/${bin}_rust" || true
    else
        echo "[PASA install] WARNING: Rust binary not found: ${bin}" >&2
    fi
done

# Build C++ components
echo "[PASA install] Building C++ components..."
BUILD_TARGETS=(
    "pasa_cpp"
    "pasa-plugins/slclust"
    "pasa-plugins/cdbtools/cdbfasta"
    "pasa-plugins/seqclean/mdust"
    "pasa-plugins/seqclean/psx"
    "pasa-plugins/seqclean/trimpoly"
)

for target in "${BUILD_TARGETS[@]}"; do
    target_path="${PASA_ROOT}/${target}"
    if [ -d "${target_path}" ] && [ -f "${target_path}/Makefile" ]; then
        echo "[PASA install] Building ${target}..."
        if (cd "${target_path}" && make); then
            echo "[PASA install] Successfully built ${target}"
        else
            echo "[PASA install] WARNING: Failed to build ${target}" >&2
            BUILD_ERRORS=$((BUILD_ERRORS + 1))
        fi
    fi
done

# Copy C++ binaries
echo "[PASA install] Installing C++ binaries..."
CPP_BINS=(
    "pasa_cpp/pasa:pasa"
    "pasa-plugins/slclust/src/slclust:slclust"
    "pasa-plugins/cdbtools/cdbfasta/cdbfasta:cdbfasta"
    "pasa-plugins/cdbtools/cdbfasta/cdbyank:cdbyank"
    "pasa-plugins/seqclean/mdust/mdust:mdust"
    "pasa-plugins/seqclean/psx/psx:psx"
    "pasa-plugins/seqclean/trimpoly/trimpoly:trimpoly"
)

for entry in "${CPP_BINS[@]}"; do
    IFS=: read -r src_path dst_name <<< "${entry}"
    src_bin="${PASA_ROOT}/${src_path}"
    if [ -x "${src_bin}" ]; then
        cp "${src_bin}" "${BIN_DIR}/${dst_name}" || true
    else
        echo "[PASA install] WARNING: Binary not found: ${src_path}" >&2
    fi
done

# Copy seqclean utilities (shell scripts and Python)
echo "[PASA install] Installing seqclean utilities..."
SEQCLEAN_UTILS=(
    "pasa-plugins/seqclean/seqclean/seqclean"
    "pasa-plugins/seqclean/seqclean/cln2qual"
    "pasa-plugins/seqclean/seqclean/bin/seqclean.psx"
)

for util in "${SEQCLEAN_UTILS[@]}"; do
    src_util="${PASA_ROOT}/${util}"
    if [ -f "${src_util}" ]; then
        cp "${src_util}" "${BIN_DIR}/" || true
    fi
done

# Copy Perl library and pipeline orchestration
echo "[PASA install] Installing Perl libraries and pipeline..."
if [ -d "${PASA_ROOT}/PerlLib" ]; then
    cp -r "${PASA_ROOT}/PerlLib" "${SRC_DIR}/" || true
fi

# Copy pipeline scripts and config
if [ -f "${PASA_ROOT}/Launch_PASA_pipeline.pl" ]; then
    cp "${PASA_ROOT}/Launch_PASA_pipeline.pl" "${SRC_DIR}/"
fi

if [ -d "${PASA_ROOT}/pasa_conf" ]; then
    cp -r "${PASA_ROOT}/pasa_conf" "${SRC_DIR}/" || true
fi

if [ -d "${PASA_ROOT}/schema" ]; then
    cp -r "${PASA_ROOT}/schema" "${SRC_DIR}/" || true
fi

# Copy main scripts directory
if [ -d "${PASA_ROOT}/scripts" ]; then
    cp -r "${PASA_ROOT}/scripts" "${SRC_DIR}/" || true
fi

# Copy misc_utilities (required for minimap2 and other alignment processing)
if [ -d "${PASA_ROOT}/misc_utilities" ]; then
    cp -r "${PASA_ROOT}/misc_utilities" "${SRC_DIR}/" || true
    echo "[PASA install] Installed misc_utilities to ${SRC_DIR}/misc_utilities"
else
    echo "[PASA install] WARNING: misc_utilities directory not found at ${PASA_ROOT}/misc_utilities" >&2
fi

# Setup TransDecoder for PASA
# PASA expects TransDecoder in its own pasa-plugins directory
# First ensure the transdecoder plugin directory exists
mkdir -p "${SRC_DIR}/pasa-plugins/transdecoder"

# Look for TransDecoder in conda/pixi standard location: opt/transdecoder/util/
TRANSDECODER_UTIL="${INSTALL_PREFIX}/opt/transdecoder/util"
if [ -f "${TRANSDECODER_UTIL}/TransDecoder.LongOrfs" ]; then
    echo "[PASA install] Creating symlinks to pixi TransDecoder from ${TRANSDECODER_UTIL}..."
    (cd "${SRC_DIR}/pasa-plugins/transdecoder" && \
     ln -sf ../../../../../opt/transdecoder/util/TransDecoder.LongOrfs TransDecoder.LongOrfs 2>/dev/null || true && \
     ln -sf ../../../../../opt/transdecoder/util/TransDecoder.Predict TransDecoder.Predict 2>/dev/null || true)
    echo "[PASA install] Created TransDecoder symlinks at ${SRC_DIR}/pasa-plugins/transdecoder"
elif [ -f "${BIN_DIR}/TransDecoder.LongOrfs" ]; then
    # Fallback: symlink from bin if available
    (cd "${SRC_DIR}/pasa-plugins/transdecoder" && \
     ln -sf ../../../bin/TransDecoder.LongOrfs TransDecoder.LongOrfs 2>/dev/null || true && \
     ln -sf ../../../bin/TransDecoder.Predict TransDecoder.Predict 2>/dev/null || true)
    echo "[PASA install] Created TransDecoder symlinks from bin"
else
    echo "[PASA install] WARNING: TransDecoder not found at ${TRANSDECODER_UTIL} or ${BIN_DIR}" >&2
fi

# Make binaries executable
echo "[PASA install] Setting executable permissions..."
chmod +x "${BIN_DIR}"/* 2>/dev/null || true

echo "[PASA install] Installation complete at ${INSTALL_PREFIX}"

if [ ${BUILD_ERRORS} -gt 0 ]; then
    echo "[PASA install] WARNING: ${BUILD_ERRORS} build(s) had errors, but installation proceeded" >&2
    exit 1
fi

exit 0
