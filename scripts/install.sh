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

# Check if build is needed (exit early only if a fully-correct install exists).
# Requires all four rust binaries under the names PASA's PerlLib probes for
# (cdbyank_rust / faidx_rust, NOT the doubled cdbyank_rust_rust), the launcher,
# and the src/bin -> ../bin symlink. If any is missing, fall through and rebuild
# so an older/misnamed install gets repaired in place.
if [ -x "${BIN_DIR}/pasa_rust" ] && [ -x "${BIN_DIR}/slclust_rust" ] \
    && [ -x "${BIN_DIR}/cdbyank_rust" ] && [ -x "${BIN_DIR}/faidx_rust" ] \
    && [ -x "${SRC_DIR}/Launch_PASA_pipeline.pl" ] && [ -e "${SRC_DIR}/bin" ]; then
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

# Copy Rust binaries.
# Map each cargo target (release/<src>) to the exact name PASA's PerlLib probes
# for on $PATH -- do NOT blindly append "_rust", because two cargo targets are
# already named with a _rust suffix. The consumers are:
#   PASA_alignment_assembler.pm -> `which pasa_rust`
#   SingleLinkageClusterer.pm    -> _which("slclust_rust")
#   CdbTools.pm                  -> _which('cdbyank_rust'), _which('faidx_rust')
# Appending "_rust" to cdbyank_rust/faidx_rust yields cdbyank_rust_rust /
# faidx_rust_rust, which those modules never look for, so PASA silently falls
# back to the C++/non-rust path. Use explicit src:dst pairs instead.
echo "[PASA install] Installing Rust binaries..."
RUST_BINS=(
    "pasa:pasa_rust"
    "slclust:slclust_rust"
    "cdbyank_rust:cdbyank_rust"
    "faidx_rust:faidx_rust"
)
for entry in "${RUST_BINS[@]}"; do
    IFS=: read -r src_name dst_name <<< "${entry}"
    SRC_BIN="${PASA_ROOT}/pasa_rust/target/release/${src_name}"
    if [ -x "${SRC_BIN}" ]; then
        cp "${SRC_BIN}" "${BIN_DIR}/${dst_name}" || true
        # Remove any double-suffixed leftover from older installs so PASA's
        # PATH probe cannot pick up a stale/misnamed copy.
        if [ "${dst_name}" != "${src_name}_rust" ]; then
            rm -f "${BIN_DIR}/${src_name}_rust"
        fi
    else
        echo "[PASA install] WARNING: Rust binary not found: ${src_name}" >&2
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

# Copy pasa-plugins directory (includes transdecoder and other bundled tools)
if [ -d "${PASA_ROOT}/pasa-plugins" ]; then
    cp -r "${PASA_ROOT}/pasa-plugins" "${SRC_DIR}/" || true
    echo "[PASA install] Installed pasa-plugins to ${SRC_DIR}/pasa-plugins"
else
    echo "[PASA install] WARNING: pasa-plugins directory not found at ${PASA_ROOT}/pasa-plugins" >&2
fi

# Copy misc_utilities (required for minimap2 and other alignment processing)
if [ -d "${PASA_ROOT}/misc_utilities" ]; then
    cp -r "${PASA_ROOT}/misc_utilities" "${SRC_DIR}/" || true
    echo "[PASA install] Installed misc_utilities to ${SRC_DIR}/misc_utilities"
else
    echo "[PASA install] WARNING: misc_utilities directory not found at ${PASA_ROOT}/misc_utilities" >&2
fi

# Make $SRC_DIR a self-contained PASAHOME. Launch_PASA_pipeline.pl runs
# `$ENV{PATH} = "$FindBin::Bin/bin:$ENV{PATH}"`, i.e. it expects the compiled
# tools under $PASAHOME/bin, but this installer keeps them in the sibling
# $INSTALL_PREFIX/bin. Symlink src/bin -> ../bin so $PASAHOME/bin resolves to the
# rust-enabled binaries whether or not $INSTALL_PREFIX/bin is on PATH. This makes
# `export PASAHOME=$INSTALL_PREFIX/src` work standalone, without relying on any
# downstream wrapper to create the link.
if [ ! -e "${SRC_DIR}/bin" ]; then
    ln -s ../bin "${SRC_DIR}/bin"
    echo "[PASA install] Linked ${SRC_DIR}/bin -> ../bin (PASAHOME=${SRC_DIR})"
fi

# Setup TransDecoder for PASA
# TransDecoder can come from: conda package, bundled in pasa-plugins submodule, or absent
# PASA's find_tool() will search PATH first, then fall back to pasa-plugins/transdecoder/

# Ensure the transdecoder plugin directory exists (for fallback lookup)
mkdir -p "${SRC_DIR}/pasa-plugins/transdecoder"

# Priority 1: Check if conda transdecoder is installed
TRANSDECODER_UTIL="${INSTALL_PREFIX}/opt/transdecoder/util"
if [ -f "${TRANSDECODER_UTIL}/TransDecoder.LongOrfs" ]; then
    echo "[PASA install] Found conda TransDecoder at ${TRANSDECODER_UTIL}"
    # Ensure it's in bin for standard PATH discovery
    if [ ! -f "${BIN_DIR}/TransDecoder.LongOrfs" ]; then
        echo "[PASA install] Symlinking conda TransDecoder to bin..."
        ln -sf "../../opt/transdecoder/util/TransDecoder.LongOrfs" "${BIN_DIR}/TransDecoder.LongOrfs" 2>/dev/null || true
        ln -sf "../../opt/transdecoder/util/TransDecoder.Predict" "${BIN_DIR}/TransDecoder.Predict" 2>/dev/null || true
    fi
elif [ -f "${BIN_DIR}/TransDecoder.LongOrfs" ]; then
    echo "[PASA install] Found TransDecoder in bin (likely from conda/pixi PATH)"
else
    # Priority 2: Check for bundled TransDecoder in git submodule
    if [ -f "${PASA_ROOT}/pasa-plugins/transdecoder/TransDecoder.LongOrfs" ]; then
        echo "[PASA install] Found bundled TransDecoder in pasa-plugins submodule"
    else
        echo "[PASA install] WARNING: TransDecoder not found (conda package or bundled)" >&2
        echo "[PASA install]   Checked: ${TRANSDECODER_UTIL}, ${BIN_DIR}, ${PASA_ROOT}/pasa-plugins/transdecoder/" >&2
    fi
fi

# Fix conda TransDecoder symlinks: conda package creates indirection
# (bin -> opt -> util) which breaks standard PATH discovery.
# Correct bin symlinks to point directly to actual scripts.
echo "[PASA install] Validating TransDecoder bin symlinks..."
for tool in TransDecoder.LongOrfs TransDecoder.Predict; do
    bin_link="${BIN_DIR}/${tool}"
    if [ -L "${bin_link}" ]; then
        # Symlink exists; check if it points to actual executable
        if ! readlink -f "${bin_link}" | xargs test -x 2>/dev/null; then
            actual_tool="${INSTALL_PREFIX}/opt/transdecoder/util/${tool}"
            if [ -x "${actual_tool}" ]; then
                echo "[PASA install] Correcting ${tool} symlink (conda indirection) ..."
                rm "${bin_link}"
                ln -s "../../opt/transdecoder/util/${tool}" "${bin_link}"
            fi
        fi
    fi
done

# Make binaries executable
echo "[PASA install] Setting executable permissions..."
chmod +x "${BIN_DIR}"/* 2>/dev/null || true

echo "[PASA install] Installation complete at ${INSTALL_PREFIX}"

if [ ${BUILD_ERRORS} -gt 0 ]; then
    echo "[PASA install] WARNING: ${BUILD_ERRORS} build(s) had errors, but installation proceeded" >&2
    exit 1
fi

exit 0
