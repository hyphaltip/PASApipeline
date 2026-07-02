#!/usr/bin/env bash
set -euo pipefail

# Use the install.sh script to build and install
mkdir -p "$PREFIX/opt/pasa-rust-3.0"
CONDA_PREFIX="$PREFIX" "$SRC_DIR/scripts/install.sh" --install-prefix "$PREFIX/opt/pasa-rust-3.0"

# Move binaries to $PREFIX/bin for easy access
mkdir -p "$PREFIX/bin"
cp "$PREFIX/opt/pasa-rust-3.0/bin"/* "$PREFIX/bin/" 2>/dev/null || true

# Create symlinks to Perl scripts in the pipeline
if [ -f "$PREFIX/opt/pasa-rust-3.0/src/Launch_PASA_pipeline.pl" ]; then
    ln -sf "$PREFIX/opt/pasa-rust-3.0/src/Launch_PASA_pipeline.pl" "$PREFIX/bin/Launch_PASA_pipeline.pl"
fi

# Set environment variable for PASA location
mkdir -p "$PREFIX/etc/conda/activate.d" "$PREFIX/etc/conda/deactivate.d"
cat > "$PREFIX/etc/conda/activate.d/pasa_rust.sh" <<'EOF'
export PASAHOME="${CONDA_PREFIX}/opt/pasa-rust-3.0/src"
EOF
cat > "$PREFIX/etc/conda/deactivate.d/pasa_rust.sh" <<'EOF'
unset PASAHOME
EOF

echo "PASA installed to $PREFIX"
echo "PASAHOME will be set to: $PREFIX/opt/pasa-rust-3.0/src"
ls -la "$PREFIX/bin" | grep -E "(pasa|slclust|cdbfasta|seqclean)" || true
