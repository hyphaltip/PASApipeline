# PASA Rust - Conda Recipe

This directory contains the conda recipe for building and distributing the Rust-optimized PASA (Program to Assemble Spliced Alignments) package.

## Building the Package

### Prerequisites

- conda-build (or mamba build)
- Git
- Rust toolchain
- C/C++ compiler (gcc/g++)
- Make

### Build Instructions

```bash
# Option 1: Build locally from this directory
cd PASApipeline/conda-recipe
conda build . -c conda-forge

# Option 2: Build using mamba (faster)
mamba build . -c conda-forge

# Option 3: Build with specific output directory
conda build . -c conda-forge --croot /path/to/build/output
```

### Testing the Build

```bash
# Test the built package
conda install --use-local pasa-rust

# Verify installation
Launch_PASA_pipeline.pl --help
pasa_rust --help  # Rust assembler
```

## Publishing to Bioconda

To publish this package to [bioconda](https://bioconda.github.io/), follow the bioconda contribution guide:

1. Fork the bioconda-recipes repository
2. Create a new recipe directory: `recipes/pasa-rust/`
3. Copy this recipe's contents
4. Submit a pull request

### Example bioconda recipe structure:
```
recipes/pasa-rust/
├── meta.yaml
├── build.sh
└── LICENSE.txt
```

## Environment Variables

When installed via conda, this package automatically sets:

- `PASAHOME` — Points to the PASA installation directory (`$CONDA_PREFIX/opt/pasa-rust-3.0/src`)
  - Set in `$CONDA_PREFIX/etc/conda/activate.d/pasa_rust.sh`
  - Unset in `$CONDA_PREFIX/etc/conda/deactivate.d/pasa_rust.sh`

## Available Binaries

### Rust components:
- `pasa_rust` — Rust-optimized transcript assembler
- `slclust_rust` — Rust-optimized clustering tool
- `cdbyank_rust` — Rust-optimized FASTA retrieval
- `faidx_rust` — Rust-optimized FASTA indexing

### C++ components:
- `pasa` — Original C++ assembler
- `slclust` — Original C++ clustering tool
- `cdbfasta` — FASTA indexing utility
- `cdbyank` — FASTA retrieval utility
- `mdust` — Sequence masking utility
- `psx` — Sequence utility
- `trimpoly` — Sequence trimming utility
- `seqclean` — Sequence cleaning pipeline
- `cln2qual` — Quality conversion utility
- `seqclean.psx` — Seqclean PSX configuration

### Perl pipeline:
- `Launch_PASA_pipeline.pl` — Main PASA pipeline orchestrator

## Notes

- This recipe clones from GitHub during build, so internet access is required
- Build time: ~15-30 minutes depending on system (includes Rust and C++ compilation)
- The package uses the `install.sh` script from the PASApipeline repository
- Architecture: Linux 64-bit only
- C/C++ compiler and Rust toolchain are required at build time but not at runtime
- Perl is required at runtime for the pipeline orchestration

## Development

To rebuild locally while developing:

```bash
# After making changes to install.sh or source code
conda build . -c conda-forge --no-anaconda-upload --force-rebuild
```

To use the locally built package:

```bash
# Install from local build
conda install --use-local pasa-rust

# Or temporarily:
conda install -c file://$(pwd)/../.. pasa-rust
```
