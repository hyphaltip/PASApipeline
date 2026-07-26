
## If wanting to build on mac osx using gcc instead of clang:
## run make like so:
##
##   make CC=gcc CXX=g++
##
## The C++ assembler (bin/pasa) is the default and is always built.  The
## optional Rust binaries are not built or installed unless asked for:
##
##   make WITH_RUST=1
##
## or build them on their own at any time with `make rust`.

WITH_RUST ?= 0

all:
	if [ ! -d bin ]; then mkdir bin; fi
	cd pasa_cpp && $(MAKE) && cp pasa ../bin/.
	cd pasa-plugins/slclust && $(MAKE) && cp src/slclust ../../bin/.
	cd pasa-plugins/cdbtools/cdbfasta && $(MAKE) && cp cdbfasta ../../../bin/. && cp cdbyank ../../../bin/.
	cd pasa-plugins/seqclean/mdust && $(MAKE) && cp mdust ../../../bin
	cd pasa-plugins/seqclean/psx && $(MAKE) && cp psx ../../../bin
	cd pasa-plugins/seqclean/trimpoly && $(MAKE) && cp trimpoly ../../../bin
	cp pasa-plugins/seqclean/seqclean/seqclean pasa-plugins/seqclean/seqclean/cln2qual pasa-plugins/seqclean/seqclean/bin/seqclean.psx ./bin
	@if [ "$(WITH_RUST)" = "1" ]; then \
	   $(MAKE) rust; \
	 else \
	   echo "Skipping optional Rust binaries (build them with: make WITH_RUST=1)"; \
	 fi

## Build optional Rust components (pasa assembler + slclust clusterer)
rust:
	if [ ! -d bin ]; then mkdir bin; fi
	cd pasa_rust && cargo build --release
	cp pasa_rust/target/release/pasa bin/pasa_rust
	cp pasa_rust/target/release/slclust bin/slclust_rust
	cp pasa_rust/target/release/cdbyank_rust bin/cdbyank_rust
	cp pasa_rust/target/release/faidx_rust bin/faidx_rust

## Run Rust test suite
rust-test:
	cd pasa_rust && cargo test --release

clean:
	cd pasa_cpp && $(MAKE) clean
	cd pasa-plugins/slclust && $(MAKE) clean
	cd pasa-plugins/cdbtools/cdbfasta && $(MAKE) clean
	cd pasa-plugins/seqclean/mdust && $(MAKE) clean
	cd pasa-plugins/seqclean/psx && $(MAKE) clean
	cd pasa-plugins/seqclean/trimpoly && $(MAKE) clean
	rm -f bin/*

###################################################################

