
## If wanting to build on mac osx using gcc instead of clang:
## run make like so:
##
##   make CC=gcc CXX=g++

all:
	if [ ! -d bin ]; then mkdir bin; fi
	cd pasa_cpp && $(MAKE) && cp pasa ../bin/.
	cd pasa-plugins/slclust && $(MAKE) && cp src/slclust ../../bin/.
	cd pasa-plugins/cdbtools/cdbfasta && $(MAKE) && cp cdbfasta ../../../bin/. && cp cdbyank ../../../bin/.
	cd pasa-plugins/seqclean/mdust && $(MAKE) && cp mdust ../../../bin
	cd pasa-plugins/seqclean/psx && $(MAKE) && cp psx ../../../bin
	cd pasa-plugins/seqclean/trimpoly && $(MAKE) && cp trimpoly ../../../bin
	cp pasa-plugins/seqclean/seqclean/seqclean pasa-plugins/seqclean/seqclean/cln2qual pasa-plugins/seqclean/seqclean/bin/seqclean.psx ./bin
	$(MAKE) rust

## Build optimized Rust components (pasa assembler + slclust clusterer)
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

