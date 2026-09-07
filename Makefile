# G_CAPSUL — build, test, benchmark.
#
#   make            build both binaries into bin/
#   make test       run the self-contained test suite
#   make bench IN=reads.fq
#   make clean
#
# Requires: g++ with C++17 and OpenMP, liblzma (-llzma), python3 for the tests.

.PHONY: all test bench clean check-deps

BIN := bin

all: check-deps $(BIN)/capsule_encode $(BIN)/capsule_decode

check-deps:
	@command -v g++ >/dev/null || { echo "FATAL: g++ not found"; exit 1; }
	@echo '#include <lzma.h>' | g++ -E -x c++ - >/dev/null 2>&1 || \
	  { echo "FATAL: liblzma headers not found (apt install liblzma-dev)"; exit 1; }

# Depend on the actual sources, or make will happily hand back a stale binary.
ENC_SRC := stages/106_inprocess.cpp $(wildcard include/*.h)
DEC_SRC := stages/capsule_decode.cpp $(wildcard include/*.h)

$(BIN)/capsule_encode: $(ENC_SRC) scripts/build106.sh | $(BIN)
	bash scripts/build106.sh $@

$(BIN)/capsule_decode: $(DEC_SRC) scripts/build_decode.sh | $(BIN)
	bash scripts/build_decode.sh $@

$(BIN):
	mkdir -p $(BIN)

test: all
	bash scripts/run_tests.sh

bench: all
	@test -n "$(IN)" || { echo "usage: make bench IN=<reads.fq>"; exit 1; }
	bash scripts/benchmark_final.sh "$(IN)" $(if $(OUT),$(OUT),)

clean:
	rm -rf $(BIN)
