# G_CAPSUL: thin wrapper over the CMake build.
#
#   make         build build/best106 (encoder) and build/capsule_decode
#   make test    build, then run ctest
#   make clean   remove build/
#
# Requires: g++ with C++17 and OpenMP, liblzma headers, cmake, python3 (for the tests).

.PHONY: all test clean

all:
	cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
	cmake --build build -j

test: all
	cd build && ctest --output-on-failure

clean:
	rm -rf build
