cmake_minimum_required(VERSION 3.16)
project(g_capsul CXX C)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()

# Mirrors scripts/build106.sh and scripts/build_decode.sh: default to a
# portable baseline instruction set (SSE4.2/POPCNT, universal since ~2009) so
# the binary isn't tied to the exact CPU it was built on. -DCAPS_BUILD_NATIVE=ON
# opts into -march=native for a build meant to run only on this machine.
option(CAPS_BUILD_NATIVE "Build with -march=native (this machine only)" OFF)
if(MSVC)
  add_compile_options(/O2)
else()
  add_compile_options(-O3)
  if(CAPS_BUILD_NATIVE)
    add_compile_options(-march=native)
  elseif(CMAKE_SYSTEM_PROCESSOR MATCHES "^(x86_64|AMD64)$")
    # -march=x86-64-v2 is x86-only; arm64/aarch64 (e.g. Apple Silicon) has no
    # such fragmentation problem, so no arch flag is needed there.
    add_compile_options(-march=x86-64-v2)
  endif()
endif()

find_package(Threads REQUIRED)
find_package(OpenMP REQUIRED)
find_package(LibLZMA REQUIRED)

set(HTSCODECS_SRC
    thirdparty/htscodecs/fqzcomp_qual.c
    thirdparty/htscodecs/utils.c)
# g++ rejects these C sources' implicit void* conversions -- compile as C, not C++.
add_library(htscodecs_obj OBJECT ${HTSCODECS_SRC})
set_source_files_properties(${HTSCODECS_SRC} PROPERTIES LANGUAGE C)
target_include_directories(htscodecs_obj PRIVATE ${CMAKE_SOURCE_DIR}/thirdparty/htscodecs)

set(PPMD_SRC
    thirdparty/ppmd/Ppmd7.c thirdparty/ppmd/Ppmd7Enc.c
    thirdparty/ppmd/Ppmd7Dec.c thirdparty/ppmd/Alloc.c
    thirdparty/ppmd/CpuArch.c)

set(FSE_SRC
    thirdparty/fse/fse_compress.c thirdparty/fse/fse_decompress.c
    thirdparty/fse/huf_compress.c thirdparty/fse/huf_decompress.c
    thirdparty/fse/entropy_common.c thirdparty/fse/hist.c
    thirdparty/fse/debug.c)

set(COMMON_INCLUDES
    ${CMAKE_SOURCE_DIR}
    ${CMAKE_SOURCE_DIR}/include
    ${CMAKE_SOURCE_DIR}/thirdparty/ppmd
    ${CMAKE_SOURCE_DIR}/thirdparty/fse
    ${CMAKE_SOURCE_DIR}/thirdparty/htscodecs)

add_executable(best106 src/encoder.cpp ${PPMD_SRC} ${FSE_SRC} $<TARGET_OBJECTS:htscodecs_obj>)
target_include_directories(best106 PRIVATE ${COMMON_INCLUDES})
target_link_libraries(best106 PRIVATE Threads::Threads OpenMP::OpenMP_CXX LibLZMA::LibLZMA)

add_executable(capsule_decode src/decoder.cpp ${PPMD_SRC} ${FSE_SRC} $<TARGET_OBJECTS:htscodecs_obj>)
target_include_directories(capsule_decode PRIVATE ${COMMON_INCLUDES})
target_link_libraries(capsule_decode PRIVATE Threads::Threads OpenMP::OpenMP_CXX LibLZMA::LibLZMA)
