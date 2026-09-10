#!/usr/bin/env bash
# Patch SPRING for GCC 13: insert <cstdint> into every COMPILED source/header.
#
# Without this the build fails with:
#     'uint8_t' does not name a type
# cascading into dozens of "no member named" errors, because GCC 13 no longer
# transitively includes <cstdint>. old_src/ is not built and is skipped.
#
#   usage: patch_spring.sh [/root/SPRING]
set -eu
SPRING_DIR="${1:-/root/SPRING}"
[ -d "$SPRING_DIR/src" ] || { echo "no such SPRING source dir: $SPRING_DIR/src" >&2; exit 1; }
cd "$SPRING_DIR/src"
#
for f in bitset_util.cpp id_compression/src/stream_model.cpp params.h spring.h \
    id_compression/src/sam_file_allocation.cpp encoder.cpp reorder_compress_streams.cpp \
    id_compression/include/id_compression.h BooPHF.h decompress.cpp main.cpp libbsc/bsc.h \
    id_compression/src/compression.cpp encoder.h pe_encode.cpp bitset_util.h reorder.h \
    util.cpp spring.cpp decompress.h reorder_compress_quality_id.cpp preprocess.cpp \
    id_compression/src/io_functions.cpp id_compression/include/sam_block.h \
    id_compression/include/Arithmetic_stream.h id_compression/src/sam_models.cpp \
    reorder_compress_quality_id.h id_compression/include/stream_model.h \
    id_compression/src/id_compression.cpp libbsc/bsc_str_array.cpp \
    id_compression/src/Arithmetic_stream.cpp qvz/src/cluster.cpp qvz/src/well.cpp \
    qvz/src/qv_compressor.cpp qvz/src/pmf.cpp qvz/src/quantizer.cpp qvz/include/qvz.h \
    qvz/src/codebook.cpp qvz/src/distortion.cpp qvz/src/lines.cpp; do
  if [ -f "$f" ] && ! grep -q '#include <cstdint>' "$f"; then
    sed -i '0,/^#include/s//#include <cstdint>\n#include/' "$f"
  fi
done
echo "patched: $(grep -rl "#include <cstdint>" . 2>/dev/null | wc -l) files carry <cstdint>"
