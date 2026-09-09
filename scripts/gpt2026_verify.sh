#!/bin/bash
# Always decode the actual archive into a new directory before comparing.
set -euo pipefail
GPT2026_ARCHIVE="$(realpath "${1:?archive}")"
GPT2026_INPUT="$(realpath "${2:?original FASTQ}")"
GPT2026_OUT="${3:?new verification directory}"
GPT2026_DEC="${GPT2026_DEC:-/tmp/gpt2026_baseline_decode}"
GPT2026_VERIFY="${GPT2026_VERIFY:-/tmp/gpt2026_fastq_verify}"
test -s "$GPT2026_ARCHIVE"
test -s "$GPT2026_INPUT"
test -x "$GPT2026_DEC"
test -x "$GPT2026_VERIFY"
if [[ -e "$GPT2026_OUT" ]]; then
  echo "Refusing to reuse verification directory: $GPT2026_OUT" >&2
  exit 2
fi
mkdir -p "$GPT2026_OUT/dumps"
"$GPT2026_DEC" "$GPT2026_ARCHIVE" "$GPT2026_OUT/dumps" "$GPT2026_OUT/reads.seq" \
  > "$GPT2026_OUT/decode.stdout" 2> "$GPT2026_OUT/decode.log"
"$GPT2026_VERIFY" "$GPT2026_ARCHIVE" "$GPT2026_OUT/reads.seq" "$GPT2026_INPUT" \
  > "$GPT2026_OUT/verify.txt" 2> "$GPT2026_OUT/verify.log"
cat "$GPT2026_OUT/verify.txt"
