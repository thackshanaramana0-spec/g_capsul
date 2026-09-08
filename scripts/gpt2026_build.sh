#!/bin/bash
# Additive builds using the shipped build scripts and their exact toolchain.
# Usage: bash scripts/gpt2026_build.sh MODE OUTPUT
# Modes: baseline, profile, owned, memory, combined, combined-tests, verify.
set -euo pipefail
GPT2026_HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
GPT2026_MODE="${1:?select a build mode}"
GPT2026_OUT="${2:?supply output binary path}"
GPT2026_CXX="$(command -v g++)"
GPT2026_TMP="$(mktemp -d /tmp/gpt2026_build.XXXXXX)"
trap 'rm -rf "$GPT2026_TMP"' EXIT
export GPT2026_CXX GPT2026_MODE GPT2026_HERE
# Only inject diagnostic defines or substitute the additive verifier source.
# Compilation/link flags, vendored C compilation and libraries remain governed
# by build106.sh/build_decode.sh. The temporary compiler shim is never installed.
cat > "$GPT2026_TMP/g++" <<'SH'
#!/bin/bash
set -euo pipefail
GPT2026_ARGS=("$@")
case "$GPT2026_MODE" in
  baseline) ;;
  profile) GPT2026_ARGS=(-DGPT2026_MAP_PROFILE "${GPT2026_ARGS[@]}") ;;
  owned) GPT2026_ARGS=(-DGPT2026_MAP_OWNED "${GPT2026_ARGS[@]}") ;;
  memory) GPT2026_ARGS=(-DGPT2026_BUFFER_OWNERSHIP "${GPT2026_ARGS[@]}") ;;
  combined|combined-tests) GPT2026_ARGS=(-DGPT2026_MAP_OWNED -DGPT2026_BUFFER_OWNERSHIP "${GPT2026_ARGS[@]}") ;;
  verify)
    for GPT2026_I in "${!GPT2026_ARGS[@]}"; do
      if [[ "${GPT2026_ARGS[$GPT2026_I]}" == "$GPT2026_HERE/stages/capsule_decode.cpp" ]]; then
        GPT2026_ARGS[$GPT2026_I]="$GPT2026_HERE/stages/gpt2026_fastq_verify.cpp"
      fi
    done ;;
  *) exit 2 ;;
esac
exec "$GPT2026_CXX" "${GPT2026_ARGS[@]}"
SH
chmod +x "$GPT2026_TMP/g++"
if [[ "$GPT2026_MODE" == combined-tests ]]; then
  GPT2026_MAP_OWNED=1 GPT2026_RELEASE_PARENT=1 GPT2026_MOVE_BODIES=1 \
    PATH="$GPT2026_TMP:$PATH" bash "$GPT2026_HERE/scripts/run_tests.sh" > "$GPT2026_OUT" 2>&1
elif [[ "$GPT2026_MODE" == verify ]]; then
  PATH="$GPT2026_TMP:$PATH" bash "$GPT2026_HERE/scripts/build_decode.sh" "$GPT2026_OUT"
else
  PATH="$GPT2026_TMP:$PATH" bash "$GPT2026_HERE/scripts/build106.sh" "$GPT2026_OUT"
fi
