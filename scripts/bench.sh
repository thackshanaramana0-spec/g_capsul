#!/bin/bash
# Per-stage benchmark harness for the reimpl prototypes.
#
# Wall-clock totals are not usable on this machine: another Claude session runs
# concurrently, and the same binary has measured 1.32 s and 1.91 s for the same
# stage. Two things make a number trustworthy here:
#
#   MIN over repetitions -- contention only ever ADDS time, so the minimum is the
#   least contaminated estimate, where the mean is dragged by whatever else ran.
#
#   PER STAGE -- a change to one stage is judged on that stage's timer, so noise
#   in the other six cannot drown it.
#
# Also reports the spread (min..max) so a result is only believed when the spread
# is small relative to the effect being claimed.
#
#   ./bench.sh <reps> <binary> [args...]
set -u
REPS=${1:?usage: bench.sh <reps> <binary> [args...]}; shift
BIN=$1; shift
TMP=$(mktemp -d); trap 'rm -rf "$TMP"' EXIT

echo "# $BIN $*"
echo "# load before: $(uptime | sed 's/.*load average: //')"
for r in $(seq 1 "$REPS"); do
    "$BIN" "$@" > "$TMP/out.$r" 2> "$TMP/err.$r"
done
echo "# load after:  $(uptime | sed 's/.*load average: //')"

python3 - "$TMP" "$REPS" <<'PY'
import re,sys,glob,os
tmp,reps=sys.argv[1],int(sys.argv[2])
runs=[]
for r in range(1,reps+1):
    st={}
    for line in open(f"{tmp}/err.{r}"):
        m=re.match(r"\s{2}(\S.*?)\s{2,}([\d.]+) s",line)
        if m: st[m.group(1).strip()]=float(m.group(2))
    runs.append(st)
keys=[k for k in runs[0]] if runs else []
print(f"{'stage':<26} {'min':>7} {'max':>7} {'spread':>7}")
tot_min=0.0
for k in keys:
    v=[r[k] for r in runs if k in r]
    if not v: continue
    tot_min+=min(v)
    sp=max(v)-min(v)
    flag="  <-- noisy" if sp>0.15*max(min(v),0.01) else ""
    print(f"{k:<26} {min(v):7.2f} {max(v):7.2f} {sp:7.2f}{flag}")
print(f"{'TOTAL (sum of mins)':<26} {tot_min:7.2f}")
# outputs must be identical across reps or the timing means nothing
outs={open(f"{tmp}/out.{r}").read() for r in range(1,reps+1)}
print(f"\noutput identical across {reps} reps: {'YES' if len(outs)==1 else 'NO -- results not comparable'}")
for line in sorted(outs)[0].strip().split("\n"): print(f"  {line}")
PY
