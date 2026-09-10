#!/usr/bin/env bash
# Confirm a fresh server has everything, BEFORE spending hours on data or runs.
#
# This is deliberately separate from scripts/benchmark_0_preflight.sh: that one
# also needs the datasets, the references, the truth sets and a working build.
# This one answers the earlier question -- "is the machine and the toolchain
# right?" -- so a setup problem is found in seconds instead of after a download.
set -u
OK=0; BAD=0
p(){ printf "  [ OK ] %-30s %s\n" "$1" "${2:-}"; OK=$((OK+1)); }
f(){ printf "  [FAIL] %-30s %s\n" "$1" "${2:-}"; BAD=$((BAD+1)); }
h(){ printf "\n== %s ==\n" "$1"; }

h "machine"
C=$(nproc); M=$(free -g | awk '/^Mem:/{print $2}'); D=$(df -BG / | awk 'NR==2{gsub("G","",$4);print $4}')
[ "$C" -ge 8 ]   && p "cores >= 8"        "$C"      || f "cores < 8"        "$C (slower, not fatal)"
[ "$M" -ge 32 ]  && p "RAM >= 32 GB"      "${M} GB" || f "RAM < 32 GB"      "${M} GB -- peak measured 19.9 GB"
[ "$D" -ge 60 ]  && p "disk >= 60 GB free" "${D} GB" || f "disk < 60 GB free" "${D} GB -- peak transient ~36 GB"

h "the /dev/stdout trap (genozip writes NO archive without it)"
[ -e /dev/stdout ] && p "/dev/stdout exists" || f "/dev/stdout MISSING" "fix: ln -sfn /proc/self/fd/1 /dev/stdout"

h "apt tools"
for t in cmake gcc g++ samtools bcftools tabix bwa mosdepth seqtk kmc; do
  v=$(command -v $t 2>/dev/null) && p "$t" "$v" || f "$t" "apt-get install -y $t"
done

h "built-from-source tools"
for t in spring genozip rtg aws prefetch fasterq-dump; do
  v=$(command -v $t 2>/dev/null) && p "$t" "$v" || f "$t" "see SETUP_FULL.md"
done

h "competitor callers / assembler (paths, not PATH)"
for pth in /root/SPAdes-4.0.0-Linux/bin/spades.py /root/DiscoSnp/run_discoSnp++.sh \
           /root/Kmer2SNP/kmer2snp.py /root/miniconda3/envs/kmer2snp_r/bin/python; do
  [ -e "$pth" ] && p "$(basename "$pth")" "$pth" || f "$(basename "$pth")" "missing: $pth"
done
[ -x /root/miniconda3/envs/kmer2snp_r/bin/python ] && {
  /root/miniconda3/envs/kmer2snp_r/bin/python -c 'import networkx' 2>/dev/null \
    && p "networkx in conda env" || f "networkx in conda env" "kmer2snp.py cannot import"; }
PATH="$HOME/DiscoSnp:$PATH" command -v run_discoSnp++.sh >/dev/null \
  && p "DiscoSNP++ reachable on PATH" || f "DiscoSNP++ not on PATH" "export PATH=\$HOME/DiscoSnp:\$PATH"

h "genozip actually produces an archive (>1MB probe -- smaller ones pass falsely)"
T=$(mktemp -d)
awk 'BEGIN{srand(3);for(i=0;i<12000;i++){s="";q="";for(j=0;j<151;j++){s=s substr("ACGT",int(rand()*4)+1,1);q=q "I"}
     printf "@p%d\n%s\n+\n%s\n",i,s,q}}' > "$T/p.fq"
if genozip --force -o "$T/p.genozip" "$T/p.fq" >/dev/null 2>"$T/err" && [ -s "$T/p.genozip" ]; then
  p "genozip round trip" "$(stat -c%s "$T/p.fq") -> $(stat -c%s "$T/p.genozip") B"
else
  f "genozip produces NO archive" "$(grep -aoiE 'LICENSE ERROR.*|Neither curl nor wget.*' "$T/err" | head -1)"
fi
rm -rf "$T"

printf "\n  passed %d, failed %d\n" "$OK" "$BAD"
[ "$BAD" -eq 0 ] && echo "  ENVIRONMENT OK -- now fetch data, then run scripts/benchmark_0_preflight.sh" \
                 || echo "  FIX THE ABOVE -- see server/SETUP_FULL.md"
exit $([ "$BAD" -eq 0 ] && echo 0 || echo 1)
