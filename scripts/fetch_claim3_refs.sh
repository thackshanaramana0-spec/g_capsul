#!/bin/bash
# ═══════════════════════════════════════════════════════════════════════════
#  Claim 3 reference genomes — fetch, verify, bwa-index.
#
#  T6b compares `capsule_decode coverage` against bwa+samtools+mosdepth, and
#  the conventional route needs a REFERENCE per dataset. These are the five
#  non-human references for the six Claim-3 datasets (the sixth, HG002, uses
#  chr20.fa, which is already present and indexed).
#
#  ACCESSIONS ARE PINNED AND STRAIN-MATCHED to the locked FASTQ accessions --
#  a reference for the wrong strain silently degrades the alignment baseline
#  and makes the comparison unfair to the competitor, which is worse than not
#  running it. Verified against NCBI Datasets, 2026-09-05:
#
#    SRR2584863  E. coli B REL606        -> GCF_000017985.1  (ASM1798v1)
#    ERR5181310  SARS-CoV-2              -> NC_045512.2      (Wuhan-Hu-1)
#    SRR29296997 H. salinarum NRC-1      -> GCF_000006805.1  (ASM680v1)
#    SRR37283774 P. falciparum 3D7       -> GCF_000002765.6  (ASM276v2)
#    DRR976266   S. cerevisiae S288C     -> GCF_000146045.2  (R64)
#
#  Idempotent: a reference that is already present and indexed is skipped, so
#  re-running costs nothing and a partial download is retried.
#
#  usage: bash scripts/fetch_claim3_refs.sh [REFS_DIR]
# ═══════════════════════════════════════════════════════════════════════════
set -u
REFS="${1:-$HOME/refs}"
mkdir -p "$REFS"
cd "$REFS"
log(){ echo "[refs] $(date '+%H:%M:%S') $*"; }
FAIL=0

# name|accession|NCBI FTP path
SETS=(
"ecoli|GCF_000017985.1|https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/017/985/GCF_000017985.1_ASM1798v1/GCF_000017985.1_ASM1798v1_genomic.fna.gz"
"sarscov2|NC_045512.2|https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/009/858/895/GCF_009858895.2_ASM985889v3/GCF_009858895.2_ASM985889v3_genomic.fna.gz"
"halobacterium|GCF_000006805.1|https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/006/805/GCF_000006805.1_ASM680v1/GCF_000006805.1_ASM680v1_genomic.fna.gz"
"pfalciparum|GCF_000002765.6|https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/002/765/GCF_000002765.6_GCA_000002765/GCF_000002765.6_GCA_000002765_genomic.fna.gz"
"scerevisiae|GCF_000146045.2|https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/146/045/GCF_000146045.2_R64/GCF_000146045.2_R64_genomic.fna.gz"
)

for entry in "${SETS[@]}"; do
    IFS='|' read -r name acc url <<< "$entry"
    fa="$REFS/c3_${name}.fa"
    echo ""
    log "── $name ($acc)"
    if [ -s "$fa" ] && [ -s "$fa.bwt" ]; then
        log "   already present and indexed ($(du -h "$fa" | cut -f1)) — skipping"
        continue
    fi
    if [ ! -s "$fa" ]; then
        log "   downloading ..."
        if ! curl -fsSL --retry 3 --retry-delay 5 "$url" -o "$fa.gz.tmp"; then
            log "   ✗ DOWNLOAD FAILED — $url"; FAIL=$((FAIL+1)); rm -f "$fa.gz.tmp"; continue
        fi
        # verify it is really gzip before trusting it (a 404 page is not)
        if ! gzip -t "$fa.gz.tmp" 2>/dev/null; then
            log "   ✗ NOT A VALID GZIP (likely a 404 page) — discarding"; FAIL=$((FAIL+1)); rm -f "$fa.gz.tmp"; continue
        fi
        gunzip -c "$fa.gz.tmp" > "$fa" && rm -f "$fa.gz.tmp"
    fi
    # verify FASTA content
    if ! head -1 "$fa" | grep -q '^>'; then
        log "   ✗ NOT A FASTA (no '>' header) — discarding"; FAIL=$((FAIL+1)); rm -f "$fa"; continue
    fi
    NSEQ=$(grep -c '^>' "$fa"); NBP=$(grep -v '^>' "$fa" | tr -d '\n' | wc -c)
    log "   ✓ $(du -h "$fa" | cut -f1)  ${NSEQ} sequence(s)  ${NBP} bp"
    log "   bwa index ..."
    if bwa index "$fa" >/dev/null 2>&1; then log "   ✓ indexed"; else log "   ✗ bwa index FAILED"; FAIL=$((FAIL+1)); fi
    samtools faidx "$fa" 2>/dev/null || true
done

# HG002's reference is chr20, already present -- index it if it is not
echo ""
log "── HG002 (chr20.fa, already present)"
if [ -s "$REFS/chr20.fa" ]; then
    if [ -s "$REFS/chr20.fa.bwt" ]; then log "   already indexed — skipping"
    else log "   bwa index (this one is 64 MB, takes ~1 min) ..."
         bwa index "$REFS/chr20.fa" >/dev/null 2>&1 && log "   ✓ indexed" || { log "   ✗ FAILED"; FAIL=$((FAIL+1)); }; fi
else log "   ✗ chr20.fa MISSING"; FAIL=$((FAIL+1)); fi

echo ""
log "════════════════════════════════════════════════════════════"
if [ "$FAIL" -eq 0 ]; then
    log "ALL CLAIM 3 REFERENCES READY"
    ls -la "$REFS"/c3_*.fa 2>/dev/null | awk '{printf "       %-34s %8.1f MB\n",$9,$5/1048576}'
else
    log "$FAIL reference(s) FAILED — T6b coverage rows for those will be skipped"
fi
log "════════════════════════════════════════════════════════════"
exit $FAIL
