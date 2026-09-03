#!/bin/bash
# Single source of truth for dataset/tool locations across every CAPSULE
# script. If the datasets move to a new path later, this is the ONE file
# to edit -- nothing else in scripts/ should hardcode a data path.
#
# Every value here is an env-var default: an existing environment variable
# always wins (":=` assigns only if unset), so a one-off override doesn't
# require editing this file either:
#   CAPSULE_DATA_DIR=/mnt/other/fastq bash scripts/run_capsule.sh 1

: "${CAPSULE_DATA_DIR:=/data/fastq}"          # Claim 1 SRA-derived FASTQ, Claim 2 GIAB pooled FASTQ
: "${CAPSULE_REFS_DIR:=$HOME/refs}"           # chr20.fa (Claim 2, Claim 3 coverage baseline)
: "${CAPSULE_TRUTH_DIR:=$HOME/giab_truth}"    # GIAB truth VCFs + confident BEDs (Claim 2)
: "${CAPSULE_OUT_DIR:=./results}"             # where every claim writes its results
: "${CAPSULE_BIN_DIR:=/tmp/capsule_bin}"      # built binaries (best106, capsule_decode)

export CAPSULE_DATA_DIR CAPSULE_REFS_DIR CAPSULE_TRUTH_DIR CAPSULE_OUT_DIR CAPSULE_BIN_DIR
