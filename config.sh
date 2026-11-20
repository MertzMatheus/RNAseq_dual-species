#!/usr/bin/env bash
# Pipeline configuration. Paths can also be overridden from the environment, e.g.
#   PROJECT_DIR=/data/my_project bash run_pipeline.sh

# ---- paths -------------------------------------------------------------------
export PROJECT_DIR="${PROJECT_DIR:-/path/to/project}"    # must contain raw_data/ and genomes/
export PROJECT_NAME="${PROJECT_NAME:-dual_rnaseq}"       # used to name the results package
export RAW_DIR="${RAW_DIR:-${PROJECT_DIR}/raw_data}"     # one sub-folder per sample (see samplesheet.tsv)
export GENOME_DIR="${GENOME_DIR:-${PROJECT_DIR}/genomes}"
export OUT_DIR="${OUT_DIR:-${PROJECT_DIR}/analysis}"
export PIPE_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
export SAMPLESHEET="${SAMPLESHEET:-${PIPE_DIR}/samplesheet.tsv}"

# ---- reference genomes ---------------------------------------------------------
export CALB_FNA="${GENOME_DIR}/C_albicans_SC5314.fna"
export SAUR_FNA="${GENOME_DIR}/S_aureus_USA300_FPR3757.fna"
# GFF3 files must use the same sequence names as the FASTA files;
# 00_get_annotations.sh downloads them from NCBI RefSeq when missing.
export CALB_GFF="${GENOME_DIR}/C_albicans_SC5314.gff"
export SAUR_GFF="${GENOME_DIR}/S_aureus_USA300_FPR3757.gff"
export CALB_ACC="GCF_000182965.3"   # C. albicans SC5314, Assembly 22, haplotype A (ASM18296v3)
export SAUR_ACC="GCF_000013465.1"   # S. aureus USA300_FPR3757, chromosome + pUSA01-03

# ---- compute -------------------------------------------------------------------
export THREADS=${THREADS:-16}              # threads per job
export PAR_JOBS=${PAR_JOBS:-2}             # samples in parallel (THREADS x PAR_JOBS <= cores)
export STAR_RAM=${STAR_RAM:-16000000000}   # bytes for STAR BAM sorting

# ---- conda environments (envs/) ------------------------------------------------
export ENV_TOOLS="rnaseq_tools"
export ENV_R="rnaseq_R"

# ---- trimming (fastp) ----------------------------------------------------------
export MIN_LEN=36
export QUAL_CUT=20

# ---- rRNA removal (SortMeRNA) --------------------------------------------------
export RRNA_FILTER="yes"   # "no" skips step 02b
export SORTMERNA_DB_DIR="${SORTMERNA_DB_DIR:-${GENOME_DIR}/sortmerna_db}"
export SORTMERNA_DB="${SORTMERNA_DB:-${SORTMERNA_DB_DIR}/smr_v4.3_default_db.fasta}"

# ---- alignment (STAR) ----------------------------------------------------------
# C. albicans introns are short; capping intron and mate-gap length prevents spurious
# splicing in the bacterial genome.
export MAX_INTRON=5000
export MAX_MULTIMAP=10
export READ_LEN_MINUS1=149   # read length - 1, for the splice-junction index

# ---- counting (featureCounts) --------------------------------------------------
export STRAND="auto"   # auto (inferred from the data), 0, 1 or 2
export MIN_MAPQ=10

# ---- differential expression ---------------------------------------------------
export PADJ_CUT=0.05
export LFC_CUT=1
export MIN_COUNT=10        # keep genes with >= MIN_COUNT reads in >= MIN_SAMPLES samples
export MIN_SAMPLES=3
export DROP_RRNA_TRNA="yes"

# ---- helpers -------------------------------------------------------------------
activate_env () {
  if command -v conda >/dev/null 2>&1; then
    eval "$(conda shell.bash hook)"; conda activate "$1"
  elif command -v micromamba >/dev/null 2>&1; then
    eval "$(micromamba shell hook -s bash)"; micromamba activate "$1"
  else
    echo "[WARN] conda/micromamba not found; assuming tools are on PATH" >&2
  fi
}
log () { echo "[$(date '+%F %T')] $*"; }
export -f activate_env log
