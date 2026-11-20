#!/usr/bin/env bash
# Step 2b: in-silico rRNA removal with SortMeRNA (SILVA 16S/18S/23S/28S + Rfam 5S/5.8S, all domains).
# A pair is removed if either mate matches rRNA (--paired_in). The removed fraction is reported in the QC
# table; it includes reads from multi-copy rRNA operons that genome-based counts miss as multi-mappers.
set -euo pipefail
source "$(dirname "$0")/config.sh"
activate_env "$ENV_TOOLS"
OUTR="$OUT_DIR/02b_rRNA_removed"; mkdir -p "$OUTR" "$OUT_DIR/logs" "$SORTMERNA_DB_DIR"

if [[ "$RRNA_FILTER" != "yes" ]]; then log "RRNA_FILTER=$RRNA_FILTER - skipping rRNA removal"; exit 0; fi

# ---- database (downloaded once, ~240 MB) ----
if [[ ! -s "$SORTMERNA_DB" ]]; then
  log "Downloading SortMeRNA rRNA databases (v4.3.4 release)"
  wget -q -O - https://github.com/biocore/sortmerna/releases/download/v4.3.4/database.tar.gz | tar xz -C "$SORTMERNA_DB_DIR"
fi
IDX="$SORTMERNA_DB_DIR/idx_$(basename "$SORTMERNA_DB" .fasta)"
if [[ ! -e "$IDX/.done" ]]; then
  log "Indexing $(basename "$SORTMERNA_DB") (once; can take ~30-60 min)"
  rm -rf "$IDX" "$SORTMERNA_DB_DIR/idx_tmp"; mkdir -p "$IDX"
  sortmerna --ref "$SORTMERNA_DB" --index 1 --idx-dir "$IDX" --workdir "$SORTMERNA_DB_DIR/idx_tmp" \
            --threads "$THREADS" > "$OUT_DIR/logs/sortmerna_index.log" 2>&1
  rm -rf "$SORTMERNA_DB_DIR/idx_tmp"; touch "$IDX/.done"
fi

one_sample () {
  local s=$1 wd="$OUTR/tmp_$1"
  if [[ -s "$OUTR/${s}_R1.nonrRNA.fq.gz" && -s "$OUTR/${s}.sortmerna.log" ]]; then log "$s rRNA-filtered - skipping"; return; fi
  log "$s: SortMeRNA"
  rm -rf "$wd"; mkdir -p "$wd"
  sortmerna --ref "$SORTMERNA_DB" --idx-dir "$IDX" --index 0 --workdir "$wd" --threads "$THREADS" \
            --reads "$OUT_DIR/02_trimmed/${s}_R1.trim.fq.gz" --reads "$OUT_DIR/02_trimmed/${s}_R2.trim.fq.gz" \
            --fastx --paired_in --out2 --num_alignments 1 \
            --aligned "$wd/${s}_rRNA" --other "$wd/${s}_nonrRNA" > "$OUT_DIR/logs/${s}.sortmerna.stdout" 2>&1 \
    || { log "[ERROR] $s: SortMeRNA failed (logs/${s}.sortmerna.stdout; often out of memory: lower PAR_JOBS)"; return 1; }
  # output file names differ between SortMeRNA versions
  local f r
  f=$(ls "$wd/${s}_nonrRNA"*fwd* 2>/dev/null | head -1); r=$(ls "$wd/${s}_nonrRNA"*rev* 2>/dev/null | head -1)
  [[ -n "$f" && -n "$r" ]] || { log "[ERROR] $s: SortMeRNA produced no output"; return 1; }
  for x in "$f:R1" "$r:R2"; do
    src=${x%:*}; rd=${x##*:}
    if [[ "$src" == *.gz ]]; then mv "$src" "$OUTR/${s}_${rd}.nonrRNA.fq.gz"
    else pigz -p "$THREADS" -c "$src" > "$OUTR/${s}_${rd}.nonrRNA.fq.gz"; fi
  done
  cp "$wd/${s}_rRNA.log" "$OUTR/${s}.sortmerna.log"
  rm -rf "$wd"   # rRNA reads and the per-sample key-value store are not kept
}
export -f one_sample
export OUTR IDX

tail -n +2 "$SAMPLESHEET" | cut -f1 | xargs -P "$PAR_JOBS" -I{} bash -c 'one_sample {}' \
  || { log "[ERROR] rRNA removal failed for at least one sample"; exit 1; }
log "rRNA removal done"
