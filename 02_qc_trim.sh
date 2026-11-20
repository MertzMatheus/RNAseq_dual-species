#!/usr/bin/env bash
# Step 2: raw read QC (FastQC), adapter/quality trimming (fastp), QC of trimmed reads (FastQC).
set -euo pipefail
source "$(dirname "$0")/config.sh"
activate_env "$ENV_TOOLS"
mkdir -p "$OUT_DIR"/{01_fastqc_raw,02_trimmed,03_fastqc_trimmed,logs}

one_sample () {
  local sample=$1 folder=$2 d="$RAW_DIR/$2"
  local r1 r2
  if [[ -s "$OUT_DIR/02_trimmed/${sample}_R1.trim.fq.gz" && -s "$OUT_DIR/02_trimmed/${sample}.fastp.json" ]]; then
    log "$sample already trimmed - skipping"; return; fi
  mapfile -t R1 < <(ls "$d"/*_1.fq.gz | sort)
  mapfile -t R2 < <(ls "$d"/*_2.fq.gz | sort)
  [[ ${#R1[@]} -eq ${#R2[@]} && ${#R1[@]} -gt 0 ]] || { echo "[ERROR] $sample: R1/R2 files not paired"; return 1; }
  if [[ ${#R1[@]} -gt 1 ]]; then           # several lanes -> concatenate
    log "$sample: merging ${#R1[@]} lanes"
    cat "${R1[@]}" > "$OUT_DIR/02_trimmed/${sample}_R1.raw.fq.gz"
    cat "${R2[@]}" > "$OUT_DIR/02_trimmed/${sample}_R2.raw.fq.gz"
    r1="$OUT_DIR/02_trimmed/${sample}_R1.raw.fq.gz"; r2="$OUT_DIR/02_trimmed/${sample}_R2.raw.fq.gz"
  else
    r1=${R1[0]}; r2=${R2[0]}
  fi
  log "$sample: FastQC raw"
  fastqc -q -t 2 -o "$OUT_DIR/01_fastqc_raw" "$r1" "$r2"
  log "$sample: fastp"
  fastp -i "$r1" -I "$r2" \
        -o "$OUT_DIR/02_trimmed/${sample}_R1.trim.fq.gz" -O "$OUT_DIR/02_trimmed/${sample}_R2.trim.fq.gz" \
        --detect_adapter_for_pe --trim_poly_g \
        --cut_right --cut_right_window_size 4 --cut_right_mean_quality "$QUAL_CUT" \
        --length_required "$MIN_LEN" --n_base_limit 5 \
        --thread "$THREADS" \
        --json "$OUT_DIR/02_trimmed/${sample}.fastp.json" --html "$OUT_DIR/02_trimmed/${sample}.fastp.html" \
        --report_title "$sample" 2> "$OUT_DIR/logs/${sample}.fastp.log"
  log "$sample: FastQC trimmed"
  fastqc -q -t 2 -o "$OUT_DIR/03_fastqc_trimmed" \
         "$OUT_DIR/02_trimmed/${sample}_R1.trim.fq.gz" "$OUT_DIR/02_trimmed/${sample}_R2.trim.fq.gz"
  rm -f "$OUT_DIR/02_trimmed/${sample}_R1.raw.fq.gz" "$OUT_DIR/02_trimmed/${sample}_R2.raw.fq.gz"
}
export -f one_sample

tail -n +2 "$SAMPLESHEET" | cut -f1,2 | \
  xargs -P "$PAR_JOBS" -L 1 bash -c 'one_sample "$0" "$1"'
log "Trimming done"
