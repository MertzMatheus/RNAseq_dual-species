#!/usr/bin/env bash
# Step 6: QC table for every step (QC_summary.xlsx) + MultiQC HTML report.
set -euo pipefail
source "$(dirname "$0")/config.sh"
activate_env "$ENV_TOOLS"
python3 "$PIPE_DIR/scripts/qc_table.py" --outdir "$OUT_DIR" --samplesheet "$SAMPLESHEET" --min-count "$MIN_COUNT"
multiqc -f -q -o "$OUT_DIR/qc/multiqc" -n multiqc_report.html \
  "$OUT_DIR/01_fastqc_raw" "$OUT_DIR/02_trimmed" "$OUT_DIR/02b_rRNA_removed" "$OUT_DIR/03_fastqc_trimmed" \
  "$OUT_DIR/04_alignment" "$OUT_DIR/05_counts" || log "[WARN] MultiQC failed (QC table is still complete)"
log "QC report: $OUT_DIR/qc/QC_summary.xlsx and $OUT_DIR/qc/multiqc/multiqc_report.html"
