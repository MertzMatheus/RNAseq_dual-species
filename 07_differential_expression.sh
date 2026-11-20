#!/usr/bin/env bash
# Step 7: DESeq2 contrasts (C1-C4), isolate-specific and shared genes, PCA, heatmaps and figures.
set -euo pipefail
source "$(dirname "$0")/config.sh"
activate_env "$ENV_R"
DE="$OUT_DIR/06_differential_expression"; mkdir -p "$DE" "$OUT_DIR/logs"
Rscript "$PIPE_DIR/scripts/differential_expression.R" \
  --counts "$OUT_DIR/05_counts/counts_matrix.tsv" --annot "$OUT_DIR/reference/gene_annotation.tsv" \
  --samples "$SAMPLESHEET" --outdir "$DE" \
  --padj "$PADJ_CUT" --lfc "$LFC_CUT" --min_count "$MIN_COUNT" --min_samples "$MIN_SAMPLES" --drop_rrna "$DROP_RRNA_TRNA" \
  2>&1 | tee "$OUT_DIR/logs/differential_expression.log"
log "Results: $DE"
