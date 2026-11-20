#!/usr/bin/env bash
# Step 9: GO BP, GO MF and KEGG over-representation analysis of every gene set from step 7.
# Requires the term tables from 08_functional_annotation.sh.
set -euo pipefail
source "$(dirname "$0")/config.sh"
activate_env "$ENV_R"
ORA="$OUT_DIR/07_functional_enrichment"; mkdir -p "$ORA" "$OUT_DIR/logs"
Rscript "$PIPE_DIR/scripts/enrichment.R" --ora_input "$OUT_DIR/06_differential_expression/ora_input.rds" \
  --term_dir "$OUT_DIR/reference/functional" --outdir "$ORA" --padj "$PADJ_CUT" \
  2>&1 | tee "$OUT_DIR/logs/enrichment.log"
log "Results: $ORA"
