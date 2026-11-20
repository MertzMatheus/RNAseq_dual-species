#!/usr/bin/env bash
# Step 10: collect the files needed for the manuscript (methods, QC, figures, tables) into one folder + zip.
#   INCLUDE_TIFF=no      leave out the 600-dpi TIFF copies of the figures
#   RCLONE_REMOTE=...    also copy the folder to an rclone remote (e.g. gdrive:project/results)
set -euo pipefail
source "$(dirname "$0")/config.sh"
NAME="${PROJECT_NAME}_results"; S="$OUT_DIR/$NAME"
DE="$OUT_DIR/06_differential_expression"; ORA="$OUT_DIR/07_functional_enrichment"
EXT=(pdf png); [[ "${INCLUDE_TIFF:-yes}" == "yes" ]] && EXT+=(tiff)
rm -rf "$S" "$S.zip"; mkdir -p "$S"/{01_methods,02_quality_control,03_figures,04_tables}

figs () {   # $1 source dir, $2 destination dir
  [[ -d "$1" ]] || return 0; mkdir -p "$2"
  for e in "${EXT[@]}"; do find "$1" -maxdepth 1 -name "*.$e" -exec cp {} "$2/" \; ; done
}
cp "$PIPE_DIR/docs/METHODS.md" "$S/01_methods/"
for f in "$OUT_DIR/qc/QC_summary.xlsx" "$OUT_DIR/qc/multiqc/multiqc_report.html" "$OUT_DIR/reference/functional/annotation_coverage.tsv"; do
  [[ -e "$f" ]] && cp "$f" "$S/02_quality_control/"
done

declare -A NM=( [Overview]=00_overview [C1]=C1_Saur_isolates_alone [C2a]=C2a_Saur_isolates_with_Calb
                [C2b]=C2b_Saur_response_to_Calb [C3]=C3_Calb_by_partner_isolate [C4]=C4_Calb_response_to_each_isolate
                [C1_vs_C2a]=C1_vs_C2a_concordance [C3xC4]=C3xC4_isolate_specific_fungal_response )
for k in "${!NM[@]}"; do
  figs "$DE/figures/$k"  "$S/03_figures/${NM[$k]}"
  figs "$ORA/figures/$k" "$S/03_figures/${NM[$k]}/enrichment"
done
[[ -d "$S/03_figures/00_overview/enrichment" ]] && mv "$S"/03_figures/00_overview/enrichment/* "$S/03_figures/00_overview/" && rmdir "$S/03_figures/00_overview/enrichment"
cp "$DE"/tables/*.xlsx "$ORA"/tables/*.xlsx "$S/04_tables/" 2>/dev/null || true
find "$S" -type d -empty -delete

cat > "$S/README.txt" <<TXT
Dual RNA-seq: C. albicans SC5314 x S. aureus (JE2, CI19, S5, CI5) biofilms
Generated $(date '+%Y-%m-%d')

01_methods/          methods text (bioinformatics)
02_quality_control/  per-sample QC table (QC_summary.xlsx), MultiQC report, GO/KEGG annotation coverage
03_figures/          PDF (vector), TIFF (600 dpi) and PNG; 7 pt Arial, 89/183 mm width
  00_overview/       PCA, sample correlation, DEG counts, expression heatmaps of all DEGs and
                     enrichment of their clusters, colour key
  <comparison>/      volcano plots, Venn/UpSet diagrams; enrichment/ holds <comparison>_ORA_summary
                     (all enriched terms) and <contrast>_ORA (dot plot + gene-term network)
04_tables/           DEG_summary.xlsx (overview), <comparison>_DE.xlsx (all genes, every contrast),
                     <comparison>_shared_specific_genes.xlsx, <comparison>_ORA.xlsx,
                     Expression_heatmap_clusters.xlsx

Comparisons (log2FC > 0 = higher in the first group; DEG: adjusted P < ${PADJ_CUT}, |log2FC| >= ${LFC_CUT})
  C1         S. aureus genes, isolates alone: all 6 isolate pairs
  C2a        S. aureus genes, the same 6 pairs in co-culture with C. albicans
  C2b        S. aureus genes, each isolate with vs without C. albicans
  C3         C. albicans genes, the 6 pairs of co-cultures (Ca+CI19 vs Ca+JE2, ...)
  C4         C. albicans genes, each co-culture vs C. albicans alone, plus the mean of the 4 co-cultures
  C1_vs_C2a  concordance of isolate differences without and with C. albicans
  C3xC4      C. albicans genes specific to one partner isolate and changed vs C. albicans alone

Isolate-specific genes (C1, C2a, C3): DE in all three contrasts of an isolate, in the same direction.
Shared/specific responses (C2b, C4): DE with all four isolates / with one isolate only.
TXT
( cd "$S" && find . -type f | sort ) > "$S/file_list.txt"
( cd "$OUT_DIR" && zip -qr "$NAME.zip" "$NAME" ) || log "[WARN] zip not available"
log "Results package: $S ($(du -sh "$S" | cut -f1))"
if [[ -n "${RCLONE_REMOTE:-}" ]]; then
  command -v rclone >/dev/null || { echo "[ERROR] rclone not found" >&2; exit 1; }
  rclone copy "$S" "$RCLONE_REMOTE/$NAME" --progress
fi
