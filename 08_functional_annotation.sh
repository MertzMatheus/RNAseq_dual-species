#!/usr/bin/env bash
# Step 8: download GO / KEGG annotations and build gene -> term tables for both species.
#   C. albicans GO : CGD GAF (GO Consortium CANAL-mod.gaf.gz; fallback: CGD website) + UniProt (taxon 237561)
#   S. aureus GO   : UniProt, strain USA300 (taxon 367830)
#   KEGG pathways  : cal (C. albicans SC5314), saa (S. aureus USA300_FPR3757; uses SAUSA300_xxxx old locus tags)
# Raw files are dated and kept, so the analysis can be reproduced with the same annotation release.
set -euo pipefail
source "$(dirname "$0")/config.sh"
activate_env "$ENV_TOOLS"
FUN="$OUT_DIR/reference/functional"; RAW="$FUN/raw"; mkdir -p "$RAW"

get () {   # $1 url  $2 dest  -> returns 1 if it fails
  [[ -s "$2" ]] && return 0
  curl -fsSL --retry 3 --max-time 1800 -o "$2.tmp" "$1" && mv "$2.tmp" "$2" && return 0
  rm -f "$2.tmp"; return 1
}

log "GO ontology"
get https://current.geneontology.org/ontology/go-basic.obo "$RAW/go-basic.obo" || \
get http://purl.obolibrary.org/obo/go/go-basic.obo "$RAW/go-basic.obo"

log "C. albicans GO annotations (CGD)"
get https://current.geneontology.org/annotations/gaf/CANAL-mod.gaf.gz "$RAW/Calb_GO.gaf.gz" || \
get http://www.candidagenome.org/download/go/gene_association.cgd.gz "$RAW/Calb_GO.gaf.gz" || \
log "[WARN] CGD GAF not available - C. albicans GO will come from UniProt only"

UNI="https://rest.uniprot.org/uniprotkb/stream?format=tsv&compressed=true&fields=accession,gene_oln,gene_orf,gene_names,go_p,go_f,go_c&query=organism_id:"
log "UniProt GO"
get "${UNI}237561" "$RAW/Calb_uniprot.tsv.gz" || log "[WARN] UniProt C. albicans download failed"
get "${UNI}367830" "$RAW/Saur_uniprot.tsv.gz" || log "[WARN] UniProt S. aureus download failed"

log "KEGG"
for pair in Calb:cal Saur:saa; do
  sp=${pair%:*}; org=${pair#*:}
  get "https://rest.kegg.jp/link/pathway/$org" "$RAW/${sp}_kegg_link_pathway.txt" || log "[WARN] KEGG $org link failed"
  get "https://rest.kegg.jp/list/pathway/$org" "$RAW/${sp}_kegg_list_pathway.txt" || true
  get "https://rest.kegg.jp/conv/uniprot/$org" "$RAW/${sp}_kegg_conv_uniprot.txt" || true
done
echo "Downloaded: $(date '+%F')  (files in $RAW)" > "$FUN/versions.txt"
ls -l "$RAW" >> "$FUN/versions.txt"

log "Building term tables"
python3 "$PIPE_DIR/scripts/build_term2gene.py" --raw "$RAW" --annot "$OUT_DIR/reference/gene_annotation.tsv" \
  --calb-gff "$CALB_GFF" --saur-gff "$SAUR_GFF" --outdir "$FUN"
log "Coverage: $FUN/annotation_coverage.tsv"
