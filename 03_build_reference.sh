#!/usr/bin/env bash
# Step 3: combined C. albicans + S. aureus reference, species-tagged annotation and STAR index.
# Mapping every sample to the concatenated genome lets reads compete between the two species
# (no cross-species mis-assignment) and gives the per-species read split used for QC.
set -euo pipefail
source "$(dirname "$0")/config.sh"
activate_env "$ENV_TOOLS"
REF="$OUT_DIR/reference"; mkdir -p "$REF/STAR_index"

log "Concatenating genomes"
cat "$CALB_FNA" "$SAUR_FNA" > "$REF/combined_Calb_Saur.fna"
samtools faidx "$REF/combined_Calb_Saur.fna"

log "Building species-tagged annotation"
python3 "$PIPE_DIR/scripts/make_annotation.py" \
  --euk-gff "$CALB_GFF" --euk-fna "$CALB_FNA" --euk-name Calb \
  --prok-gff "$SAUR_GFF" --prok-fna "$SAUR_FNA" --prok-name Saur \
  --outdir "$REF"
awk '$2=="Calb"{print $1}' "$REF/contig_species.tsv" > "$REF/Calb.contigs"
awk '$2=="Saur"{print $1}' "$REF/contig_species.tsv" > "$REF/Saur.contigs"

# genomeSAindexNbases = min(14, log2(GenomeLength)/2 - 1)
GLEN=$(awk '{s+=$2} END{print s}' "$REF/combined_Calb_Saur.fna.fai")
SAIDX=$(python3 -c "import math; print(min(14, int(math.log2($GLEN)/2 - 1)))")
log "Genome length $GLEN bp -> genomeSAindexNbases $SAIDX"

SJ_OPT=()
if [[ -s "$REF/calb_exons.gtf" ]]; then
  SJ_OPT=(--sjdbGTFfile "$REF/calb_exons.gtf" --sjdbOverhang "${READ_LEN_MINUS1:-149}")
fi
STAR --runMode genomeGenerate --runThreadN "$THREADS" \
     --genomeDir "$REF/STAR_index" --genomeFastaFiles "$REF/combined_Calb_Saur.fna" \
     --genomeSAindexNbases "$SAIDX" "${SJ_OPT[@]}" \
     --outFileNamePrefix "$REF/STAR_index/"
log "Reference ready: $REF"
