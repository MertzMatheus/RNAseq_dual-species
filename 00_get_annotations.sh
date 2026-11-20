#!/usr/bin/env bash
# Step 0: make sure each genome has a GFF3 annotation with matching sequence names.
# Downloads the RefSeq genome + GFF3 from NCBI when the GFF is missing.
set -euo pipefail
source "$(dirname "$0")/config.sh"
activate_env "$ENV_TOOLS"
mkdir -p "$GENOME_DIR/ncbi_download"

fetch () {   # $1 accession  $2 target gff  $3 label
  local acc=$1 gff=$2 lab=$3 d="$GENOME_DIR/ncbi_download/$1"
  if [[ -s "$gff" ]]; then log "$lab: GFF found ($gff)"; return; fi
  log "$lab: GFF not found - downloading $acc from NCBI"
  mkdir -p "$d"
  if command -v datasets >/dev/null 2>&1; then
    datasets download genome accession "$acc" --include genome,gff3 --filename "$d/$acc.zip"
    unzip -o -q "$d/$acc.zip" -d "$d"
    cp "$d"/ncbi_dataset/data/"$acc"/genomic.gff "$gff"
    cp "$d"/ncbi_dataset/data/"$acc"/*_genomic.fna "$d/${lab}_NCBI.fna"
  else   # FTP fallback
    case $acc in
      GCF_000182965.3) base="https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/182/965/GCF_000182965.3_ASM18296v3/GCF_000182965.3_ASM18296v3" ;;
      GCF_000013465.1) base="https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/013/465/GCF_000013465.1_ASM1346v1/GCF_000013465.1_ASM1346v1" ;;
      *) echo "Unknown accession $acc - download the GFF3 manually"; exit 1 ;;
    esac
    wget -q -O - "${base}_genomic.gff.gz" | gunzip -c > "$gff"
    wget -q -O - "${base}_genomic.fna.gz" | gunzip -c > "$d/${lab}_NCBI.fna"
  fi
}

fetch "$CALB_ACC" "$CALB_GFF" "C_albicans"
fetch "$SAUR_ACC" "$SAUR_GFF" "S_aureus"

# ---- check that FASTA and GFF use the same contig names ----
check () {
  local fna=$1 gff=$2 lab=$3
  comm -12 <(grep '^>' "$fna" | sed 's/^>//; s/ .*//' | sort -u) \
           <(grep -v '^#' "$gff" | cut -f1 | sort -u) > /tmp/shared_$$.txt
  local n_fa n_sh
  n_fa=$(grep -c '^>' "$fna"); n_sh=$(wc -l < /tmp/shared_$$.txt)
  echo -e "$lab\tFASTA contigs: $n_fa\tshared with GFF: $n_sh"
  if [[ $n_sh -eq 0 ]]; then
    echo "[ERROR] $lab: FASTA and GFF contig names do not match."
    echo "        FASTA: $(grep -m3 '^>' "$fna" | tr '\n' ' ')"
    echo "        GFF  : $(grep -v '^#' "$gff" | cut -f1 | sort -u | head -3 | tr '\n' ' ')"
    echo "        -> use the NCBI FASTA downloaded next to the GFF (ncbi_download/), or a GFF from the"
    echo "           same source as your FASTA (e.g. CGD for CGD-style 'Ca22chr1A_C_albicans_SC5314' names)."
    exit 1
  fi
  rm -f /tmp/shared_$$.txt
}
check "$CALB_FNA" "$CALB_GFF" "C_albicans"
# C. albicans must be a HAPLOID representation (haplotype A, 8 chromosomes [+ mtDNA]).
# A diploid A+B FASTA (CGD) would make nearly every read a multi-mapper and they would be discarded.
if grep '^>' "$CALB_FNA" | grep -qiE 'chr[0-9R]B|haplotype B|_B_C_albicans' || [[ $(grep -c '^>' "$CALB_FNA") -gt 10 ]]; then
  echo "[ERROR] $CALB_FNA looks like a diploid (haplotype A + B) assembly. Use haplotype A only"
  echo "        (NCBI GCF_000182965.3 / ASM18296v3, or CGD 'A22 haplotype A' chromosomes)."; exit 1
fi
check "$SAUR_FNA" "$SAUR_GFF" "S_aureus"
log "Annotations OK"
