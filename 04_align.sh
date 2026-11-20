#!/usr/bin/env bash
# Step 4: STAR alignment to the combined genome + per-species read accounting.
set -euo pipefail
source "$(dirname "$0")/config.sh"
activate_env "$ENV_TOOLS"
REF="$OUT_DIR/reference"; ALN="$OUT_DIR/04_alignment"; mkdir -p "$ALN" "$OUT_DIR/logs"

# input: rRNA-depleted reads (step 02b) or trimmed reads
if [[ "$RRNA_FILTER" == "yes" ]]; then IN_PREFIX="$OUT_DIR/02b_rRNA_removed/"; IN_SUFFIX=".nonrRNA.fq.gz"
else IN_PREFIX="$OUT_DIR/02_trimmed/"; IN_SUFFIX=".trim.fq.gz"; fi
export IN_PREFIX IN_SUFFIX

one_sample () {
  local s=$1 bam="$ALN/$1.Aligned.sortedByCoord.out.bam"
  if [[ -s "$bam.bai" && -s "$ALN/$s.species_counts.tsv" ]]; then log "$s aligned - skipping"; return; fi
  log "$s: STAR"
  STAR --runThreadN "$THREADS" --genomeDir "$REF/STAR_index" \
       --readFilesIn "${IN_PREFIX}${s}_R1${IN_SUFFIX}" "${IN_PREFIX}${s}_R2${IN_SUFFIX}" \
       --readFilesCommand zcat \
       --outFileNamePrefix "$ALN/$s." \
       --outSAMtype BAM SortedByCoordinate --limitBAMsortRAM "$STAR_RAM" \
       --outSAMattributes NH HI AS nM NM MD \
       --outFilterType BySJout \
       --outFilterMultimapNmax "$MAX_MULTIMAP" \
       --alignIntronMin 20 --alignIntronMax "$MAX_INTRON" --alignMatesGapMax "$MAX_INTRON" \
       --alignSJoverhangMin 8 --alignSJDBoverhangMin 1 \
       --outSAMunmapped Within > "$OUT_DIR/logs/$s.star.log" 2>&1
  samtools index -@ "$THREADS" "$bam"
  samtools flagstat -@ "$THREADS" "$bam" > "$ALN/$s.flagstat.txt"
  samtools idxstats "$bam" > "$ALN/$s.idxstats.txt"

  # read PAIRS (first mate, primary) per species: all mapped / uniquely mapped (MAPQ 255)
  {
    echo -e "sample\tspecies\tpairs_mapped\tpairs_unique\tmismatch_rate"
    for sp in Calb Saur; do
      mapfile -t C < "$REF/$sp.contigs"
      all=$(samtools view -c -@ 4 -f 0x40 -F 0x904 "$bam" "${C[@]}")
      uni=$(samtools view -c -@ 4 -f 0x40 -F 0x904 -q 255 "$bam" "${C[@]}")
      # error (mismatch) rate of uniquely mapped reads on this species' contigs; reflects isolate divergence vs reference
      err=$(samtools view -u -q 255 -F 0x904 "$bam" "${C[@]}" | samtools stats - 2>/dev/null | \
            awk -F'\t' '$1=="SN" && $2=="error rate:"{print $3}')
      echo -e "$s\t$sp\t$all\t$uni\t${err:-NA}"
    done
  } > "$ALN/$s.species_counts.tsv"
}
export -f one_sample
export REF ALN

tail -n +2 "$SAMPLESHEET" | cut -f1 | xargs -P "$PAR_JOBS" -I{} bash -c 'one_sample {}' \
  || { log "[ERROR] alignment failed for at least one sample"; exit 1; }
log "Alignment done"
