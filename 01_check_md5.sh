#!/usr/bin/env bash
# Step 1: verify file integrity (MD5.txt delivered with each sample) and list the FASTQ files.
set -euo pipefail
source "$(dirname "$0")/config.sh"
mkdir -p "$OUT_DIR/qc"
out="$OUT_DIR/qc/md5_check.tsv"
echo -e "sample\tfile\texpected_md5\tobserved_md5\tstatus" > "$out"

tail -n +2 "$SAMPLESHEET" | cut -f1,2 | while read -r sample folder; do
  d="$RAW_DIR/$folder"
  [[ -d "$d" ]] || { echo -e "$sample\t-\t-\t-\tFOLDER_MISSING" >> "$out"; continue; }
  md5file=$(ls "$d"/MD5.txt "$d"/md5.txt "$d"/*MD5*.txt 2>/dev/null | head -1 || true)
  for fq in "$d"/*_1.fq.gz "$d"/*_2.fq.gz; do
    [[ -e "$fq" ]] || continue
    f=$(basename "$fq")
    obs=$(md5sum "$fq" | cut -d' ' -f1)
    exp="NA"; st="NO_MD5_ENTRY"
    if [[ -n "$md5file" ]]; then
      exp=$(awk -v f="$f" '{n=$NF; sub(".*/","",n); if(n==f){print $1; exit}}' "$md5file")
      [[ -z "$exp" ]] && exp="NA"
      if [[ "$exp" != "NA" ]]; then [[ "$exp" == "$obs" ]] && st="OK" || st="MISMATCH"; fi
    fi
    echo -e "$sample\t$f\t$exp\t$obs\t$st" >> "$out"
  done
done

log "MD5 summary:"; cut -f5 "$out" | tail -n +2 | sort | uniq -c
if grep -qP "\tMISMATCH$|\tFOLDER_MISSING$" "$out"; then
  echo "[ERROR] corrupted or missing files - see $out (re-copy those files before continuing)"; exit 1
fi
