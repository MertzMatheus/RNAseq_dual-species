#!/usr/bin/env bash
# Step 5: gene-level counting with featureCounts on the combined annotation.
# Strandedness is inferred from the data: the library is counted as unstranded (0),
# forward (1) and reverse (2); the orientation with clearly more assigned pairs wins.
set -euo pipefail
source "$(dirname "$0")/config.sh"
activate_env "$ENV_TOOLS"
REF="$OUT_DIR/reference"; ALN="$OUT_DIR/04_alignment"; CNT="$OUT_DIR/05_counts"; mkdir -p "$CNT"

mapfile -t SAMPLES < <(tail -n +2 "$SAMPLESHEET" | cut -f1)
BAMS=(); for s in "${SAMPLES[@]}"; do BAMS+=("$ALN/$s.Aligned.sortedByCoord.out.bam"); done

fc () {   # $1 strand  $2 output
  featureCounts -T "$THREADS" -F SAF -a "$REF/combined.saf" -o "$2" \
     -p --countReadPairs -B -C -s "$1" -Q "$MIN_MAPQ" \
     "${BAMS[@]}" > "$2.log" 2>&1
}

if [[ "$STRAND" == "auto" ]]; then
  for st in 0 1 2; do log "featureCounts -s $st"; fc $st "$CNT/fc_s$st.txt"; done
  # per-sample assigned pairs for each orientation
  python3 - "$CNT" "${SAMPLES[@]}" <<'EOF'
import sys
cnt, samples = sys.argv[1], sys.argv[2:]
def assigned(f):
    for line in open(f):
        if line.startswith("Assigned"):
            return [int(x) for x in line.split("\t")[1:]]
a = {s: assigned(f"{cnt}/fc_s{s}.txt.summary") for s in (0, 1, 2)}
tot1, tot2 = sum(a[1]), sum(a[2])
ratio = (tot1 + 1) / (tot2 + 1)
strand = 1 if ratio > 3 else (2 if ratio < 1/3 else 0)
with open(f"{cnt}/strandedness.tsv", "w") as fh:
    fh.write("sample\tassigned_s0\tassigned_s1\tassigned_s2\tfrac_forward\n")
    for i, s in enumerate(samples):
        fwd = a[1][i] / max(a[1][i] + a[2][i], 1)
        fh.write(f"{s}\t{a[0][i]}\t{a[1][i]}\t{a[2][i]}\t{fwd:.3f}\n")
    fh.write(f"#chosen_strand\t{strand}\tforward/reverse ratio={ratio:.3f}\n")
print(strand)
EOF
  CHOSEN=$(awk '$1=="#chosen_strand"{print $2}' "$CNT/strandedness.tsv")
  log "Inferred strandedness: -s $CHOSEN  (0=unstranded, 1=forward, 2=reverse; see strandedness.tsv)"
  cp "$CNT/fc_s$CHOSEN.txt" "$CNT/featureCounts_counts.txt"
  cp "$CNT/fc_s$CHOSEN.txt.summary" "$CNT/featureCounts_counts.txt.summary"
else
  fc "$STRAND" "$CNT/featureCounts_counts.txt"
  echo -e "#chosen_strand\t$STRAND\tset in config" > "$CNT/strandedness.tsv"
fi

# clean matrix: gene_id + one column per sample
python3 - "$CNT/featureCounts_counts.txt" "$CNT/counts_matrix.tsv" "${SAMPLES[@]}" <<'EOF'
import sys
src, dst, samples = sys.argv[1], sys.argv[2], sys.argv[3:]
with open(src) as fi, open(dst, "w") as fo:
    for line in fi:
        if line.startswith("#"): continue
        f = line.rstrip("\n").split("\t")
        if f[0] == "Geneid":
            fo.write("gene_id\t" + "\t".join(samples) + "\n")
        else:
            fo.write(f[0] + "\t" + "\t".join(f[6:]) + "\n")
EOF
log "Count matrix: $CNT/counts_matrix.tsv"
