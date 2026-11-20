#!/usr/bin/env bash
# Run all steps in order. Finished samples are skipped, so the script can be relaunched safely.
#   nohup bash run_pipeline.sh > run_pipeline.log 2>&1 &
set -euo pipefail
cd "$(dirname "$0")"
for step in 00_get_annotations.sh 01_check_md5.sh 02_qc_trim.sh 02b_remove_rRNA.sh 03_build_reference.sh \
            04_align.sh 05_count.sh 06_qc_report.sh 07_differential_expression.sh \
            08_functional_annotation.sh 09_enrichment.sh 10_package_results.sh; do
  echo "== $step"
  bash "$step"
done
