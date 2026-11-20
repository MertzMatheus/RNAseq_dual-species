# Dual RNA-seq of *Candida albicans* and *Staphylococcus aureus* biofilms

Analysis pipeline for paired-end RNA-seq of single- and dual-species biofilms of *C. albicans* SC5314 and four *S. aureus* isolates (JE2, CI19, S5, CI5). Reads from both organisms are aligned to a combined reference, counted per gene and analysed per species with DESeq2, followed by GO/KEGG over-representation analysis and publication figures.

## Study design

27 libraries, three biological replicates per group:

| Group | Organisms | Samples |
|---|---|---|
| `Calb_mono` | *C. albicans* | 3 |
| `JE2_mono`, `CI19_mono`, `S5_mono`, `CI5_mono` | *S. aureus* isolate | 3 each |
| `Calb_JE2`, `Calb_CI19`, `Calb_S5`, `Calb_CI5` | *C. albicans* + isolate | 3 each |

Comparisons (log2FC > 0 means higher in the first group):

| ID | Genes | Contrasts |
|---|---|---|
| C1 | *S. aureus* | isolates grown alone: all 6 isolate pairs |
| C2a | *S. aureus* | the same 6 pairs in co-culture with *C. albicans* |
| C2b | *S. aureus* | each isolate with vs without *C. albicans* (4) |
| C3 | *C. albicans* | the 6 pairs of co-cultures, i.e. by partner isolate |
| C4 | *C. albicans* | each co-culture vs *C. albicans* alone (4), plus the mean of the 4 co-cultures |

Derived analyses:

- isolate-specific genes (C1, C2a, C3): differentially expressed in all three contrasts involving an isolate, in the same direction;
- shared and isolate-specific responses (C2b, C4): intersections of the four per-isolate DEG lists, up and down separately;
- C1 vs C2a: concordance of isolate differences with and without *C. albicans*;
- C3 x C4: fungal genes that are specific to one partner isolate and also change relative to *C. albicans* alone;
- expression heatmaps of all DEGs per species, with hierarchical clustering of genes.

## Requirements

- Linux, bash, conda/mamba (or micromamba)
- ~4 GB RAM for STAR (combined genome is ~17 Mb); 16 cores recommended
- Internet access for steps 00, 02b (first run) and 08

```bash
mamba env create -f envs/rnaseq_tools.yml   # fastp, FastQC, SortMeRNA, STAR, SAMtools, Subread, MultiQC
mamba env create -f envs/rnaseq_R.yml       # R, DESeq2, ashr, ggplot2, patchwork, ggraph, VennDiagram
```

## Input layout

```
<PROJECT_DIR>/
├── raw_data/<sample>/<sample>_1.fq.gz, <sample>_2.fq.gz, MD5.txt
└── genomes/
    ├── C_albicans_SC5314.fna         # GCF_000182965.3, Assembly 22 haplotype A
    └── S_aureus_USA300_FPR3757.fna   # GCF_000013465.1
```

GFF3 annotations are downloaded from NCBI RefSeq by step 00 if absent. Samples, groups and folder names are defined in `samplesheet.tsv`. All parameters are in `config.sh`; paths can be set there or through the environment (`PROJECT_DIR`, `RAW_DIR`, `GENOME_DIR`, `OUT_DIR`, `SORTMERNA_DB_DIR`).

## Usage

```bash
export PROJECT_DIR=/path/to/project
nohup bash run_pipeline.sh > run_pipeline.log 2>&1 &
```

Steps can also be run individually; sample-level steps skip samples that are already complete.

| Step | Script | Description | Output (`$OUT_DIR`) |
|---|---|---|---|
| 00 | `00_get_annotations.sh` | download RefSeq GFF3, check FASTA/GFF sequence names | `genomes/*.gff` |
| 01 | `01_check_md5.sh` | verify FASTQ checksums | `qc/md5_check.tsv` |
| 02 | `02_qc_trim.sh` | FastQC, fastp trimming, FastQC | `01_fastqc_raw/`, `02_trimmed/`, `03_fastqc_trimmed/` |
| 02b | `02b_remove_rRNA.sh` | SortMeRNA rRNA removal (SILVA + Rfam) | `02b_rRNA_removed/` |
| 03 | `03_build_reference.sh` | combined genome, species-tagged annotation, STAR index | `reference/` |
| 04 | `04_align.sh` | STAR alignment, per-species read split, mismatch rate | `04_alignment/` |
| 05 | `05_count.sh` | featureCounts, strandedness inferred from the data | `05_counts/counts_matrix.tsv` |
| 06 | `06_qc_report.sh` | per-sample QC table for every step, MultiQC | `qc/QC_summary.xlsx` |
| 07 | `07_differential_expression.sh` | DESeq2 (C1-C4), derived gene sets, figures | `06_differential_expression/` |
| 08 | `08_functional_annotation.sh` | GO (CGD, UniProt) and KEGG (`cal`, `saa`) term tables | `reference/functional/` |
| 09 | `09_enrichment.sh` | over-representation analysis, enrichment figures | `07_functional_enrichment/` |
| 10 | `10_package_results.sh` | methods, QC, figures and tables in one folder + zip | `<PROJECT_NAME>_results/` |

On a SLURM cluster, wrap a step with `sbatch -c 32 --mem 64G --wrap "bash 04_align.sh"`.

## Methods summary

- **Combined reference.** All libraries, including single-species controls, are aligned to *C. albicans* + *S. aureus* together, so each read is assigned to the organism it fits best; reads in the wrong organism in single-species samples are reported as cross-species signal (flagged above 1%).
- **Species-specific normalisation.** The fungus:bacteria read ratio varies between dual-species biofilms, so counts are split by organism and each organism is normalised on its own genes.
- **rRNA.** Removed before alignment with SortMeRNA (the five *S. aureus* rRNA operons would otherwise appear only as multi-mappers); annotated rRNA/tRNA genes are also excluded from DE.
- **Statistics.** DESeq2 Wald tests, ashr-shrunken log2FC reported alongside MLE; DEG = BH-adjusted P < 0.05 and |log2FC| >= 1.
- **Enrichment.** One-sided hypergeometric test per gene set and database (GO BP, GO MF, KEGG), background = tested genes with at least one annotation, terms of 5-500 genes, BH correction; redundant terms (Jaccard >= 0.7) are hidden in figures but kept in tables.

The full methods text, with references, is in [`docs/METHODS.md`](docs/METHODS.md).

## Outputs

`06_differential_expression/`

- `tables/DEG_summary.xlsx`: DEG counts per contrast and per derived set
- `tables/<comparison>_DE.xlsx`: all genes for every contrast (log2FC, shrunken log2FC, P, group means)
- `tables/<comparison>_shared_specific_genes.xlsx`: Venn gene lists
- `tables/Expression_heatmap_clusters.xlsx`: heatmap cluster membership and mean z-score per group
- `figures/`: PCA, sample correlation, expression heatmaps, volcano plots, Venn and UpSet diagrams, concordance plot

`07_functional_enrichment/`

- `tables/<comparison>_ORA.xlsx`, `ORA_summary.xlsx`, `ORA_all_tested_terms.csv`
- `figures/<comparison>/<comparison>_ORA_summary`: all enriched terms of a comparison in one dot matrix
- `figures/<comparison>/<contrast>_ORA`: dot plot and gene-term network per contrast

Figures are written as PDF (vector, embedded fonts), 600-dpi TIFF and PNG, using Arial (or Liberation Sans) at 7 pt and 89/183 mm widths. Colours are defined in `scripts/fig_style.R`.

## Caveats

- All *S. aureus* isolates are aligned to USA300_FPR3757. Genes absent from that reference are not quantified, and genes absent from an isolate look down-regulated; the `zero_reads_in` column flags genes with no reads in every replicate of a group, and the per-isolate mismatch rate is reported in the QC table.
- Check the number of *S. aureus* read pairs in the dual-species libraries (`QC_summary.xlsx`); fewer than ~2 M non-rRNA pairs limits power for C2a and C2b.

## Testing

`test/simulate_test_data.py` builds a small project with known fold changes, isolate divergence, a multi-lane sample and a contaminated sample; `test/make_test_functional_annotation.py` writes offline GO/KEGG files with a planted pathway so steps 08-09 run without internet.

```bash
python3 test/simulate_test_data.py /tmp/test_project
PROJECT_DIR=/tmp/test_project bash run_pipeline.sh   # run make_test_functional_annotation.py before step 08
```

`test/example_output/` contains figures produced from simulated data.

## Software

fastp, FastQC, SortMeRNA, STAR, SAMtools, Subread (featureCounts), MultiQC, R/Bioconductor (DESeq2, ashr, ggplot2, patchwork, ggraph, igraph, VennDiagram, openxlsx). Versions are pinned in `envs/`; `sessionInfo()` is written next to the R outputs.
