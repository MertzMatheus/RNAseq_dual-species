# Methods – RNA-seq data analysis

## Read processing and quality control

Raw paired-end reads from the 27 libraries were checked for integrity against the MD5 checksums provided by the sequencing facility, and read quality was assessed with FastQC v0.12.1 (Andrews, 2010) before and after processing. Adapter sequences were detected and removed with fastp v0.23.4 (Chen et al., 2018), which also trimmed poly-G tails and low-quality 3′ ends using a 4-nt sliding window with a mean quality threshold of Q20; read pairs in which either read was shorter than 36 nt or contained more than five ambiguous bases were discarded.

Residual ribosomal RNA reads were removed in silico with SortMeRNA v4.3.7 (Kopylova et al., 2012) against its default v4.3 database, which comprises bacterial, archaeal and eukaryotic small- and large-subunit rRNA sequences from SILVA (Quast et al., 2013) and 5S and 5.8S rRNA sequences from Rfam (Kalvari et al., 2021). A read pair was discarded when either mate matched an rRNA sequence, and the fraction removed was recorded for each library as a measure of rRNA-depletion efficiency. Quality metrics from every processing step were aggregated with MultiQC v1.35 (Ewels et al., 2016) and are reported in Supplementary Table [x].

## Read alignment and gene quantification

Because dual-species libraries contain transcripts from both organisms, all reads were aligned to a single composite reference comprising the *Candida albicans* SC5314 genome (Assembly 22, haplotype A; RefSeq GCF_000182965.3; Skrzypek et al., 2017) and the *Staphylococcus aureus* USA300_FPR3757 chromosome and plasmids pUSA01–pUSA03 (RefSeq GCF_000013465.1; Diep et al., 2006), the latter serving as the reference for all four *S. aureus* isolates. Mapping every library, including single-species controls, to this composite genome allowed each read to be assigned to the organism it matched best and provided a direct measure of cross-species signal. Alignment was performed with STAR v2.7.11b (Dobin et al., 2013) using an index supplemented with annotated *C. albicans* splice junctions; maximum intron length and mate gap were restricted to 5,000 nt to prevent spurious splicing within the bacterial genome, and reads aligning to more than 10 loci were discarded. Uniquely aligned read pairs on each organism's sequences, together with the per-base mismatch rate of each *S. aureus* isolate relative to USA300_FPR3757, were obtained with SAMtools v1.20 (Danecek et al., 2021).

Read pairs were assigned to genes with featureCounts from Subread v2.0.6 (Liao et al., 2014) using the RefSeq GFF3 annotations, summarising exons per gene for *C. albicans* and counting over annotated gene coordinates for *S. aureus*. Only uniquely aligned pairs (MAPQ ≥ 10) with both mates on the same reference sequence were counted, and pairs overlapping more than one gene were excluded. Library strandedness was inferred empirically by comparing assignment rates obtained under unstranded, forward and reverse settings; the libraries were [reverse-stranded], and the corresponding counts were used for all downstream analyses.

## Differential expression analysis

Differential expression was analysed in R v4.3.3 (R Core Team, 2023) with DESeq2 v1.42.0 (Love et al., 2014). Because the relative contribution of fungal and bacterial transcripts differs between dual-species biofilms, the count matrix was partitioned by organism and each organism was analysed with its own genes only, so that size factors were estimated within species. Residual counts on annotated rRNA and tRNA genes were removed, and genes with at least 10 reads in at least three samples were retained. Each comparison was fitted as a separate model with experimental group as the single factor (design ~ group), restricted to the samples involved, and contrasts were evaluated with Wald tests.

Differences between *S. aureus* isolates were assessed with all six pairwise contrasts (CI19 versus JE2, S5 versus JE2, CI5 versus JE2, CI19 versus S5, CI19 versus CI5 and S5 versus CI5), using *S. aureus* genes from the single-species biofilms (comparison 1) and, separately, from the dual-species biofilms (comparison 2a). The response of each isolate to *C. albicans* was evaluated by comparing the dual-species biofilm with the single-species biofilm of the same isolate (comparison 2b). For the fungus, *C. albicans* genes were compared between the four dual-species biofilms using the same six pairwise contrasts (comparison 3), and each dual-species biofilm was compared with *C. albicans* grown alone (comparison 4); the average response of *C. albicans* to the four isolates was additionally estimated with a numeric contrast of the mean of the four dual-species groups against *C. albicans* alone.

P values were adjusted for multiple testing with the Benjamini–Hochberg procedure (Benjamini and Hochberg, 1995), and genes with an adjusted P < 0.05 and an absolute log2 fold change ≥ 1 were considered differentially expressed. Log2 fold changes shrunk with the adaptive shrinkage estimator ashr v2.2-63 (Stephens, 2017) are reported alongside the maximum-likelihood estimates. Because a single reference genome was used for all *S. aureus* isolates, genes without reads in any replicate of a group were flagged as potentially absent from that isolate rather than transcriptionally repressed.

## Shared and isolate-specific genes

For the pairwise comparisons (1, 2a and 3), each contrast was oriented relative to a focal isolate, and a gene was defined as isolate-specific when it was differentially expressed in all three contrasts involving that isolate in the same direction, that is, consistently higher or lower than in every other isolate. The three contrasts of each isolate were displayed as Venn diagrams, and the overlaps among the six pairwise sets of differentially expressed genes as UpSet plots (Lex et al., 2014). For comparisons 2b and 4, the up- and down-regulated genes of the four isolate-specific contrasts were intersected separately: genes differentially expressed with all four isolates were defined as the shared response, and genes differentially expressed with only one isolate as isolate-specific responses.

To determine whether differences between isolates persist in the presence of *C. albicans*, the log2 fold changes of each pair of isolates in single- and dual-species biofilms were compared by Spearman correlation, and genes were classified according to whether they were differentially expressed in one or both conditions. Finally, *C. albicans* genes that were isolate-specific in comparison 3 and differentially expressed in the same direction relative to *C. albicans* alone in comparison 4 were defined as isolate-specific fungal responses.

## Functional enrichment analysis

Gene Ontology (GO; Gene Ontology Consortium, 2023) and KEGG pathway (Kanehisa and Goto, 2000) annotations were retrieved on [date]. For *C. albicans*, GO annotations curated by the Candida Genome Database (Skrzypek et al., 2017) were combined with UniProtKB annotations for strain SC5314 (UniProt Consortium, 2023), and KEGG pathways were obtained for organism code cal. For *S. aureus*, GO annotations were obtained from UniProtKB for strain USA300, and KEGG pathways for USA300_FPR3757 (organism code saa) were linked to RefSeq genes through their original locus tags. Identifiers were reconciled through locus tags, CGD identifiers, Assembly 22 systematic names and orf19 aliases, retaining only identifiers that mapped to a single gene. Annotations carrying a NOT qualifier were discarded, all evidence codes were retained, and GO annotations were propagated to ancestral terms through the is_a and part_of relationships of the GO ontology (release [version]); root terms and KEGG global overview maps were excluded. GO Biological Process annotations covered [X]% of *C. albicans* and [Y]% of *S. aureus* genes.

Over-representation analysis was performed for the up- and down-regulated genes of every contrast, for the isolate-specific gene sets, for the shared and isolate-specific responses and for the isolate-specific fungal responses. The background of each gene set comprised the genes tested in the corresponding contrast or comparison that carried at least one annotation in the database under consideration. Enrichment of GO Biological Process and Molecular Function terms and of KEGG pathways containing 5–500 background genes was assessed with a one-sided hypergeometric test, P values were adjusted with the Benjamini–Hochberg procedure, and terms with an adjusted P < 0.05 and at least two genes from the set were considered enriched. To reduce redundancy in figures, a significant term was omitted when its enriched genes overlapped by at least 70% (Jaccard index) with those of a more significant term. For each contrast, the five most significant non-redundant terms per database and direction were displayed as a dot plot (gene ratio, number of genes and adjusted P) together with a gene–term network of the five most significant terms. For each comparison, the enriched terms of all its contrasts and gene sets were summarised in a single dot matrix, and the same analysis was applied to the gene clusters defined in the expression heatmaps.

## Data visualisation

For principal component analysis, counts were transformed with the DESeq2 variance-stabilising transformation and the 500 most variable genes of each organism were used; the same genes were used to compute Pearson correlations between samples. In volcano plots, the five up-regulated and five down-regulated DEGs with the largest absolute fold changes were labelled. To compare overall expression between groups, the variance-stabilised expression of all genes that were differentially expressed in at least one contrast (*S. aureus*: C1, C2a and C2b; *C. albicans*: C3 and C4) was centred and scaled per gene (z-score) and shown as heatmaps; genes were clustered hierarchically (Euclidean distance, Ward's criterion; Ward, 1963) and the dendrogram was cut into six expression clusters. Figures were produced in R with ggplot2 (Wickham, 2016) and patchwork (Pedersen, 2024b), gene–term networks with ggraph (Pedersen, 2024a) and igraph (Csardi and Nepusz, 2006), and Venn diagrams with VennDiagram (Chen and Boutros, 2011). A fixed palette was used throughout, in which each *S. aureus* isolate has the same colour in all figures and single- and dual-species biofilms are distinguished by symbol.

## Data and code availability

Raw sequencing reads have been deposited in the NCBI Sequence Read Archive under BioProject [accession]. All analysis scripts are available at [repository URL].

## References

Andrews S (2010). FastQC: a quality control tool for high throughput sequence data. Babraham Bioinformatics. https://www.bioinformatics.babraham.ac.uk/projects/fastqc/

Benjamini Y, Hochberg Y (1995). Controlling the false discovery rate: a practical and powerful approach to multiple testing. *J R Stat Soc Series B* 57, 289–300.

Chen H, Boutros PC (2011). VennDiagram: a package for the generation of highly-customizable Venn and Euler diagrams in R. *BMC Bioinformatics* 12, 35.

Chen S, Zhou Y, Chen Y, Gu J (2018). fastp: an ultra-fast all-in-one FASTQ preprocessor. *Bioinformatics* 34, i884–i890.

Csardi G, Nepusz T (2006). The igraph software package for complex network research. *InterJournal Complex Systems* 1695.

Danecek P, Bonfield JK, Liddle J, et al. (2021). Twelve years of SAMtools and BCFtools. *GigaScience* 10, giab008.

Diep BA, Gill SR, Chang RF, et al. (2006). Complete genome sequence of USA300, an epidemic clone of community-acquired meticillin-resistant *Staphylococcus aureus*. *Lancet* 367, 731–739.

Dobin A, Davis CA, Schlesinger F, et al. (2013). STAR: ultrafast universal RNA-seq aligner. *Bioinformatics* 29, 15–21.

Ewels P, Magnusson M, Lundin S, Käller M (2016). MultiQC: summarize analysis results for multiple tools and samples in a single report. *Bioinformatics* 32, 3047–3048.

Gene Ontology Consortium; Aleksander SA, Balhoff J, Carbon S, et al. (2023). The Gene Ontology knowledgebase in 2023. *Genetics* 224, iyad031.

Kalvari I, Nawrocki EP, Ontiveros-Palacios N, et al. (2021). Rfam 14: expanded coverage of metagenomic, viral and microRNA families. *Nucleic Acids Res* 49, D192–D200.

Kanehisa M, Goto S (2000). KEGG: Kyoto Encyclopedia of Genes and Genomes. *Nucleic Acids Res* 28, 27–30.

Kopylova E, Noé L, Touzet H (2012). SortMeRNA: fast and accurate filtering of ribosomal RNAs in metatranscriptomic data. *Bioinformatics* 28, 3211–3217.

Lex A, Gehlenborg N, Strobelt H, Vuillemot R, Pfister H (2014). UpSet: visualization of intersecting sets. *IEEE Trans Vis Comput Graph* 20, 1983–1992.

Liao Y, Smyth GK, Shi W (2014). featureCounts: an efficient general purpose program for assigning sequence reads to genomic features. *Bioinformatics* 30, 923–930.

Love MI, Huber W, Anders S (2014). Moderated estimation of fold change and dispersion for RNA-seq data with DESeq2. *Genome Biol* 15, 550.

Pedersen TL (2024a). ggraph: an implementation of grammar of graphics for graphs and networks. R package. https://CRAN.R-project.org/package=ggraph

Pedersen TL (2024b). patchwork: the composer of plots. R package. https://CRAN.R-project.org/package=patchwork

Quast C, Pruesse E, Yilmaz P, et al. (2013). The SILVA ribosomal RNA gene database project: improved data processing and web-based tools. *Nucleic Acids Res* 41, D590–D596.

R Core Team (2023). R: a language and environment for statistical computing. R Foundation for Statistical Computing, Vienna, Austria.

Skrzypek MS, Binkley J, Binkley G, Miyasato SR, Simison M, Sherlock G (2017). The Candida Genome Database (CGD): incorporation of Assembly 22, systematic identifiers and visualization of high throughput sequencing data. *Nucleic Acids Res* 45, D592–D596.

Stephens M (2017). False discovery rates: a new deal. *Biostatistics* 18, 275–294.

UniProt Consortium (2023). UniProt: the Universal Protein Knowledgebase in 2023. *Nucleic Acids Res* 51, D523–D531.

Ward JH (1963). Hierarchical grouping to optimize an objective function. *Journal of the American Statistical Association* 58, 236–244.

Wickham H (2016). *ggplot2: Elegant Graphics for Data Analysis*, 2nd edn. Springer, New York.
