#!/usr/bin/env Rscript
# Differential expression for dual-species biofilms (C. albicans SC5314 x S. aureus isolates).
#
#   C1   S. aureus genes    isolates grown alone            6 pairwise isolate contrasts
#   C2a  S. aureus genes    isolates with C. albicans       the same 6 pairs in co-culture
#   C2b  S. aureus genes    each isolate with vs without C. albicans (4 contrasts)
#   C3   C. albicans genes  co-cultures                     the same 6 pairs, by partner isolate
#   C4   C. albicans genes  each co-culture vs C. albicans alone (4) + mean of the 4 co-cultures vs alone
#
# Isolate-specific genes (C1, C2a, C3): DE in all three contrasts involving an isolate, same direction.
# Shared/specific responses (C2b, C4): intersections of the four per-isolate DEG lists, up and down separately.
# Also: C1 vs C2a log2FC concordance, C3 x C4 isolate-specific fungal responses, PCA, sample correlation and
# expression heatmaps. Each species is normalised on its own genes. Gene sets for ORA go to ora_input.rds.
suppressPackageStartupMessages({
  library(optparse); library(DESeq2); library(ggplot2); library(ggrepel); library(openxlsx); library(patchwork)
})
.script_dir <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])))
source(file.path(.script_dir, "fig_style.R"))

opt <- parse_args(OptionParser(option_list = list(
  make_option("--counts"), make_option("--annot"), make_option("--samples"), make_option("--outdir"),
  make_option("--padj", type = "double", default = 0.05), make_option("--lfc", type = "double", default = 1),
  make_option("--min_count", type = "integer", default = 10), make_option("--min_samples", type = "integer", default = 3),
  make_option("--drop_rrna", default = "yes"))))
PADJ <- opt$padj; LFC <- opt$lfc; set.seed(1)
TAB <- file.path(opt$outdir, "tables"); FIG <- file.path(opt$outdir, "figures")
dir.create(TAB, recursive = TRUE, showWarnings = FALSE); dir.create(FIG, recursive = TRUE, showWarnings = FALSE)

# ------------------------------------------------------------------ data ----
counts <- as.matrix(read.delim(opt$counts, row.names = 1, check.names = FALSE))
annot <- read.delim(opt$annot, quote = "", colClasses = "character", na.strings = ""); annot[is.na(annot)] <- ""
rownames(annot) <- annot$gene_id
ss <- read.delim(opt$samples, colClasses = "character"); rownames(ss) <- ss$sample
counts <- counts[, ss$sample]; annot <- annot[rownames(counts), ]
use_gene <- if (tolower(opt$drop_rrna) == "yes") !grepl("rRNA|tRNA", annot$biotype, ignore.case = TRUE) else rep(TRUE, nrow(annot))
ANN <- c("gene_id", "name", "old_locus_tag", "product")
gene_label <- function(ids) ifelse(annot[ids, "name"] != "", annot[ids, "name"],
                                   ifelse(annot[ids, "old_locus_tag"] != "", annot[ids, "old_locus_tag"], ids))
SPN <- c(Calb = "C. albicans", Saur = "S. aureus")
PAIRS <- list(c("CI19", "JE2"), c("S5", "JE2"), c("CI5", "JE2"), c("CI19", "S5"), c("CI19", "CI5"), c("S5", "CI5"))

# ------------------------------------------------------------ DESeq2 core ----
fit <- function(sp, groups) {
  smp <- ss$sample[ss$group %in% groups]
  g <- annot$gene_id[annot$species == sp & use_gene]
  cts <- round(counts[g, smp, drop = FALSE])
  cd <- ss[smp, ]; cd$group <- factor(cd$group, levels = groups)
  keep <- rowSums(cts >= opt$min_count) >= opt$min_samples
  cat(sprintf("  %s, %d samples: %d/%d genes kept\n", SPN[sp], length(smp), sum(keep), length(keep)))
  DESeq(DESeqDataSetFromMatrix(cts[keep, ], cd, ~group), quiet = TRUE)
}
zero_in <- function(dds) {
  raw <- counts(dds); z <- sapply(levels(dds$group), function(g) rowSums(raw[, dds$group == g, drop = FALSE]) == 0)
  apply(z, 1, function(x) paste(colnames(z)[x], collapse = ","))
}
contrast_table <- function(dds, contrast, groups_shown, sp) {
  r <- results(dds, contrast = contrast, alpha = PADJ)
  sh <- suppressMessages(lfcShrink(dds, contrast = contrast, res = r, type = "ashr", quiet = TRUE))
  nc <- counts(dds, normalized = TRUE)
  gm <- sapply(groups_shown, function(g) rowMeans(nc[, dds$group == g, drop = FALSE])); colnames(gm) <- paste0("mean_norm_", groups_shown)
  d <- data.frame(annot[rownames(r), ANN], baseMean = r$baseMean, log2FoldChange = r$log2FoldChange, lfcSE = r$lfcSE,
                  log2FC_shrunk = sh$log2FoldChange, pvalue = r$pvalue, padj = r$padj, check.names = FALSE)
  d$regulation <- ifelse(!is.na(d$padj) & d$padj < PADJ & abs(d$log2FoldChange) >= LFC, ifelse(d$log2FoldChange > 0, "Up", "Down"), "NS")
  d <- cbind(d, round(gm[rownames(r), , drop = FALSE], 2))
  if (sp == "Saur") d$zero_reads_in <- zero_in(dds)[rownames(r)]
  d[order(d$padj, -abs(d$log2FoldChange), na.last = TRUE), ]
}

# --------------------------------------------------------------- figures ----
volcano <- function(t, title) {
  d <- t[!is.na(t$padj), ]; d$y <- -log10(pmax(d$padj, 1e-300))
  d$regulation <- factor(d$regulation, levels = c("NS", "Down", "Up")); d <- d[order(d$regulation), ]
  up <- d[d$regulation == "Up", ]; dn <- d[d$regulation == "Down", ]
  top <- rbind(head(up[order(-up$log2FoldChange), ], 5), head(dn[order(dn$log2FoldChange), ], 5))   # largest fold changes
  nu <- sum(d$regulation == "Up"); nd <- sum(d$regulation == "Down")
  ggplot(d, aes(log2FoldChange, y)) +
    geom_vline(xintercept = c(-LFC, LFC), linetype = "22", linewidth = lw(0.4), colour = "grey55") +
    geom_hline(yintercept = -log10(PADJ), linetype = "22", linewidth = lw(0.4), colour = "grey55") +
    geom_point(aes(colour = regulation), size = 0.5, stroke = 0, alpha = 0.85) +
    scale_colour_manual(values = REG_COL, guide = "none") +
    geom_text_repel(data = top, aes(label = gene_label(gene_id)), size = fs(TINY), family = FONT, fontface = "italic",
                    segment.size = lw(0.3), min.segment.length = 0, max.overlaps = Inf, box.padding = 0.2) +
    labs(title = title, subtitle = sprintf("%d up, %d down", nu, nd),
         x = expression(log[2]~fold~change), y = expression(-log[10]~adjusted~italic(P))) + theme_pub()
}
volcano_panel <- function(tabs, labels, file, ncol = 3) {
  ps <- lapply(names(tabs), function(n) volcano(tabs[[n]], labels[[n]]))
  nr <- ceiling(length(ps) / ncol)
  save_fig(wrap_plots(ps, ncol = ncol) + tags_pub(), file, W2, 58 * nr)
}

# --------------------------------------------------------- comparisons ----
run_pairwise <- function(id, sp, grp, lab_fun) {
  groups <- unname(sapply(STRAINS, grp)); dds <- fit(sp, groups)
  tabs <- list(); labels <- list()
  for (p in PAIRS) {
    tag <- paste0(p[1], "_vs_", p[2])
    tabs[[tag]] <- contrast_table(dds, c("group", grp(p[1]), grp(p[2])), c(grp(p[1]), grp(p[2])), sp)
    labels[[tag]] <- lab_fun(p)
  }
  list(id = id, sp = sp, dds = dds, tabs = tabs, labels = labels, type = "pairwise")
}
run_response <- function(id, sp, groups, base_grp, grp, lab_fun, extra = FALSE) {
  dds <- fit(sp, groups); tabs <- list(); labels <- list()
  for (s in STRAINS) {
    tag <- paste0(s, "_withPartner_vs_alone")
    tabs[[tag]] <- contrast_table(dds, c("group", grp(s), base_grp(s)), c(base_grp(s), grp(s)), sp)
    labels[[tag]] <- lab_fun(s)
  }
  out <- list(id = id, sp = sp, dds = dds, tabs = tabs, labels = labels, type = "response")
  if (extra) {   # mean of the four co-cultures vs C. albicans alone (numeric contrast on the group coefficients)
    rn <- resultsNames(dds); w <- ifelse(grepl("^group_", rn), 1 / length(STRAINS), 0)
    out$extra <- contrast_table(dds, w, groups, sp)
  }
  out
}
cat("\n== Fitting models ==\n")
R <- list()
R$C1  <- run_pairwise("C1", "Saur", function(s) paste0(s, "_mono"), function(p) paste(p[1], "vs", p[2]))
R$C2a <- run_pairwise("C2a", "Saur", function(s) paste0("Calb_", s), function(p) paste(p[1], "vs", p[2]))
R$C2b <- run_response("C2b", "Saur", c(paste0(STRAINS, "_mono"), paste0("Calb_", STRAINS)),
                      function(s) paste0(s, "_mono"), function(s) paste0("Calb_", s), function(s) paste(s, "+ Ca vs", s))
R$C3  <- run_pairwise("C3", "Calb", function(s) paste0("Calb_", s), function(p) paste0("Ca+", p[1], " vs Ca+", p[2]))
R$C4  <- run_response("C4", "Calb", c("Calb_mono", paste0("Calb_", STRAINS)),
                      function(s) "Calb_mono", function(s) paste0("Calb_", s), function(s) paste0("Ca+", s, " vs Ca"), extra = TRUE)
FOLDER <- c(C1 = "C1_Saur_isolates_alone", C2a = "C2a_Saur_isolates_with_Calb", C2b = "C2b_Saur_response_to_Calb",
            C3 = "C3_Calb_by_partner_isolate", C4 = "C4_Calb_response_to_each_isolate")
TITLE <- c(C1 = "'C1'~~italic('S. aureus')~'isolates, alone'",
           C2a = "'C2a'~~italic('S. aureus')~'isolates, with'~italic('C. albicans')",
           C2b = "'C2b'~~italic('S. aureus')~'with vs without'~italic('C. albicans')",
           C3 = "'C3'~~italic('C. albicans')~'by partner isolate'",
           C4 = "'C4'~~italic('C. albicans')~'with each isolate vs alone'")

# ------------------------------------------------------- write DE tables ----
book <- function(sheets, file) {
  wb <- createWorkbook(); hs <- createStyle(textDecoration = "bold", fgFill = "#E7E6E6", wrapText = TRUE)
  for (n in names(sheets)) {
    sn <- substr(gsub("[\\[\\]\\*\\?/\\\\:]", "_", n), 1, 31); addWorksheet(wb, sn)
    writeData(wb, sn, sheets[[n]], headerStyle = hs); freezePane(wb, sn, firstRow = TRUE, firstCol = TRUE)
    setColWidths(wb, sn, cols = seq_len(max(1, ncol(sheets[[n]]))), widths = 14)
    if ("regulation" %in% names(sheets[[n]])) {
      ci <- which(names(sheets[[n]]) == "regulation"); nr <- nrow(sheets[[n]]) + 1
      conditionalFormatting(wb, sn, cols = ci, rows = 2:nr, type = "contains", rule = "Up", style = createStyle(bgFill = "#F4CCCC"))
      conditionalFormatting(wb, sn, cols = ci, rows = 2:nr, type = "contains", rule = "Down", style = createStyle(bgFill = "#CFE2F3"))
    }
  }
  saveWorkbook(wb, file, overwrite = TRUE)
}
SUMMARY <- list()
for (k in names(R)) {
  x <- R[[k]]; tabs <- x$tabs; if (!is.null(x$extra)) tabs[["Mean_of_4_vs_alone"]] <- x$extra
  labs <- c(x$labels, if (!is.null(x$extra)) list(Mean_of_4_vs_alone = "Ca+Sa (mean of 4) vs Ca"))
  s <- do.call(rbind, lapply(names(tabs), function(n) data.frame(Comparison = k, Contrast = labs[[n]], Sheet = n,
        Genes_tested = sum(!is.na(tabs[[n]]$padj)), Up = sum(tabs[[n]]$regulation == "Up"), Down = sum(tabs[[n]]$regulation == "Down"))))
  SUMMARY[[k]] <- s
  book(c(list(Summary = s), tabs), file.path(TAB, paste0(FOLDER[k], "_DE.xlsx")))
  volcano_panel(x$tabs, x$labels, file.path(FIG, k, paste0(k, "_volcano")), ncol = if (x$type == "pairwise") 3 else 2)
  if (!is.null(x$extra)) save_fig(volcano(x$extra, "Ca+Sa (mean of 4 co-cultures) vs Ca"), file.path(FIG, k, paste0(k, "_extra_volcano")), W1, 70)
}

# ----------------------------------------- isolate-specific genes (pairwise) ----
oriented <- function(x, focal) {   # the three contrasts of one isolate, oriented so positive = higher in the focal isolate
  out <- list()
  for (p in PAIRS) if (focal %in% p) {
    t <- x$tabs[[paste0(p[1], "_vs_", p[2])]]; other <- setdiff(p, focal); sg <- if (p[1] == focal) 1 else -1
    reg <- t$regulation; if (sg < 0) reg <- ifelse(reg == "Up", "Down", ifelse(reg == "Down", "Up", reg))
    out[[other]] <- data.frame(gene_id = t$gene_id, lfc = sg * t$log2FoldChange, padj = t$padj, reg = reg)
  }
  out
}
SETS <- list()   # gene sets passed to enrichment.R
add_set <- function(id, comparison, kind, label, direction, genes, universe, sp, lfc, group = NA) {
  SETS[[id]] <<- list(id = id, comparison = comparison, kind = kind, label = label, direction = direction, group = group,
                      genes = unique(genes), universe = universe, species = sp, lfc = lfc)
}
isolate_specific <- function(x) {
  k <- x$id; universe <- rownames(x$dds); rows <- list(); venns <- list(); spec <- list()
  for (s in STRAINS) {
    o <- oriented(x, s); others <- names(o)
    for (dir in c("Up", "Down")) {
      sets <- setNames(lapply(o, function(d) d$gene_id[d$reg == dir]), others)
      core <- Reduce(intersect, sets); spec[[paste(s, dir)]] <- core
      hl <- if (dir == "Up") "Higher" else "Lower"
      ttl <- if (k == "C3") paste0(hl, " with ", s, " than with:") else paste0(hl, " in ", s, " than in:")
      sub_ <- sprintf("%s-specific (centre): %d genes", s, length(core))
      venns[[paste(s, dir)]] <- wrap_elements(full = venn_grob(sets, STRAIN_COL[others], ttl, STRAIN_COL[s], sub_))
      if (length(core)) {
        m <- do.call(cbind, lapply(others, function(y) { d <- o[[y]][match(core, o[[y]]$gene_id), ]
                 setNames(data.frame(round(d$lfc, 3), signif(d$padj, 3)), paste0(c("log2FC_vs_", "padj_vs_"), y)) }))
        rows[[paste(s, dir)]] <- data.frame(Isolate = s, Class = paste0(s, "-specific (", tolower(if (dir == "Up") "higher" else "lower"), ")"),
                                            annot[core, ANN], m, mean_log2FC = round(rowMeans(m[, grep("^log2FC", names(m))]), 3), check.names = FALSE)
      }
      lf <- setNames(rowMeans(sapply(o, function(d) d$lfc[match(universe, d$gene_id)]), na.rm = TRUE), universe)
      add_set(paste0(k, "|", s, "_specific_", dir), k, "isolate_specific", paste0(s, "-specific ", if (dir == "Up") "(higher)" else "(lower)"),
              dir, core, universe, x$sp, lf, s)
    }
  }
  what <- if (k == "C3") "C. albicans genes compared between co-cultures" else if (k == "C2a") "S. aureus genes, isolates in co-culture with C. albicans" else "S. aureus genes, isolates grown alone"
  cap <- paste0(what, ". Each diagram belongs to one isolate (title colour); each circle holds the genes differentially expressed between that isolate\n",
                "and one other isolate (circle colour). Top row: genes higher in the isolate; bottom row: genes lower. The centre of a diagram holds the\n",
                "isolate-specific genes, i.e. higher (or lower) than in all three other isolates; their number is given under each title.\n",
                "Gene lists: ", FOLDER[k], "_shared_specific_genes.xlsx (one sheet per isolate).")
  fig <- wrap_plots(venns[paste(STRAINS, "Up")], nrow = 1) / wrap_plots(venns[paste(STRAINS, "Down")], nrow = 1) +
    fig_title(QUESTION[k], cap) + plot_annotation(tag_levels = "a") & theme(plot.tag = element_text(size = 8, face = "bold", family = FONT))
  save_fig(fig, file.path(FIG, k, paste0(k, "_venn_isolate_specific")), W2, 130)
  pair_sets <- setNames(lapply(x$tabs, function(t) t$gene_id[t$regulation != "NS"]), unlist(x$labels))
  up <- upset_pub(pair_sets, ylab = "DEGs")
  if (!is.null(up)) save_fig(up + fig_title(QUESTION[k], paste0("Overlap of the DEGs of the six isolate pairs. Top bars: genes in each combination of pairs ",
      "(filled dots below); left bars: DEGs per pair.\nA gene shared by the three pairs of one isolate is usually specific to that isolate (see the Venn diagrams).")),
      file.path(FIG, k, paste0(k, "_upset_pairs")), W2, 90)
  # every gene that is DE in at least one pair
  ids <- unique(unlist(pair_sets))
  allg <- data.frame(annot[ids, ANN], check.names = FALSE)
  for (n in names(x$tabs)) { t <- x$tabs[[n]][match(ids, x$tabs[[n]]$gene_id), ]
    allg[[paste0("reg_", n)]] <- t$regulation; allg[[paste0("log2FC_", n)]] <- round(t$log2FoldChange, 3) }
  allg$n_pairs_DE <- rowSums(allg[, grep("^reg_", names(allg))] != "NS", na.rm = TRUE)
  cls <- rep("", length(ids)); for (n in names(spec)) cls[ids %in% spec[[n]]] <- paste0(cls[ids %in% spec[[n]]], sub(" Up", "-specific (higher)", sub(" Down", "-specific (lower)", n)), "; ")
  allg$isolate_specific <- sub("; $", "", cls); allg <- allg[order(-allg$n_pairs_DE), ]
  summ <- data.frame(Isolate = STRAINS, Specific_higher = sapply(STRAINS, function(s) length(spec[[paste(s, "Up")]])),
                     Specific_lower = sapply(STRAINS, function(s) length(spec[[paste(s, "Down")]])))
  per_iso <- lapply(setNames(STRAINS, STRAINS), function(s) { r <- rows[intersect(paste(s, c("Up", "Down")), names(rows))]
    if (length(r)) do.call(rbind, r) else data.frame(Note = paste("No", s, "-specific genes")) })
  list(summary = summ, per_isolate = per_iso, all = allg, spec = spec)
}

# ------------------------------------ shared / specific responses (C2b, C4) ----
venn_regions <- function(sets) {   # region label of each gene in a 4-set Venn
  ids <- unique(unlist(sets)); M <- sapply(sets, function(s) ids %in% s); if (is.null(dim(M))) M <- matrix(M, nrow = 1)
  lab <- apply(M, 1, function(r) { n <- names(sets)[r]
    if (length(n) == length(sets)) "Shared by all 4" else if (length(n) == 1) paste("Specific to", n) else paste("Shared by", paste(n, collapse = " + ")) })
  data.frame(gene_id = ids, region = lab, n_isolates = rowSums(M))
}
response_sets <- function(x) {
  k <- x$id; universe <- rownames(x$dds); venns <- list(); regs <- list()
  lfcm <- sapply(x$tabs, function(t) t$log2FoldChange[match(universe, t$gene_id)]); colnames(lfcm) <- STRAINS
  for (dir in c("Up", "Down")) {
    sets <- setNames(lapply(x$tabs, function(t) t$gene_id[t$regulation == dir]), STRAINS)
    ttl <- if (dir == "Up") "Up-regulated" else "Down-regulated"
    rg <- venn_regions(sets); rg$direction <- dir; regs[[dir]] <- rg
    ttl2 <- if (k == "C2b") paste0(ttl, " with C. albicans") else paste0(ttl, " with S. aureus")
    sub_ <- sprintf("Shared by all 4 (centre): %d  |  specific to one isolate: %d", sum(rg$region == "Shared by all 4"), sum(rg$n_isolates == 1))
    venns[[dir]] <- wrap_elements(full = venn_grob(sets, STRAIN_COL[STRAINS], ttl2, "black", sub_))
    shared <- rg$gene_id[rg$region == "Shared by all 4"]
    add_set(paste0(k, "|shared_", dir), k, "shared_response", paste("Shared by all 4", tolower(ttl)), dir, shared, universe, x$sp,
            setNames(rowMeans(lfcm, na.rm = TRUE), universe))
    for (s in STRAINS) add_set(paste0(k, "|", s, "_specific_", dir), k, "specific_response", paste0(s, "-specific, ", tolower(ttl)), dir,
                               rg$gene_id[rg$region == paste("Specific to", s)], universe, x$sp, setNames(lfcm[, s], universe), s)
  }
  cap <- if (k == "C2b") paste0("Each circle = one S. aureus isolate: genes changed in that isolate when grown with C. albicans (vs alone).\n",
      "Centre = response shared by all four isolates; regions inside one circle only = isolate-specific responses.\nGene lists: ", FOLDER[k], "_shared_specific_genes.xlsx.")
    else paste0("Each circle = one partner isolate: C. albicans genes changed in co-culture with that isolate (vs C. albicans alone).\n",
      "Centre = response to all four isolates; regions inside one circle only = response specific to one isolate.\nGene lists: ", FOLDER[k], "_shared_specific_genes.xlsx.")
  save_fig(wrap_plots(venns, nrow = 1) + fig_title(QUESTION[k], cap) + plot_annotation(tag_levels = "a") &
             theme(plot.tag = element_text(size = 8, face = "bold", family = FONT)),
           file.path(FIG, k, paste0(k, "_venn_response")), W2, 92)
  rg <- do.call(rbind, regs)
  detail <- cbind(rg, annot[rg$gene_id, ANN[-1]],
                  do.call(cbind, lapply(names(x$tabs), function(n) { t <- x$tabs[[n]][match(rg$gene_id, x$tabs[[n]]$gene_id), ]
                    setNames(data.frame(round(t$log2FoldChange, 3), signif(t$padj, 3)), paste0(c("log2FC_", "padj_"), sub("_withPartner_vs_alone", "", n))) })))
  detail <- detail[order(detail$direction, -detail$n_isolates, detail$region), ]
  summ <- as.data.frame.matrix(table(detail$region, detail$direction)); summ <- cbind(Region = rownames(summ), summ)
  sheets <- list(Summary = summ, Shared_all4 = detail[detail$region == "Shared by all 4", ])
  for (s in STRAINS) sheets[[paste0("Specific_", s)]] <- detail[detail$region == paste("Specific to", s), ]
  sheets$All_DEGs <- detail
  book(sheets, file.path(TAB, paste0(FOLDER[k], "_shared_specific_genes.xlsx")))
  summ
}

cat("\n== Isolate-specific and shared genes ==\n")
SPEC <- list()
for (k in c("C1", "C2a", "C3")) {
  sp <- isolate_specific(R[[k]]); SPEC[[k]] <- sp
  book(c(list(Summary = sp$summary), setNames(sp$per_isolate, paste0(STRAINS, "_specific")), list(DEGs_any_pair = sp$all)),
       file.path(TAB, paste0(FOLDER[k], "_shared_specific_genes.xlsx")))
  cat(" ", k, "isolate-specific (higher/lower):", paste(sprintf("%s %d/%d", sp$summary$Isolate, sp$summary$Specific_higher, sp$summary$Specific_lower), collapse = ", "), "\n")
}
RESP <- list(); for (k in c("C2b", "C4")) RESP[[k]] <- response_sets(R[[k]])

# up/down sets of every contrast for ORA
for (k in names(R)) {
  x <- R[[k]]; tabs <- x$tabs; labs <- x$labels
  if (!is.null(x$extra)) { tabs$Mean_of_4_vs_alone <- x$extra; labs$Mean_of_4_vs_alone <- "Ca+Sa (mean of 4) vs Ca" }
  for (n in names(tabs)) {
    t <- tabs[[n]]; U <- t$gene_id[!is.na(t$padj)]; lf <- setNames(t$log2FoldChange, t$gene_id)
    for (dir in c("Up", "Down")) add_set(paste0(k, "|", n, "|", dir), k, "contrast", labs[[n]], dir, t$gene_id[t$regulation == dir], U, x$sp, lf, n)
  }
}

# -------------------------------------------- C1 vs C2a concordance ----
cat("\n== Concordance C1 (alone) vs C2a (with C. albicans) ==\n")
cc <- list(); cnt <- list()
for (n in names(R$C1$tabs)) {
  a <- R$C1$tabs[[n]]; b <- R$C2a$tabs[[n]]; g <- intersect(a$gene_id[!is.na(a$padj)], b$gene_id[!is.na(b$padj)])
  ra <- a$regulation[match(g, a$gene_id)]; rb <- b$regulation[match(g, b$gene_id)]
  cls <- ifelse(ra != "NS" & rb != "NS", ifelse(ra == rb, "Both, same direction", "Both, opposite direction"),
         ifelse(ra != "NS", "Only alone", ifelse(rb != "NS", "Only with C. albicans", "NS")))
  d <- data.frame(Pair = R$C1$labels[[n]], annot[g, ANN], log2FC_alone = a$log2FoldChange[match(g, a$gene_id)],
                  padj_alone = a$padj[match(g, a$gene_id)], log2FC_with_Ca = b$log2FoldChange[match(g, b$gene_id)],
                  padj_with_Ca = b$padj[match(g, b$gene_id)], class = cls, check.names = FALSE)
  cc[[n]] <- d
  cnt[[n]] <- data.frame(Pair = R$C1$labels[[n]], t(as.matrix(table(factor(cls, levels = names(CONC_COL))))),
                         Spearman_r = round(cor(d$log2FC_alone, d$log2FC_with_Ca, method = "spearman", use = "complete"), 3), check.names = FALSE)
}
ccd <- do.call(rbind, cc); ccd$Pair <- factor(ccd$Pair, levels = unlist(R$C1$labels))
ccd$class <- factor(ccd$class, levels = names(CONC_COL)); ccd <- ccd[order(ccd$class), ]
rlab <- do.call(rbind, cnt); rlab$Pair <- factor(rlab$Pair, levels = levels(ccd$Pair))
pc <- ggplot(ccd, aes(log2FC_alone, log2FC_with_Ca, colour = class)) +
  geom_hline(yintercept = 0, linewidth = lw(0.3), colour = "grey80") + geom_vline(xintercept = 0, linewidth = lw(0.3), colour = "grey80") +
  geom_abline(linetype = "22", linewidth = lw(0.4), colour = "grey55") + geom_point(size = 0.6, stroke = 0, alpha = 0.85) +
  geom_text(data = rlab, aes(x = -Inf, y = Inf, label = sprintf("rho == %.2f", Spearman_r)), parse = TRUE, inherit.aes = FALSE,
            hjust = -0.15, vjust = 1.4, size = fs(SMALL), family = FONT) +
  facet_wrap(~Pair, ncol = 3) + scale_colour_manual(values = CONC_COL, name = NULL, limits = names(CONC_COL), drop = FALSE,
    labels = c("NS" = "Not DE", "Both, same direction" = "DE in both, same direction", "Both, opposite direction" = "DE in both, opposite direction",
               "Only alone" = "DE only alone", "Only with C. albicans" = "DE only with C. albicans")) +
  guides(colour = guide_legend(override.aes = list(size = 1.6))) +
  labs(x = expression(log[2]~fold~change~"(isolates alone)"), y = expression(log[2]~fold~change~"(with"~italic("C. albicans")*")")) +
  theme_pub() + theme(legend.position = "bottom", panel.border = element_rect(fill = NA, linewidth = lw(0.5)), axis.line = element_blank())
save_fig(pc, file.path(FIG, "C1_vs_C2a", "C1_vs_C2a_concordance"), W2, 120)
book(list(Counts = do.call(rbind, cnt), Genes = do.call(rbind, cc)), file.path(TAB, "C1_vs_C2a_concordance.xlsx"))

# ------------------------------------------------ C3 x C4: isolate-specific fungal responses ----
cat("\n== C3 x C4 isolate-specific fungal responses ==\n")
comb <- list()
for (s in STRAINS) {
  t4 <- R$C4$tabs[[paste0(s, "_withPartner_vs_alone")]]
  for (dir in c("Up", "Down")) {
    g <- intersect(SPEC$C3$spec[[paste(s, dir)]], t4$gene_id[t4$regulation == dir])
    lab <- paste0(s, "-specific ", if (dir == "Up") "induction" else "repression")
    if (length(g)) comb[[paste(s, dir)]] <- data.frame(Isolate = s, Class = lab, annot[g, ANN],
      log2FC_vs_Ca_alone = round(t4$log2FoldChange[match(g, t4$gene_id)], 3), padj_vs_Ca_alone = signif(t4$padj[match(g, t4$gene_id)], 3), check.names = FALSE)
    add_set(paste0("C3xC4|", s, "_", dir), "C3xC4", "C3xC4", lab, dir, g, rownames(R$C3$dds), "Calb", setNames(t4$log2FoldChange, t4$gene_id), s)
  }
}
c34 <- if (length(comb)) do.call(rbind, comb) else data.frame()
book(list(Fungal_isolate_specific = c34), file.path(TAB, "C3xC4_isolate_specific_fungal_response.xlsx"))

# ---------------------------------------------------------------- PCA ----
cat("\n== PCA ==\n")
pca_vsd <- function(sp, groups) {   # size factors + VST; no model fit needed
  smp <- ss$sample[ss$group %in% groups]; g <- annot$gene_id[annot$species == sp & use_gene]
  cts <- round(counts[g, smp, drop = FALSE]); keep <- rowSums(cts >= opt$min_count) >= opt$min_samples
  cd <- ss[smp, ]; cd$group <- factor(cd$group, levels = groups)
  dds <- estimateSizeFactors(DESeqDataSetFromMatrix(cts[keep, ], cd, ~group))
  if (nrow(dds) >= 1500) vst(dds, blind = TRUE) else varianceStabilizingTransformation(dds, blind = TRUE)
}
pca_plot <- function(v, title, sp, labelled = FALSE) {
  m <- assay(v); m <- m[order(-apply(m, 1, var))[seq_len(min(500, nrow(m)))], ]
  pc <- prcomp(t(m)); ve <- round(100 * pc$sdev^2 / sum(pc$sdev^2), 1)
  iso <- as.character(v$isolate); mono <- v$culture == "mono"
  if (sp == "Saur") {          # S. aureus: isolate names; circle = alone, triangle = with C. albicans
    lev <- STRAINS; grp <- iso; cols <- STRAIN_COL[STRAINS]
    leg <- parse(text = paste0("'", STRAINS, "'"))
    shp <- c("Alone" = 16, "With C. albicans" = 17); shp_lab <- parse(text = c("Alone", "With~italic('C. albicans')"))
    ltitle <- "Isolate"
  } else {                     # C. albicans: alone + each co-culture, coloured by the partner isolate
    lev <- c("none", STRAINS); grp <- ifelse(iso == "none", "none", iso)
    cols <- setNames(c(STRAIN_COL[["C. albicans alone"]], STRAIN_COL[STRAINS]), lev)
    leg <- parse(text = c("italic('C. albicans')~alone", paste0("italic('C. albicans')~'+ ", STRAINS, "'")))
    shp <- c("Alone" = 16, "With S. aureus" = 17); shp_lab <- parse(text = c("Alone", "With~italic('S. aureus')"))
    ltitle <- "Biofilm"
  }
  d <- data.frame(PC1 = pc$x[, 1], PC2 = pc$x[, 2], sample = colnames(v), grp = factor(grp, levels = lev),
                  shape = factor(ifelse(mono, names(shp)[1], names(shp)[2]), levels = names(shp)))
  p <- ggplot(d, aes(PC1, PC2, colour = grp, shape = shape)) + geom_point(size = 2, stroke = 0.3) +
    scale_colour_manual(values = cols, limits = lev, labels = leg, name = ltitle) +
    scale_shape_manual(values = shp, limits = names(shp), labels = shp_lab, name = NULL) +
    guides(colour = guide_legend(order = 1, ncol = if (sp == "Saur") 4 else 2, override.aes = list(shape = 16, size = 2)),
           shape = guide_legend(order = 2, override.aes = list(colour = "black", size = 2))) +
    labs(title = title, x = paste0("PC1 (", ve[1], "%)"), y = paste0("PC2 (", ve[2], "%)")) + theme_pub() +
    theme(panel.border = element_rect(fill = NA, linewidth = lw(0.5)), axis.line = element_blank(),
          legend.position = "bottom", legend.box = "vertical", legend.spacing.y = unit(0.5, "mm"))
  if (labelled) p <- p + geom_text_repel(aes(label = sample), size = fs(TINY), family = FONT, show.legend = FALSE, segment.size = lw(0.3))
  p
}
ttl_s <- expression(italic("S. aureus")~genes); ttl_c <- expression(italic("C. albicans")~genes)
vs_s <- pca_vsd("Saur", c(paste0(STRAINS, "_mono"), paste0("Calb_", STRAINS))); vs_c <- pca_vsd("Calb", c("Calb_mono", paste0("Calb_", STRAINS)))
for (lb in c(FALSE, TRUE)) {
  ps <- pca_plot(vs_s, ttl_s, "Saur", lb); pc_ <- pca_plot(vs_c, ttl_c, "Calb", lb)
  f <- (ps | pc_) + tags_pub()      # each panel keeps its own legend
  save_fig(f, file.path(FIG, "Overview", if (lb) "Fig_PCA_labelled" else "Fig_PCA"), W2, 95)
  if (!lb) { save_fig(ps, file.path(FIG, "Overview", "PCA_Saur"), W1, 90); save_fig(pc_, file.path(FIG, "Overview", "PCA_Calb"), W1, 90) }
}

# ------------------------------------------------------ expression heatmaps ----
cat("\n== Expression heatmaps ==\n")
grp_label <- function(g, sp) { g <- as.character(g)
  if (sp == "Saur") ifelse(grepl("_mono$", g), sub("_mono$", "", g), paste(sub("^Calb_", "", g), "+ Ca"))
  else ifelse(g == "Calb_mono", "Ca alone", paste("Ca +", sub("^Calb_", "", g))) }
sample_order <- function(v, sp) {
  cd <- as.data.frame(colData(v)); cd$glab <- factor(grp_label(cd$group, sp), levels = unique(grp_label(levels(cd$group), sp)))
  cd <- cd[order(cd$group, cd$replicate), ]; cd
}
ann_panel <- function(cd) {   # strain + biofilm colour strip above a heatmap
  st <- ifelse(as.character(cd$isolate) == "none", "C. albicans alone", as.character(cd$isolate))
  bf <- ifelse(cd$culture == "mono", "Single-species", "Dual-species")
  a <- rbind(data.frame(sample = cd$sample, row = "Strain", fill = unname(STRAIN_COL[st]), glab = cd$glab),
             data.frame(sample = cd$sample, row = "Biofilm", fill = unname(BIOFILM_COL[bf]), glab = cd$glab))
  a$sample <- factor(a$sample, levels = cd$sample); a$row <- factor(a$row, levels = c("Biofilm", "Strain"))
  ggplot(a, aes(sample, row, fill = fill)) + geom_tile(colour = "white", linewidth = lw(0.3)) + scale_fill_identity() +
    facet_grid(. ~ glab, scales = "free_x", space = "free_x") + scale_x_discrete(expand = c(0, 0)) + labs(x = NULL, y = NULL) + theme_pub() +
    theme(axis.text.x = element_blank(), axis.ticks = element_blank(), axis.line = element_blank(), axis.text.y = element_text(size = SMALL),
          strip.text = element_text(size = SMALL, face = "bold"), panel.spacing.x = unit(0.6, "mm"))
}
deg_heatmap <- function(v, sp, genes, title, file, k_max = 6) {
  cd <- sample_order(v, sp); m <- assay(v)[intersect(genes, rownames(v)), cd$sample, drop = FALSE]
  if (nrow(m) > 2000) m <- m[order(-apply(m, 1, var))[1:2000], ]
  if (nrow(m) < 4) return(NULL)
  z <- t(scale(t(m))); z[is.na(z)] <- 0; z <- pmax(pmin(z, 2.5), -2.5)
  hc <- hclust(dist(z), "ward.D2"); go <- hc$labels[hc$order]; k <- min(k_max, nrow(z) - 1)
  cl0 <- cutree(hc, k); cl <- setNames(match(cl0, unique(cl0[go])), names(cl0))      # clusters numbered top to bottom
  long <- data.frame(gene = factor(rep(rownames(z), ncol(z)), levels = rev(go)), sample = factor(rep(colnames(z), each = nrow(z)), levels = cd$sample),
                     z = as.vector(z)); long$glab <- cd[as.character(long$sample), "glab"]
  pH <- ggplot(long, aes(sample, gene, fill = z)) + geom_raster() +
    scale_fill_gradientn(colours = DIVERGING, limits = c(-2.5, 2.5), name = "Expression\n(row z-score)") +
    facet_grid(. ~ glab, scales = "free_x", space = "free_x") + scale_y_discrete(expand = c(0, 0)) + scale_x_discrete(expand = c(0, 0)) +
    labs(x = NULL, y = paste(nrow(z), "DEGs")) + theme_pub() +
    theme(axis.text.y = element_blank(), axis.ticks = element_blank(), axis.line = element_blank(), strip.text = element_blank(),
          axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = TINY), panel.spacing.x = unit(0.6, "mm")) +
    guides(fill = guide_colourbar(barheight = unit(18, "mm"), barwidth = unit(2.5, "mm")))
  pos <- tapply(seq_along(go), cl[go], median)
  cb <- data.frame(gene = factor(go, levels = rev(go)), cl = factor(cl[go]))
  pC <- ggplot(cb, aes(1, gene, fill = cl)) + geom_raster() + scale_fill_manual(values = rep(c("#7F7F7F", "#BFBFBF"), length.out = k), guide = "none") +
    annotate("text", x = 1, y = length(go) - pos + 1, label = names(pos), size = fs(SMALL), fontface = "bold", family = FONT) +
    scale_y_discrete(expand = c(0, 0)) + scale_x_continuous(expand = c(0, 0)) + theme_void()
  cap <- paste0("Z-scored variance-stabilised expression of all ", nrow(z), " ", SPN[sp], " genes differentially expressed in at least one contrast",
                " (columns = samples, grouped by condition;\nrows = genes, clustered (Ward) into ", k, " expression clusters, numbered at left). ",
                "Gene lists per cluster: Expression_heatmap_clusters.xlsx.")
  fig <- wrap_plots(A = ann_panel(cd), B = pH, C = pC, design = "#A\nCB", widths = c(0.03, 1), heights = c(0.07, 1)) + fig_title(title, cap)
  save_fig(fig, file, W2, 165)
  gm <- sapply(levels(cd$glab), function(g) rowMeans(z[go, cd$glab == g, drop = FALSE]))
  list(tab = data.frame(annot[go, ANN], cluster = cl[go], round(gm, 2), check.names = FALSE), cl = cl, universe = rownames(v))
}
cor_panel <- function(v, sp, title) {
  cd <- sample_order(v, sp); m <- assay(v); m <- m[order(-apply(m, 1, var))[seq_len(min(500, nrow(m)))], cd$sample]
  r <- cor(m); long <- data.frame(x = factor(rep(colnames(r), each = nrow(r)), levels = cd$sample), y = factor(rep(rownames(r), ncol(r)), levels = rev(cd$sample)),
                                  r = as.vector(r))
  long$gx <- cd[as.character(long$x), "glab"]; long$gy <- cd[as.character(long$y), "glab"]
  ggplot(long, aes(x, y, fill = r)) + geom_raster() +
    scale_fill_gradient(low = "#FFFFFF", high = "#252525", name = "Pearson r", breaks = scales::breaks_pretty(3), labels = scales::label_number(accuracy = 0.01)) +
    facet_grid(gy ~ gx, scales = "free", space = "free", switch = "y") + scale_x_discrete(expand = c(0, 0)) + scale_y_discrete(expand = c(0, 0)) +
    labs(x = NULL, y = NULL, title = title) + theme_pub() +
    theme(axis.text = element_blank(), axis.ticks = element_blank(), axis.line = element_blank(), panel.spacing = unit(0.4, "mm"),
          strip.text.x = element_text(size = TINY + 0.5, angle = 90, hjust = 0), strip.text.y.left = element_text(size = TINY + 0.5, angle = 0, hjust = 1),
          strip.placement = "outside", legend.position = "bottom") +
    guides(fill = guide_colourbar(barwidth = unit(30, "mm"), barheight = unit(2.5, "mm")))
}
degs_of <- function(ks) unique(unlist(lapply(ks, function(k) { x <- R[[k]]; tl <- c(x$tabs, if (!is.null(x$extra)) list(x$extra))
  unlist(lapply(tl, function(t) t$gene_id[t$regulation != "NS"])) })))
HM <- list(Saur = deg_heatmap(vs_s, "Saur", degs_of(c("C1", "C2a", "C2b")), "S. aureus - expression of all differentially expressed genes",
                              file.path(FIG, "Overview", "Fig_Heatmap_DEGs_Saur")),
           Calb = deg_heatmap(vs_c, "Calb", degs_of(c("C3", "C4")), "C. albicans - expression of all differentially expressed genes",
                              file.path(FIG, "Overview", "Fig_Heatmap_DEGs_Calb")))
book(lapply(Filter(Negate(is.null), HM), `[[`, "tab"), file.path(TAB, "Expression_heatmap_clusters.xlsx"))
for (sp in names(Filter(Negate(is.null), HM))) for (i in sort(unique(HM[[sp]]$cl)))
  add_set(paste0("HM_", sp, "|cluster_", i), paste0("HM_", sp), "heatmap_cluster", paste("Cluster", i), "Cluster",
          names(HM[[sp]]$cl)[HM[[sp]]$cl == i], HM[[sp]]$universe, sp, setNames(rep(0, length(HM[[sp]]$universe)), HM[[sp]]$universe), paste("Cluster", i))
fc <- (cor_panel(vs_s, "Saur", expression(italic("S. aureus")~genes)) | cor_panel(vs_c, "Calb", expression(italic("C. albicans")~genes))) +
  fig_title("Overall similarity of the transcriptomes (sample-to-sample correlation)",
            "Pearson correlation of variance-stabilised expression (500 most variable genes of each species); samples grouped by condition.") +
  plot_annotation(tag_levels = "a") & theme(plot.tag = element_text(size = 8, face = "bold", family = FONT))
save_fig(fc, file.path(FIG, "Overview", "Fig_Sample_correlation"), W2, 105)

# ------------------------------------------------------ DEG counts ----
sm <- do.call(rbind, SUMMARY)
short <- function(k, s) if (k %in% c("C2b", "C4")) sub("_withPartner_vs_alone", "", sub("Mean_of_4_vs_alone", "Mean of 4", s)) else gsub("_vs_", " vs ", s)
sm$label <- mapply(short, sm$Comparison, sm$Sheet)
long <- rbind(data.frame(sm[, c("Comparison", "label")], dir = "Up", n = sm$Up), data.frame(sm[, c("Comparison", "label")], dir = "Down", n = -sm$Down))
long$Comparison <- factor(long$Comparison, levels = names(R), labels = TITLE[names(R)])
long$label <- factor(long$label, levels = unique(sm$label))
pd <- ggplot(long, aes(label, n, fill = dir)) + geom_col(width = 0.7) + geom_hline(yintercept = 0, linewidth = lw(0.5)) +
  geom_text(aes(label = abs(n), vjust = ifelse(n >= 0, -0.3, 1.3)), size = fs(TINY), family = FONT) +
  scale_fill_manual(values = c(Up = UP, Down = DOWN), name = NULL) +
  scale_y_continuous(labels = abs, expand = expansion(mult = 0.15)) +
  facet_wrap(~Comparison, scales = "free_x", labeller = label_parsed, ncol = 3) +
  labs(x = NULL, y = "Differentially expressed genes") + theme_pub() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "bottom")
save_fig(pd, file.path(FIG, "Overview", "Fig_DEG_counts"), W2, 110)
save_fig(colour_key_fig(), file.path(FIG, "Overview", "Colour_key"), W1 + 40, 90)

# ---------------------------------------------------------- summary ----
spec_sum <- do.call(rbind, lapply(names(SPEC), function(k) cbind(Comparison = k, SPEC[[k]]$summary)))
resp_sum <- do.call(rbind, lapply(names(RESP), function(k) cbind(Comparison = k, RESP[[k]])))
readme <- data.frame(Item = c("C1", "C2a", "C2b", "C3", "C4", "Extra", "DEG", "Isolate-specific", "Shared/specific response", "C1 vs C2a", "C3 x C4"),
  Description = c("S. aureus genes, isolates grown alone: 6 pairwise contrasts. log2FC > 0 = higher in the first isolate.",
                  "S. aureus genes, isolates in co-culture with C. albicans: the same 6 pairs.",
                  "S. aureus genes: each isolate with C. albicans vs the same isolate alone (4 contrasts).",
                  "C. albicans genes, co-cultures: the same 6 pairs, comparing the partner isolates.",
                  "C. albicans genes: C. albicans with each isolate vs C. albicans alone (4 contrasts).",
                  "C. albicans genes: mean of the four co-cultures vs C. albicans alone (sheet Mean_of_4_vs_alone in C4).",
                  sprintf("Benjamini-Hochberg adjusted P < %.2g and |log2 fold change| >= %g (DESeq2 Wald test).", PADJ, LFC),
                  "C1, C2a, C3: a gene DE in all three pairwise contrasts of an isolate, in the same direction (higher or lower than every other isolate).",
                  "C2b, C4: 4-set Venn of the per-isolate DEG lists, up and down separately; shared = DE with all 4 isolates, specific = DE with one isolate only.",
                  "For each isolate pair, log2FC alone vs with C. albicans; classes show whether a difference persists in co-culture.",
                  "C. albicans genes that are isolate-specific in C3 and DE in the same direction vs C. albicans alone in C4."))
book(list(README = readme, DEG_counts = sm[, c("Comparison", "Contrast", "Genes_tested", "Up", "Down")],
          Isolate_specific_counts = spec_sum, Response_region_counts = resp_sum, C1_vs_C2a_counts = do.call(rbind, cnt)),
     file.path(TAB, "DEG_summary.xlsx"))
saveRDS(list(sets = SETS, padj = PADJ, lfc = LFC, gene_label = setNames(gene_label(annot$gene_id), annot$gene_id),
             species = setNames(annot$species, annot$gene_id), folders = FOLDER, titles = TITLE), file.path(opt$outdir, "ora_input.rds"))
writeLines(capture.output(sessionInfo()), file.path(opt$outdir, "R_sessionInfo.txt"))
cat("\nDEG counts:\n"); print(sm[, c("Comparison", "Contrast", "Up", "Down")], row.names = FALSE)
cat("\nGene sets for ORA:", length(SETS), "\nDone:", opt$outdir, "\n")
