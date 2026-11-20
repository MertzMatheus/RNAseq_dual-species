#!/usr/bin/env Rscript
# Over-representation analysis (GO BP, GO MF, KEGG) of the gene sets produced by differential_expression.R:
# up/down genes of every contrast, isolate-specific genes (C1, C2a, C3), shared and isolate-specific
# responses (C2b, C4), C3 x C4 fungal responses and the expression-heatmap clusters.
#
# One-sided hypergeometric test (equivalent to clusterProfiler::enricher) with Benjamini-Hochberg correction.
# Background: genes tested in the corresponding contrast with >= 1 annotation in the database.
# Figures: per contrast a dot plot + gene-term network; per comparison one summary dot matrix.
suppressPackageStartupMessages({ library(optparse); library(ggplot2); library(openxlsx); library(patchwork); library(ggraph); library(igraph) })
.script_dir <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])))
source(file.path(.script_dir, "fig_style.R"))
opt <- parse_args(OptionParser(option_list = list(
  make_option("--ora_input"), make_option("--term_dir"), make_option("--outdir"),
  make_option("--padj", type = "double", default = 0.05), make_option("--min_size", type = "integer", default = 5),
  make_option("--max_size", type = "integer", default = 500), make_option("--min_overlap", type = "integer", default = 2),
  make_option("--redundancy", type = "double", default = 0.7), make_option("--top", type = "integer", default = 5))))
set.seed(42)
IN <- readRDS(opt$ora_input); SETS <- IN$sets; GL <- IN$gene_label
TAB <- file.path(opt$outdir, "tables"); FIG <- file.path(opt$outdir, "figures")
dir.create(TAB, recursive = TRUE, showWarnings = FALSE); dir.create(FIG, recursive = TRUE, showWarnings = FALSE)
DBS <- c(GO_BP = "GO BP", GO_MF = "GO MF", KEGG = "KEGG")

# ------------------------------------------------------------ annotations ----
T2G <- list()
for (sp in c("Calb", "Saur")) {
  f <- file.path(opt$term_dir, paste0(sp, "_term2gene.tsv.gz"))
  t <- read.delim(gzfile(f), quote = "", colClasses = "character"); t <- t[t$database %in% names(DBS), ]
  T2G[[sp]] <- lapply(split(t, t$database), function(d) list(sets = split(d$gene_id, d$term_id),
                                                             names = tapply(d$term_name, d$term_id, `[`, 1)))
}

# ------------------------------------------------------------------- ORA ----
ora <- function(genes, universe, sp, db) {
  a <- T2G[[sp]][[db]]; if (is.null(a)) return(NULL)
  sets <- lapply(a$sets, function(s) intersect(s, universe))
  U <- unique(unlist(sets)); L <- intersect(genes, U)
  if (length(L) < opt$min_overlap) return(NULL)
  sz <- lengths(sets); sets <- sets[sz >= opt$min_size & sz <= opt$max_size]
  k <- vapply(sets, function(s) sum(s %in% L), 0L); sets <- sets[k > 0]; k <- k[k > 0]
  if (!length(sets)) return(NULL)
  M <- lengths(sets); N <- length(U); n <- length(L)
  p <- phyper(k - 1, M, N - M, n, lower.tail = FALSE)
  hit <- lapply(sets, function(s) intersect(s, L))
  r <- data.frame(database = unname(DBS[db]), term_id = names(sets), term_name = unname(a$names[names(sets)]), count = k,
                  list_size = n, term_size = M, background = N, GeneRatio = paste0(k, "/", n), BgRatio = paste0(M, "/", N),
                  gene_ratio = k / n, fold_enrichment = round((k / n) / (M / N), 3), pvalue = p, padj = p.adjust(p, "BH"),
                  genes = vapply(hit, paste, "", collapse = ";"),
                  gene_names = vapply(hit, function(h) paste(GL[h], collapse = ";"), ""), stringsAsFactors = FALSE)
  r <- r[order(r$padj, r$pvalue), ]
  r$significant <- r$padj < opt$padj & r$count >= opt$min_overlap
  r$representative <- FALSE; kept <- list()
  for (i in which(r$significant)) {
    g <- hit[[r$term_id[i]]]
    if (!any(vapply(kept, function(k2) length(intersect(g, k2)) / length(union(g, k2)) >= opt$redundancy, TRUE))) {
      r$representative[i] <- TRUE; kept[[length(kept) + 1]] <- g }
  }
  r
}
cat("Running ORA on", length(SETS), "gene sets\n")
RES <- list()
for (id in names(SETS)) {
  s <- SETS[[id]]
  rr <- do.call(rbind, lapply(names(DBS), function(db) ora(s$genes, s$universe, s$species, db)))
  if (!is.null(rr) && nrow(rr)) RES[[id]] <- cbind(set_id = id, comparison = s$comparison, kind = s$kind, set = s$label,
                                                  direction = s$direction, group = s$group, rr, stringsAsFactors = FALSE)
}
ALL <- do.call(rbind, RES)
SIG <- if (is.null(ALL)) data.frame() else ALL[ALL$significant, ]

# --------------------------------------------------------------- tables ----
book <- function(sheets, file) {
  wb <- createWorkbook(); hs <- createStyle(textDecoration = "bold", fgFill = "#E7E6E6", wrapText = TRUE)
  for (n in names(sheets)) { sn <- substr(n, 1, 31); addWorksheet(wb, sn); writeData(wb, sn, sheets[[n]], headerStyle = hs)
    freezePane(wb, sn, firstRow = TRUE); setColWidths(wb, sn, cols = seq_len(max(1, ncol(sheets[[n]]))), widths = 14) }
  saveWorkbook(wb, file, overwrite = TRUE)
}
summ <- do.call(rbind, lapply(names(SETS), function(id) { s <- SETS[[id]]
  data.frame(comparison = s$comparison, kind = s$kind, set = s$label, direction = s$direction, set_id = id, genes = length(s$genes),
             GO_BP = sum(SIG$set_id == id & SIG$database == "GO BP"), GO_MF = sum(SIG$set_id == id & SIG$database == "GO MF"),
             KEGG = sum(SIG$set_id == id & SIG$database == "KEGG")) }))
cols_out <- c("set", "direction", "database", "term_id", "term_name", "count", "GeneRatio", "BgRatio", "fold_enrichment",
              "pvalue", "padj", "representative", "gene_names", "genes")
for (k in unique(summ$comparison)) {
  sk <- SIG[SIG$comparison == k, ]
  sheets <- list(Summary = summ[summ$comparison == k, ])
  for (kd in unique(summ$kind[summ$comparison == k])) {
    d <- sk[sk$kind == kd, cols_out]; sheets[[switch(kd, contrast = "Contrasts", isolate_specific = "Isolate_specific",
      shared_response = "Shared_response", specific_response = "Isolate_specific_response", C3xC4 = "C3xC4",
      heatmap_cluster = "Expression_clusters", kd)]] <- d
  }
  fn <- if (k %in% names(IN$folders)) IN$folders[[k]] else k
  book(sheets, file.path(TAB, paste0(fn, "_ORA.xlsx")))
}
if (!is.null(ALL)) write.csv(ALL, file.path(TAB, "ORA_all_tested_terms.csv"), row.names = FALSE)
readme <- data.frame(Item = c("Test", "Background", "Significance", "Redundancy", "Databases", "Sets"),
  Description = c("One-sided hypergeometric test (over-representation), Benjamini-Hochberg correction per gene set and database.",
                  "Genes tested in the contrast (or comparison, for derived sets) with at least one annotation in the database.",
                  sprintf("Adjusted P < %.2g and at least %d genes of the set in the term; terms with %d-%d background genes.", opt$padj, opt$min_overlap, opt$min_size, opt$max_size),
                  sprintf("representative = FALSE when the term's genes overlap >= %.0f%% (Jaccard) with a more significant term; such terms are kept in tables, hidden in figures.", 100 * opt$redundancy),
                  "GO Biological Process and Molecular Function (propagated to ancestor terms), KEGG pathways.",
                  "Up/Down DEGs of every contrast; isolate-specific genes; shared and isolate-specific responses; C3 x C4 isolate-specific fungal responses; expression-heatmap clusters."))
book(list(README = readme, Summary = summ), file.path(TAB, "ORA_summary.xlsx"))

# --------------------------------------------------------------- figures ----
# One visual language for every enrichment figure: colour = direction (red up, blue down; grey = expression cluster),
# point size = genes in the term, shading = significance (-log10 adjusted P).
DIR_COL <- c(Up = UP, Down = DOWN, Cluster = BAR_COL)
sig_alpha <- function() scale_alpha_continuous(range = c(0.3, 1), limits = c(-log10(opt$padj), 6), oob = scales::squish,
                                               name = expression(-log[10]~italic(P)[adj]))
term_lab <- function(x) cut_lab(x, 50)
pick_terms <- function(r, top = opt$top) {    # top non-redundant terms per database and direction
  r <- r[r$significant & r$representative, ]
  if (!nrow(r)) return(r)
  do.call(rbind, lapply(split(r, list(r$database, r$direction), drop = TRUE), function(d) head(d[order(d$padj), ], top)))
}
order_terms <- function(d) factor(term_lab(d$term_name), levels = rev(unique(term_lab(d$term_name[order(factor(d$database, levels = DBS), d$padj)]))))
dot_panel <- function(d) {      # one contrast: gene ratio per term
  d$database <- factor(d$database, levels = DBS); d$y <- order_terms(d); d$direction <- factor(d$direction, levels = c("Up", "Down"))
  ggplot(d, aes(gene_ratio, y)) +
    geom_point(aes(size = count, colour = direction, alpha = -log10(padj)), show.legend = TRUE) +
    scale_colour_manual(values = DIR_COL, name = NULL, drop = FALSE) + sig_alpha() + scale_size_area(max_size = 3.5, name = "Genes") +
    scale_x_continuous(limits = c(0, NA), expand = expansion(mult = c(0, 0.12))) +
    facet_grid(database ~ ., scales = "free_y", space = "free_y") + labs(x = "Gene ratio (DEGs in term / DEGs)", y = NULL) + theme_pub() +
    theme(strip.text.y = element_text(angle = 0, hjust = 0), axis.text.y = element_text(size = SMALL),
          panel.grid.major.y = element_line(linewidth = lw(0.25), colour = "grey92"), legend.position = "bottom", legend.box = "vertical") +
    guides(colour = guide_legend(order = 1, override.aes = list(size = 2.5, alpha = 1)), size = guide_legend(order = 2),
           alpha = guide_legend(order = 3, override.aes = list(size = 2.5)))
}
net_panel <- function(d, lfc, n_terms = 5, max_genes = 50) {
  d <- head(d[order(d$padj), ], n_terms)
  ed <- do.call(rbind, lapply(seq_len(nrow(d)), function(i) data.frame(from = d$term_id[i], to = strsplit(d$genes[i], ";")[[1]])))
  gn <- unique(ed$to); if (length(gn) > max_genes) { gn <- names(sort(abs(lfc[gn]), decreasing = TRUE))[1:max_genes]; ed <- ed[ed$to %in% gn, ] }
  nodes <- rbind(data.frame(name = d$term_id, type = "term", label = wrap_lab(cut_lab(d$term_name, 60), 22), db = as.character(d$database), lfc = NA, size = d$count),
                 data.frame(name = gn, type = "gene", label = GL[gn], db = NA, lfc = unname(lfc[gn]), size = 1))
  g <- graph_from_data_frame(ed, vertices = nodes, directed = FALSE)
  m <- max(abs(nodes$lfc), 1, na.rm = TRUE); set.seed(7)
  ggraph(g, layout = "fr") + geom_edge_link(colour = "grey75", edge_width = lw(0.4)) +
    geom_node_point(aes(filter = type == "gene", fill = lfc), shape = 21, size = 1.8, stroke = 0.2, colour = "grey30") +
    geom_node_point(aes(filter = type == "term", colour = db, size = size), shape = 15) +
    geom_node_text(aes(filter = type == "gene", label = label), size = fs(TINY), family = FONT, fontface = "italic", repel = TRUE,
                   segment.size = lw(0.2), max.overlaps = 30) +
    geom_node_label(aes(filter = type == "term", label = label), size = fs(SMALL), family = FONT, fontface = "bold", repel = TRUE,
                    label.size = lw(0.3), label.padding = unit(0.6, "mm"), fill = "white", alpha = 0.9) +
    scale_fill_gradientn(colours = DIVERGING, limits = c(-m, m), name = expression(log[2]*FC)) +
    scale_colour_manual(values = DB_COL, name = "Term") + scale_size_area(max_size = 5, name = "Genes in term") +
    theme_void(base_family = FONT, base_size = BASE) +
    theme(legend.title = element_text(size = BASE), legend.text = element_text(size = BASE - 1), legend.key.size = unit(3, "mm"),
          legend.position = "bottom", legend.box = "vertical", plot.margin = margin(2, 2, 2, 2, "mm")) +
    guides(colour = guide_legend(order = 1, override.aes = list(size = 3)), size = guide_legend(order = 2),
           fill = guide_colourbar(order = 3, barwidth = unit(18, "mm"), barheight = unit(2.5, "mm")))
}
composite <- function(k, n, label, file) {   # a dot plot  b gene-term network
  ids <- paste0(k, "|", n, "|", c("Up", "Down")); r <- do.call(rbind, RES[intersect(ids, names(RES))])
  if (is.null(r)) return(FALSE)
  d <- pick_terms(r); if (!nrow(d)) return(FALSE)
  h <- max(85, 45 + 3.6 * nrow(d))
  fig <- (dot_panel(d) | net_panel(d, SETS[[ids[1]]]$lfc)) + plot_layout(widths = c(1, 1.15)) +
    fig_title(paste(k, "-", label), "a, top 5 enriched terms per database and direction (red = among up-regulated, blue = among down-regulated genes).  b, genes of the 5 most significant terms, coloured by log2 fold change.") +
    plot_annotation(tag_levels = "a") & theme(plot.tag = element_text(size = 8, face = "bold", family = FONT))
  save_fig(fig, file, W2, h); TRUE
}
summary_panel <- function(r, xvar, xlev, title, size_lim = NULL, dirs = c("Up", "Down")) {   # terms x gene sets; up and down side by side in each column
  if (is.null(r) || !nrow(r)) return(NULL)
  d <- pick_terms(r); if (!nrow(d)) return(NULL)
  keep <- unique(d$term_id); d <- r[r$significant & r$term_id %in% keep, ]
  d$database <- factor(d$database, levels = DBS); d$y <- order_terms(d)
  d$xi <- match(d[[xvar]], xlev) + c(Up = -0.17, Down = 0.17, Cluster = 0)[d$direction]
  ggplot(d, aes(xi, y)) + geom_point(aes(size = count, colour = direction, alpha = -log10(padj)), show.legend = TRUE) +
    scale_colour_manual(values = DIR_COL, name = NULL, limits = dirs) + sig_alpha() +
    scale_size_area(max_size = 3.5, name = "Genes", limits = size_lim, breaks = scales::breaks_pretty(4)) +
    scale_x_continuous(breaks = seq_along(xlev), labels = xlev, limits = c(0.5, length(xlev) + 0.5), expand = c(0, 0)) +
    facet_grid(database ~ ., scales = "free_y", space = "free_y") + labs(x = NULL, y = NULL, title = title) + theme_pub() +
    theme(axis.text.x = element_text(angle = 40, hjust = 1), axis.text.y = element_text(size = SMALL), strip.text.y = element_text(angle = 0, hjust = 0),
          panel.grid.major = element_line(linewidth = lw(0.25), colour = "grey92"), panel.border = element_rect(fill = NA, linewidth = lw(0.3), colour = "grey60"),
          axis.line = element_blank()) +
    guides(colour = guide_legend(order = 1, override.aes = list(size = 2.5, alpha = 1)), size = guide_legend(order = 2), alpha = guide_legend(order = 3, override.aes = list(size = 2.5)))
}
n_rows <- function(p) if (is.null(p)) 0 else length(unique(p$data$y))

cat("Drawing figures\n")
made <- 0
for (k in c("C1", "C2a", "C2b", "C3", "C4")) {
  cs <- Filter(function(s) s$comparison == k && s$kind == "contrast" && s$direction == "Up", SETS)
  for (s in cs) made <- made + composite(k, s$group, s$label, file.path(FIG, k, paste0(k, "_", s$group, "_ORA")))
  is_k <- grepl(paste0("^", k, "\\|"), names(RES))
  rc <- do.call(rbind, RES[is_k & vapply(names(RES), function(i) SETS[[i]]$kind == "contrast", TRUE)])
  rd <- do.call(rbind, RES[is_k & vapply(names(RES), function(i) SETS[[i]]$kind != "contrast", TRUE)])
  if (!is.null(rd)) rd$grp <- ifelse(is.na(rd$group), "Shared by all 4", rd$group)
  xl <- if (k %in% c("C2b", "C4")) c("Shared by all 4", STRAINS) else STRAINS
  sl <- c(1, max(c(2, rc$count[rc$significant], rd$count[rd$significant]), na.rm = TRUE))
  pa <- summary_panel(rc, "set", vapply(cs, `[[`, "", "label"), "Enriched terms per contrast", sl)
  pb <- summary_panel(rd, "grp", xl, if (k %in% c("C2b", "C4")) "Enriched terms of the shared and isolate-specific responses" else "Enriched terms of the isolate-specific genes", sl)
  ps <- Filter(Negate(is.null), list(pa, pb))
  if (length(ps)) {
    hh <- vapply(ps, n_rows, 0)
    fig <- wrap_plots(ps, ncol = 1) + plot_layout(heights = pmax(hh, 2), guides = "collect") +
      fig_title(QUESTION[k], "Colour = direction (red: terms enriched among up-regulated or higher genes; blue: among down-regulated or lower genes); size = genes in the term; shading = significance.") +
      plot_annotation(tag_levels = "a") & theme(plot.tag = element_text(size = 8, face = "bold", family = FONT))
    save_fig(fig, file.path(FIG, k, paste0(k, "_ORA_summary")), W2, 40 + 3.8 * sum(hh) + 25 * length(ps))
  }
}
for (hm in c("HM_Saur", "HM_Calb")) {
  rh <- do.call(rbind, RES[grepl(paste0("^", hm, "\\|"), names(RES))])
  if (is.null(rh)) next
  xl <- unique(vapply(Filter(function(s) s$comparison == hm, SETS), `[[`, "", "label"))
  p <- summary_panel(rh, "set", xl, NULL, dirs = "Cluster")
  if (!is.null(p)) save_fig(p + fig_title(paste(if (hm == "HM_Saur") "S. aureus" else "C. albicans", "- enriched terms of the expression clusters (Fig_Heatmap_DEGs)")),
                            file.path(FIG, "Overview", paste0("Fig_Heatmap_clusters_ORA_", sub("HM_", "", hm))), W2, 40 + 3.8 * n_rows(p))
}
r34 <- do.call(rbind, RES[grepl("^C3xC4\\|", names(RES))])
if (!is.null(r34)) { p <- summary_panel(r34, "group", STRAINS, NULL)
  if (!is.null(p)) save_fig(p + fig_title("C3 x C4 - enriched terms of the isolate-specific fungal responses"), file.path(FIG, "C3xC4", "C3xC4_ORA_summary"), W2, 40 + 3.8 * n_rows(p)) }
cat("Composite figures:", made, "\n")
writeLines(capture.output(sessionInfo()), file.path(opt$outdir, "R_sessionInfo.txt"))
cat("Significant terms:", nrow(SIG), "\nDone:", opt$outdir, "\n")
