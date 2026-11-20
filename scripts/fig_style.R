# Shared figure style, sourced by differential_expression.R and enrichment.R.
#
# Journal layout: Arial (or the metric-compatible Liberation Sans), 7 pt text, 8 pt panel letters,
# 0.5 pt lines, 89 mm (single column) or 183 mm (double column) wide. Each figure is written as
# PDF (vector, embedded fonts), TIFF (600 dpi, LZW) and PNG (300 dpi). All colours are defined here.
suppressPackageStartupMessages({
  library(ggplot2); library(patchwork); library(grid); library(VennDiagram); library(systemfonts); library(ragg)
})
futile.logger::flog.threshold(futile.logger::ERROR, name = "VennDiagramLogger")

# ---------------------------------------------------------------- font ----
setup_font <- function() {
  have <- function(f) f %in% systemfonts::system_fonts()$family
  for (f in c("Arial", "Liberation Sans")) if (have(f)) return(f)
  # no Arial: register the Liberation Sans files shipped with the fontLiberation package
  src <- list.files(system.file(package = "fontLiberation"), "^LiberationSans-(Regular|Bold|Italic|BoldItalic)\\.ttf$",
                    recursive = TRUE, full.names = TRUE)
  if (length(src)) {
    dst <- file.path(Sys.getenv("HOME"), ".local", "share", "fonts"); dir.create(dst, recursive = TRUE, showWarnings = FALSE)
    file.copy(src, dst, overwrite = FALSE)
    system("fc-cache -f >/dev/null 2>&1"); systemfonts::reset_font_cache()
    if (have("Liberation Sans")) return("Liberation Sans")
  }
  warning("Neither Arial nor Liberation Sans is available: figures use the default sans-serif font")
  "sans"
}
FONT <- setup_font()
message("Figure font: ", FONT)

# ---------------------------------------------------------------- sizes ----
BASE <- 7                      # pt, all text
SMALL <- 6                     # pt, dense labels (axis text of large panels)
TINY <- 5                      # pt, gene labels in networks / heatmaps (journal minimum)
W1 <- 89; W2 <- 183            # mm, single / double column
MAX_H <- 240                   # mm, full page height
lw <- function(pt) pt / ggplot2::.pt          # line width given in points
fs <- function(pt) pt / ggplot2::.pt          # geom_text size given in points

# ---------------------------------------------------------------- colours ----
STRAINS <- c("JE2", "CI19", "S5", "CI5")
STRAIN_COL <- c("C. albicans alone" = "#4D4D4D", JE2 = "#E69F00", CI19 = "#009E73", S5 = "#0c05e6", CI5 = "#ff7003")
BIOFILM_SHAPE <- c("Single-species" = 16, "Dual-species" = 17)
BIOFILM_COL <- c("Single-species" = "#D9D9D9", "Dual-species" = "#525252")
UP <- "#B2182B"; DOWN <- "#2166AC"; NS_COL <- "#CCCCCC"
REG_COL <- c(Up = UP, Down = DOWN, NS = NS_COL)
DIVERGING <- c(DOWN, "#F7F7F7", UP)
DB_COL <- c("GO BP" = "#525252", "GO MF" = "#969696", "KEGG" = "#D9D9D9")       # term nodes (neutral, genes carry red/blue)
CONC_COL <- c("NS" = NS_COL, "Both, same direction" = "#000000", "Both, opposite direction" = "#7F7F7F",
              "Only alone" = "#8C510A", "Only with C. albicans" = "#01665E")
BAR_COL <- "#4D4D4D"

# ---------------------------------------------------------------- theme ----
theme_pub <- function(base = BASE) {
  theme_classic(base_size = base, base_family = FONT, base_line_size = lw(0.5), base_rect_size = lw(0.5)) +
    theme(text = element_text(colour = "black"),
          axis.text = element_text(size = base - 1, colour = "black"),
          axis.title = element_text(size = base),
          axis.line = element_line(linewidth = lw(0.5), colour = "black"),
          axis.ticks = element_line(linewidth = lw(0.5), colour = "black"),
          axis.ticks.length = unit(1, "mm"),
          legend.title = element_text(size = base), legend.text = element_text(size = base - 1),
          legend.key.size = unit(3, "mm"), legend.spacing = unit(1, "mm"), legend.margin = margin(0, 0, 0, 0),
          legend.background = element_blank(),
          strip.background = element_blank(), strip.text = element_text(size = base, face = "bold"),
          plot.title = element_text(size = base, face = "bold", hjust = 0, margin = margin(b = 2)),
          plot.subtitle = element_text(size = base - 1, hjust = 0),
          plot.tag = element_text(size = 8, face = "bold", family = FONT),
          plot.margin = margin(2, 2, 2, 2, "mm"))
}
tags_pub <- function() plot_annotation(tag_levels = "a") & theme(plot.tag = element_text(size = 8, face = "bold", family = FONT))

# ---------------------------------------------------------------- saving ----
save_fig <- function(p, base, width = W2, height = 120) {
  dir.create(dirname(base), recursive = TRUE, showWarnings = FALSE)
  height <- min(height, MAX_H)
  ggsave(paste0(base, ".pdf"), p, width = width, height = height, units = "mm", device = cairo_pdf, family = FONT)
  ggsave(paste0(base, ".tiff"), p, width = width, height = height, units = "mm", dpi = 600,
         device = ragg::agg_tiff, compression = "lzw", bg = "white")
  ggsave(paste0(base, ".png"), p, width = width, height = height, units = "mm", dpi = 300, device = ragg::agg_png, bg = "white")
  invisible(base)
}

# ---------------------------------------------------------------- helpers ----
strain_of <- function(isolate) ifelse(as.character(isolate) == "none", "C. albicans alone", as.character(isolate))
wrap_lab <- function(x, n = 40) vapply(x, function(s) paste(strwrap(s, n), collapse = "\n"), "")
cut_lab <- function(x, n = 55) ifelse(nchar(x) > n, paste0(substr(x, 1, n - 1), "…"), x)

# Venn diagram (2-4 sets) as a grob, in the publication font; cols named by set
venn_grob <- function(sets, cols, title = NULL, title_col = "black", subtitle = NULL) {
  sets <- sets[names(sets)]
  if (length(sets) < 2) return(nullGrob())
  cc <- unname(cols[names(sets)])
  pos <- switch(as.character(length(sets)), "2" = c(-20, 20), "3" = c(-25, 25, 180), "4" = c(-15, 15, 0, 0))
  dst <- switch(as.character(length(sets)), "2" = c(0.05, 0.05), "3" = c(0.07, 0.07, 0.06), "4" = c(0.22, 0.22, 0.11, 0.11))
  g <- venn.diagram(sets, filename = NULL, fill = cc, col = cc, alpha = 0.3, lwd = 1, cex = 7 / 12,
                    fontfamily = FONT, cat.fontfamily = FONT, cat.cex = 6.5 / 12, cat.col = cc, cat.fontface = "bold",
                    cat.pos = pos, cat.dist = dst, margin = 0.08, disable.logging = TRUE,
                    euler.d = FALSE, scaled = FALSE)   # always the classic layout, never proportional circles
  g <- gTree(children = gList(g))
  if (!is.null(title)) g <- title_grob(g, title, title_col, subtitle)
  g
}
title_grob <- function(g, title, col, subtitle = NULL) {
  lay <- grid.layout(3, 1, heights = unit(c(4.5, if (is.null(subtitle)) 0 else 4, 1), c("mm", "mm", "null")))
  fr <- frameGrob(layout = lay)
  fr <- placeGrob(fr, textGrob(title, gp = gpar(fontfamily = FONT, fontsize = BASE, fontface = "bold", col = col)), row = 1)
  if (!is.null(subtitle))
    fr <- placeGrob(fr, textGrob(subtitle, gp = gpar(fontfamily = FONT, fontsize = SMALL, col = "grey25")), row = 2)
  placeGrob(fr, g, row = 3)
}
# biological question addressed by each comparison (figure titles)
QUESTION <- c(C1 = "C1  How do the isolates differ from each other, grown alone?",
              C2a = "C2a  How do the isolates differ from each other in co-culture?",
              C2b = "C2b  Which bacterial genes respond to C. albicans?",
              C3 = "C3  Does C. albicans behave differently depending on the partner isolate?",
              C4 = "C4  How does C. albicans respond to each isolate?")
fig_title <- function(title, caption = NULL) plot_annotation(title = title, caption = caption,
  theme = theme(plot.title = element_text(size = 8, face = "bold", family = FONT),
                plot.caption = element_text(size = SMALL, family = FONT, hjust = 0, colour = "grey25", lineheight = 1.1)))

# UpSet plot built with ggplot (same fonts/sizes as every other panel); cols named by set (optional)
upset_pub <- function(sets, cols = NULL, max_bars = 20, ylab = "Genes", min_size = 1) {
  sets <- sets[lengths(sets) > 0]
  if (length(sets) < 2) return(NULL)
  if (is.null(cols)) cols <- setNames(rep(BAR_COL, length(sets)), names(sets))
  all <- unique(unlist(sets))
  M <- vapply(sets, function(s) all %in% s, logical(length(all)))
  if (is.null(dim(M))) M <- matrix(M, nrow = 1, dimnames = list(NULL, names(sets)))
  key <- apply(M, 1, function(r) paste(as.integer(r), collapse = ""))
  tab <- sort(table(key), decreasing = TRUE); tab <- tab[tab >= min_size]; tab <- head(tab, max_bars)
  if (!length(tab)) return(NULL)
  inter <- data.frame(key = names(tab), n = as.integer(tab), x = seq_along(tab))
  lev <- names(sets)[order(lengths(sets))]
  mat <- do.call(rbind, lapply(seq_len(nrow(inter)), function(i)
    data.frame(x = i, set = factor(names(sets), levels = lev), inn = strsplit(inter$key[i], "")[[1]] == "1")))
  seg <- do.call(rbind, lapply(split(mat[mat$inn, ], mat$x[mat$inn]), function(d)
    data.frame(x = d$x[1], lo = min(as.integer(d$set)), hi = max(as.integer(d$set)))))
  sz <- data.frame(set = factor(names(sets), levels = lev), n = lengths(sets))
  p_bar <- ggplot(inter, aes(x, n)) + geom_col(fill = BAR_COL, width = 0.7) +
    geom_text(aes(label = n), vjust = -0.4, size = fs(SMALL), family = FONT) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.15))) + scale_x_continuous(expand = expansion(add = 0.6)) +
    labs(x = NULL, y = paste(ylab, "in intersection")) + theme_pub() +
    theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.line.x = element_blank())
  mat$y <- as.integer(mat$set); nl <- length(lev)
  labs_df <- data.frame(y = seq_len(nl), lab = lev, col = unname(cols[lev]))
  p_mat <- ggplot(mat, aes(x, y)) +
    geom_segment(data = seg, aes(x = x, xend = x, y = lo, yend = hi), inherit.aes = FALSE, linewidth = lw(1), colour = BAR_COL) +
    geom_point(aes(colour = inn), size = 1.6) + scale_colour_manual(values = c(`TRUE` = BAR_COL, `FALSE` = "#E0E0E0"), guide = "none") +
    scale_x_continuous(expand = expansion(add = c(0.6, 0.6))) + scale_y_continuous(limits = c(0.5, nl + 0.5), expand = c(0, 0)) +
    coord_cartesian(clip = "off") + labs(x = NULL, y = NULL) + theme_pub() +
    theme(axis.text = element_blank(), axis.ticks = element_blank(), axis.line = element_blank(),
          plot.margin = margin(2, 2, 2, 2, "mm"))
  sz$y <- as.integer(sz$set)
  p_size <- ggplot(sz, aes(n, y, fill = set)) + geom_col(width = 0.6, orientation = "y") +
    scale_fill_manual(values = cols, guide = "none") +
    scale_y_continuous(limits = c(0.5, nl + 0.5), expand = c(0, 0), breaks = seq_len(nl), labels = lev, position = "right") +
    scale_x_reverse(expand = expansion(mult = c(0.15, 0))) + labs(x = paste(ylab, "per set"), y = NULL) + theme_pub() +
    theme(axis.text.y.right = element_text(size = SMALL, face = "bold", hjust = 0.5), axis.ticks.y = element_blank(), axis.line.y = element_blank())
  (plot_spacer() | p_bar) / (p_size | p_mat) + plot_layout(widths = c(1.4, 4), heights = c(2.2, 1))
}

# one-page colour legend
colour_key_fig <- function() {
  blk <- function(v, panel) data.frame(panel = panel, label = names(v), col = unname(v), i = seq_along(v))
  d <- rbind(blk(STRAIN_COL, "Strain"), blk(REG_COL, "Direction of change"), blk(DB_COL, "Database (network nodes)"),
             blk(CONC_COL, "Concordance classes"))
  d$panel <- factor(d$panel, levels = unique(d$panel))
  ggplot(d, aes(0, -i)) + geom_tile(aes(fill = col), width = 0.5, height = 0.8) + scale_fill_identity() +
    geom_text(aes(x = 0.4, label = paste0(label, "   ", col)), hjust = 0, size = fs(BASE), family = FONT) +
    facet_grid(panel ~ ., scales = "free_y", space = "free_y", switch = "y") + xlim(-0.3, 4) + theme_void(base_family = FONT) +
    theme(strip.text.y.left = element_text(size = BASE, face = "bold", angle = 0, hjust = 1), strip.placement = "outside") +
    labs(title = "Colour key",
         subtitle = "PCA: circle = grown alone, triangle = with the partner species. Heatmaps and enrichment: blue = down, red = up.") +
    theme(plot.title = element_text(size = 8, face = "bold"), plot.subtitle = element_text(size = BASE))
}
