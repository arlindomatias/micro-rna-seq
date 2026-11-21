# ============================================================
# SCRIPT 00: shared setup (sourced by every R step of the miRNA projects)
# All dataset-specific values come from the project's config.yaml (`snakemake@config`).
# Scripts run through the Snakemake `script:` directive, which creates the `snakemake` object.
# ============================================================
if (!exists("snakemake")) {
  stop("Run this script through Snakemake (`snakemake --use-conda --cores 1`).")
}

# R messages and output go to the rule log.
log_con <- file(snakemake@log[[1]], open = "wt")
sink(log_con)
sink(log_con, type = "message")

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(ggplot2)
})

cfg <- snakemake@config

for (d in c("results", "plots", "data", "cache")) dir.create(d, showWarnings = FALSE)

options(timeout = 600, lifecycle_verbosity = "quiet")
set.seed(cfg$seed %||% 123)
if (!interactive()) pdf(NULL)

ref_level <- cfg$groups$reference
trt_level <- cfg$groups$treatment
group_colors <- unlist(cfg$groups$colors)

# DE criterion: FDR < de$padj; when no miRNA reaches it, the unadjusted p-value with the stricter cutoff
# de$pvalue_fallback (default 0.01) is used instead (exploratory; report.py applies the same rule).
# de$criterion: "auto" (default, rule above), "pvalue" (always the unadjusted p < de$pvalue_fallback) or
# "padj" (always FDR).
deg_criterion <- function(tt) {
  mode <- cfg$de$criterion %||% "auto"
  stopifnot(mode %in% c("auto", "pvalue", "padj"))
  if (mode == "pvalue") return(list(col = "pvalue", cut = cfg$de$pvalue_fallback %||% 0.01))
  if (mode == "padj" || any(tt$padj < cfg$de$padj, na.rm = TRUE)) {
    list(col = "padj", cut = cfg$de$padj)
  } else {
    list(col = "pvalue", cut = cfg$de$pvalue_fallback %||% 0.01)
  }
}

# Rows of a DE table that pass the criterion and the fold-change cutoff.
de_hits <- function(tt) {
  rule <- deg_criterion(tt)
  tt[!is.na(tt[[rule$col]]) & tt[[rule$col]] < rule$cut & abs(tt$log2FoldChange) > cfg$de$lfc, ]
}

# Saves a plot without stopping the pipeline when one figure fails; a placeholder keeps the declared output.
save_plot <- function(path, expr, width = 7, height = 5, dpi = 300) {
  ok <- tryCatch({
    p <- expr
    ggsave(path, plot = p, width = width, height = height, dpi = dpi)
    TRUE
  }, error = function(e) {
    message("Plot failed (", path, "): ", conditionMessage(e))
    FALSE
  })
  if (!ok) placeholder_plot(path, "Plot not available (see the rule log)")
  invisible(ok)
}

# Same for base-graphics plots (e.g. edgeR plotBCV / plotQLDisp): `expr` draws on an open png device.
save_base_plot <- function(path, expr, width = 6.5, height = 4.5, res = 300) {
  png(path, width = width, height = height, units = "in", res = res)
  ok <- tryCatch({ expr; TRUE }, error = function(e) {
    message("Plot failed (", path, "): ", conditionMessage(e))
    FALSE
  })
  dev.off()
  if (!ok) placeholder_plot(path, "Plot not available (see the rule log)")
  invisible(ok)
}

placeholder_plot <- function(path, text) {
  p <- ggplot() + annotate("text", x = 0, y = 0, label = text, size = 5) + theme_void()
  ggsave(path, plot = p, width = 6, height = 2, dpi = 100)
}

# Writes a table with the given columns and no rows (keeps declared outputs when a step has no results).
write_empty <- function(path, columns, what = "results") {
  message("No ", what, ": writing an empty table to ", path)
  write_csv(as_tibble(setNames(replicate(length(columns), character(), simplify = FALSE), columns)), path)
}
