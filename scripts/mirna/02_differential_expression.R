# ============================================================
# SCRIPT 02: filtering, TMM, QC plots and edgeR quasi-likelihood DE for every model in config `models:`
# Each model has its own formula and tested coefficient; non-estimable design columns (e.g. an interaction
# whose main effect is absorbed by the patient term) are dropped before fitting.
# Outputs: results/de_<model>.csv, results/de_summary.tsv, results/candidates.tsv, results/log2_cpm.csv,
# data/fits.rds and the QC/DE plots.
# Positive log2FoldChange = higher in the treatment (later) timepoint, or, for an interaction slope, a larger
# increase per unit of the covariate.
# ============================================================
source(file.path(snakemake@scriptdir, "00_common.R"))
suppressPackageStartupMessages({
  library(edgeR)
  library(limma)
  library(pheatmap)
  library(ggrepel)
  library(patchwork)
})
out <- snakemake@output
models <- cfg$models
primary <- cfg$primary_model
stopifnot(primary %in% names(models))
lfc_cut <- cfg$de$lfc
candidates <- unlist(cfg$candidates)

dge_all <- readRDS(snakemake@input$dge)
coldata <- dge_all$samples
coldata$timepoint <- factor(coldata$timepoint, levels = c(ref_level, trt_level))
coldata$patient <- factor(coldata$patient)

### 1. Filter and normalize -------------------------------------------------------
min_samples <- cfg$filter$min_samples %||% min(table(coldata$timepoint))
keep <- rowSums(dge_all$counts >= cfg$filter$min_count) >= min_samples
message("Filter: >= ", cfg$filter$min_count, " reads in >= ", min_samples, " samples -> ",
        sum(keep), " of ", nrow(dge_all), " miRNAs")
dge <- dge_all[keep, , keep.lib.sizes = FALSE]
dge <- normLibSizes(dge, method = "TMM")
message("TMM factors: ", paste(sprintf("%s=%.2f", colnames(dge), dge$samples$norm.factors), collapse = ", "))

logcpm <- cpm(dge, log = TRUE, prior.count = 2)
write_csv(data.frame(mirna = rownames(logcpm), round(logcpm, 4), check.names = FALSE), out$logcpm)
cpm_all <- cpm(dge_all)          # every miRNA, library-size CPM (for candidates that are filtered out)
cpm_tmm <- cpm(dge)

sample_label <- setNames(paste0(coldata$patient, " ", coldata$timepoint), coldata$sample_id)

### 2. Models --------------------------------------------------------------------
fit_model <- function(name, m) {
  design <- model.matrix(as.formula(m$formula), data = coldata)
  drop <- nonEstimable(design)
  if (length(drop)) {
    message(name, ": dropping non-estimable columns ", paste(drop, collapse = ", "))
    design <- design[, !colnames(design) %in% drop, drop = FALSE]
  }
  if (!m$coef %in% colnames(design)) {
    stop(name, ": coefficient ", m$coef, " not in the design (", paste(colnames(design), collapse = ", "), ")")
  }
  message(name, ": ", m$formula, " | ", ncol(design), " coefficients, residual df ", nrow(design) - ncol(design))
  y <- estimateDisp(dge, design, robust = TRUE)
  fit <- glmQLFit(y, design, robust = TRUE)
  test <- glmQLFTest(fit, coef = m$coef)
  tab <- topTags(test, n = Inf, sort.by = "none")$table
  tt <- tibble(
    mirna          = rownames(tab),
    logCPM         = tab$logCPM,
    log2FoldChange = tab$logFC,
    F              = tab$F,
    stat           = sign(tab$logFC) * sqrt(tab$F),
    pvalue         = tab$PValue,
    padj           = tab$FDR,
    !!paste0("mean_cpm_", ref_level) := rowMeans(cpm_tmm[, coldata$timepoint == ref_level, drop = FALSE]),
    !!paste0("mean_cpm_", trt_level) := rowMeans(cpm_tmm[, coldata$timepoint == trt_level, drop = FALSE]),
    detected       = rowSums(dge$counts > 0)
  ) |> arrange(pvalue)
  rule <- deg_criterion(tt)
  hits <- de_hits(tt)
  if (rule$col == "pvalue") message(name, ": no miRNA with FDR < ", cfg$de$padj,
                                    " (min ", signif(min(tt$padj), 3), "); exploratory p < ", rule$cut)
  summary <- tibble(
    model = name, formula = m$formula, coef = m$coef, samples = ncol(dge), mirnas_tested = nrow(tt),
    common_bcv = round(sqrt(y$common.dispersion), 3), prior_df = round(median(fit$df.prior), 2),
    padj_cutoff = cfg$de$padj, de_criterion = rule$col, de_cutoff = rule$cut, lfc_cutoff = lfc_cut,
    de_up = sum(hits$log2FoldChange > 0), de_down = sum(hits$log2FoldChange < 0),
    padj_only = sum(tt$padj < cfg$de$padj, na.rm = TRUE), min_padj = signif(min(tt$padj), 3),
    p_below_fallback = sum(tt$pvalue < (cfg$de$pvalue_fallback %||% 0.01)),
    expected_by_chance = round(nrow(tt) * (cfg$de$pvalue_fallback %||% 0.01), 1)
  )
  print(as.data.frame(summary))
  list(design = design, disp = y, fit = fit, table = tt, summary = summary)
}

fits <- setNames(lapply(names(models), function(n) fit_model(n, models[[n]])), names(models))
saveRDS(fits, out$fit)
for (i in seq_along(fits)) {
  path <- out$tables[[i]]
  stopifnot(basename(path) == paste0("de_", names(fits)[i], ".csv"))
  write_csv(fits[[i]]$table, path)
}
de_summary <- bind_rows(lapply(fits, `[[`, "summary"))
write_tsv(de_summary, out$summary)

### 3. Candidate miRNAs from the paper --------------------------------------------
# For every model: the edgeR result when the miRNA passed the filter; plus a paired Wilcoxon test on
# log2(CPM + 1) over every miRNA (independent of the filter), as a non-parametric check.
pairs <- coldata |> dplyr::select(sample_id, patient, timepoint) |>
  pivot_wider(names_from = timepoint, values_from = sample_id) |>
  filter(!is.na(.data[[ref_level]]), !is.na(.data[[trt_level]]))
cand_tab <- bind_rows(lapply(candidates, function(mir) {
  in_data <- mir %in% rownames(dge_all)
  v <- if (in_data) log2(cpm_all[mir, ] + 1) else NULL
  d <- if (in_data) v[pairs[[trt_level]]] - v[pairs[[ref_level]]] else NA
  base <- tibble(
    mirna = mir, in_data = in_data, tested = mir %in% rownames(dge),
    reads_total = if (in_data) sum(dge_all$counts[mir, ]) else NA,
    detected = if (in_data) sum(dge_all$counts[mir, ] > 0) else NA,
    samples_ge_min_count = if (in_data) sum(dge_all$counts[mir, ] >= cfg$filter$min_count) else NA,
    patients_up = if (in_data) sum(d > 0) else NA,
    patients_down = if (in_data) sum(d < 0) else NA,
    median_paired_log2_change = if (in_data) median(d) else NA,
    paired_wilcoxon_p = if (in_data) suppressWarnings(wilcox.test(d, exact = FALSE)$p.value) else NA
  )
  per_model <- bind_cols(lapply(names(fits), function(n) {
    r <- fits[[n]]$table[fits[[n]]$table$mirna == mir, ]
    tibble(!!paste0(n, "_log2FC") := if (nrow(r)) r$log2FoldChange else NA_real_,
           !!paste0(n, "_pvalue") := if (nrow(r)) r$pvalue else NA_real_,
           !!paste0(n, "_padj")   := if (nrow(r)) r$padj else NA_real_)
  }))
  bind_cols(base, per_model)
}))
write_tsv(cand_tab, out$candidates)
print(as.data.frame(cand_tab))

### 4. QC plots ---------------------------------------------------------------------
message("Creating QC plots")
tp_colors <- group_colors
lib_df <- tibble(sample_id = colnames(dge_all), label = sample_label[colnames(dge_all)],
                 timepoint = coldata$timepoint, reads = colSums(dge_all$counts),
                 detected = colSums(dge_all$counts > 0))

# <figure name="library_sizes">
# Reads and detected miRNAs per library
p1 <- ggplot(lib_df, aes(x = reorder(label, reads), y = reads / 1e6, fill = timepoint)) +
  geom_col() + coord_flip() + scale_fill_manual(values = tp_colors) +
  labs(x = NULL, y = "Reads (millions)", fill = "Timepoint") + theme_bw(base_size = 11)
p2 <- ggplot(lib_df, aes(x = reorder(label, reads), y = detected, fill = timepoint)) +
  geom_col() + coord_flip() + scale_fill_manual(values = tp_colors) +
  labs(x = NULL, y = "miRNAs with ≥ 1 read", fill = "Timepoint") + theme_bw(base_size = 11)
save_plot(out$library, p1 + p2 + plot_layout(guides = "collect"), width = 10, height = 5)
# </figure>

# <figure name="top_mirnas">
# Share of each library taken by the 10 most abundant miRNAs (over all samples)
top10 <- names(sort(rowSums(cpm_all), decreasing = TRUE))[1:10]
share <- as.data.frame(t(dge_all$counts[top10, ]) / colSums(dge_all$counts)) |>
  rownames_to_column("sample_id") |>
  mutate(other = 1 - rowSums(across(all_of(top10))), label = sample_label[sample_id]) |>
  pivot_longer(c(all_of(top10), other), names_to = "mirna", values_to = "share") |>
  mutate(mirna = factor(sub("^hsa-", "", mirna), levels = c(sub("^hsa-", "", top10), "other")))
save_plot(out$top_mirnas,
          ggplot(share, aes(x = label, y = share, fill = mirna)) +
            geom_col() + coord_flip() +
            scale_fill_manual(values = c(scales::hue_pal()(10), "grey80")) +
            scale_y_continuous(labels = scales::percent) +
            labs(x = NULL, y = "Share of reads", fill = "miRNA") + theme_bw(base_size = 11),
          width = 8, height = 5.5)
# </figure>

# <figure name="hemolysis">
# Hemolysis ratio: log2((CPM miR-451a + 1) / (CPM miR-23a-3p + 1)), paired lines per patient
num <- cfg$qc$hemolysis_numerator
den <- cfg$qc$hemolysis_denominator
if (all(c(num, den) %in% rownames(cpm_all))) {
  hem <- tibble(sample_id = colnames(cpm_all), patient = coldata$patient, timepoint = coldata$timepoint,
                ratio = log2((cpm_all[num, ] + 1) / (cpm_all[den, ] + 1)),
                den_zero = dge_all$counts[den, ] == 0)
  save_plot(out$hemolysis,
            ggplot(hem, aes(x = timepoint, y = ratio, group = patient)) +
              geom_line(colour = "grey60") +
              geom_point(aes(colour = timepoint, shape = den_zero), size = 3) +
              geom_text_repel(data = hem[hem$timepoint == trt_level, ], aes(label = patient),
                              nudge_x = 0.15, size = 3, direction = "y") +
              scale_colour_manual(values = tp_colors) +
              scale_shape_manual(values = c(`FALSE` = 16, `TRUE` = 4),
                                 labels = c(`FALSE` = "detected", `TRUE` = "0 reads")) +
              labs(x = NULL, y = sprintf("log2 (%s + 1) / (%s + 1), CPM", sub("^hsa-", "", num),
                                         sub("^hsa-", "", den)),
                   colour = "Timepoint", shape = sub("^hsa-", "", den)) +
              theme_bw(base_size = 12),
            width = 6.5, height = 5)
} else {
  placeholder_plot(out$hemolysis, "Hemolysis miRNAs not in the count table")
}
# </figure>

# <figure name="density_logcpm">
# Density of log2 CPM (TMM) of the tested miRNAs
dens <- as.data.frame(logcpm) |> rownames_to_column("mirna") |>
  pivot_longer(-mirna, names_to = "sample_id", values_to = "logcpm") |>
  mutate(timepoint = coldata[sample_id, "timepoint"])
save_plot(out$density,
          ggplot(dens, aes(x = logcpm, group = sample_id, colour = timepoint)) +
            geom_density() + scale_colour_manual(values = tp_colors) +
            labs(x = "log2 CPM (TMM)", y = "Density", colour = "Timepoint") + theme_bw(base_size = 12),
          width = 7, height = 4.5)
# </figure>

# <figure name="mds">
# MDS of log2 CPM (limma plotMDS on the tested miRNAs), the two draws of a patient joined
mds <- plotMDS(logcpm, top = min(500, nrow(logcpm)), plot = FALSE)
mds_df <- tibble(sample_id = colnames(logcpm), x = mds$x, y = mds$y,
                 patient = coldata$patient, timepoint = coldata$timepoint)
save_plot(out$mds,
          ggplot(mds_df, aes(x = x, y = y)) +
            geom_line(aes(group = patient), colour = "grey70") +
            geom_point(aes(colour = timepoint), size = 3.5) +
            geom_text_repel(aes(label = patient), size = 3) +
            scale_colour_manual(values = tp_colors) +
            labs(x = sprintf("Leading logFC dim 1 (%.0f%%)", 100 * mds$var.explained[1]),
                 y = sprintf("Leading logFC dim 2 (%.0f%%)", 100 * mds$var.explained[2]),
                 colour = "Timepoint") +
            theme_bw(base_size = 12),
          width = 7, height = 5.5)
# </figure>

# <figure name="sample_correlation">
# Pearson correlation of log2 CPM centred per miRNA (so the shared plasma profile does not dominate)
centred <- logcpm - rowMeans(logcpm)
cc <- cor(centred)
dimnames(cc) <- list(sample_label[rownames(cc)], sample_label[colnames(cc)])
ann <- data.frame(timepoint = coldata$timepoint, row.names = sample_label[coldata$sample_id])
save_plot(out$correlation,
          pheatmap(cc, annotation_row = ann, annotation_col = ann,
                   annotation_colors = list(timepoint = tp_colors),
                   color = colorRampPalette(c("#3B4CC0", "white", "#B40426"))(101),
                   breaks = seq(-1, 1, length.out = 102), display_numbers = TRUE, number_format = "%.1f",
                   fontsize_number = 6, silent = TRUE)$gtable,
          width = 9, height = 8)
# </figure>

# <figure name="bcv">
# Biological coefficient of variation (primary model)
save_base_plot(out$bcv, {
  plotBCV(fits[[primary]]$disp, main = paste(primary, "model"))
})
# </figure>

# <figure name="ql_dispersion">
# Quasi-likelihood dispersions (primary model)
save_base_plot(out$ql, {
  plotQLDisp(fits[[primary]]$fit, main = paste(primary, "model"))
})
# </figure>

### 5. DE plots -------------------------------------------------------------------
for (i in seq_along(fits)) {
  n <- names(fits)[i]
  tt <- fits[[n]]$table
  rule <- deg_criterion(tt)
  hits <- de_hits(tt)
  p_line <- if (rule$col == "pvalue") rule$cut else max(tt$pvalue[tt$padj < rule$cut])

  # <figure name="volcano">
  # Volcano: log2FC vs -log10 p; DE miRNAs coloured, top 10 and the paper's candidates labelled
  vd <- tt |> mutate(
    class = case_when(mirna %in% hits$mirna & log2FoldChange > 0 ~ "up",
                      mirna %in% hits$mirna & log2FoldChange < 0 ~ "down", TRUE ~ "ns"),
    lab = ifelse(mirna %in% c(head(tt$mirna, 10), candidates), sub("^hsa-", "", mirna), NA))
  save_plot(out$volcano[[i]],
            ggplot(vd, aes(x = log2FoldChange, y = -log10(pvalue))) +
              geom_point(aes(colour = class), size = 1.8, alpha = 0.8) +
              geom_text_repel(aes(label = lab), size = 3, max.overlaps = 30, na.rm = TRUE) +
              geom_vline(xintercept = c(-lfc_cut, lfc_cut), linetype = 3, colour = "grey50") +
              geom_hline(yintercept = -log10(p_line), linetype = 3, colour = "grey50") +
              scale_colour_manual(values = c(up = "#CC0000", down = "#3399FF", ns = "grey70"),
                                  breaks = c("up", "down")) +
              labs(x = "log2 fold change", y = "-log10 p-value", colour = NULL,
                   subtitle = sprintf("%s (%s): %d up, %d down at %s < %g and |log2FC| > %g", n,
                                      models[[n]]$formula, sum(vd$class == "up"), sum(vd$class == "down"),
                                      ifelse(rule$col == "padj", "FDR", "p"), rule$cut, lfc_cut)) +
              theme_bw(base_size = 12),
            width = 7.5, height = 5.5)
  # </figure>

  # <figure name="pvalue_histogram">
  # Histogram of p-values
  save_plot(out$pvalues[[i]],
            ggplot(tt, aes(x = pvalue)) +
              geom_histogram(breaks = seq(0, 1, 0.05), fill = "grey60", colour = "white") +
              geom_hline(yintercept = nrow(tt) / 20, linetype = 2, colour = "#CC0000") +
              labs(x = "p-value", y = "miRNAs", subtitle = sprintf("%s; dashed = uniform expectation", n)) +
              theme_bw(base_size = 12),
            width = 6, height = 4)
  # </figure>
}

# <figure name="heatmap">
# Top miRNAs of the primary model (by p-value): log2 CPM centred within each patient
tt_p <- fits[[primary]]$table
top <- head(tt_p$mirna, cfg$plots$heatmap_max_mirnas %||% 40)
within <- logcpm[top, , drop = FALSE]
for (p in levels(coldata$patient)) {
  s <- coldata$sample_id[coldata$patient == p]
  within[, s] <- within[, s] - rowMeans(within[, s, drop = FALSE])
}
ord <- coldata$sample_id[order(coldata$timepoint, coldata$patient)]
within <- within[, ord]
colnames(within) <- sample_label[ord]
rownames(within) <- sub("^hsa-", "", rownames(within))
lim <- max(abs(within))
save_plot(out$heatmap,
          pheatmap(within, cluster_cols = FALSE, annotation_col = ann,
                   annotation_colors = list(timepoint = tp_colors),
                   color = colorRampPalette(c("#3399FF", "white", "#CC0000"))(101),
                   breaks = seq(-lim, lim, length.out = 102), fontsize_row = 7, silent = TRUE,
                   gaps_col = sum(coldata$timepoint == ref_level))$gtable,
          width = 8, height = 2 + 0.14 * length(top))
# </figure>

# <figure name="model_comparison">
# Primary model against each other model: -log10 p of the same miRNA
others <- setdiff(names(fits), primary)
cmp <- bind_rows(lapply(others, function(n) {
  inner_join(fits[[primary]]$table |> dplyr::select(mirna, p_primary = pvalue),
             fits[[n]]$table |> dplyr::select(mirna, p_other = pvalue), by = "mirna") |>
    mutate(model = n)
}))
if (nrow(cmp)) {
  save_plot(out$model_comparison,
            ggplot(cmp, aes(x = -log10(p_primary), y = -log10(p_other))) +
              geom_abline(linetype = 2, colour = "grey50") +
              geom_point(aes(colour = mirna %in% candidates), size = 1.6, alpha = 0.8) +
              geom_text_repel(data = cmp[cmp$mirna %in% candidates, ], aes(label = sub("^hsa-", "", mirna)),
                              size = 3) +
              scale_colour_manual(values = c(`FALSE` = "grey40", `TRUE` = "#CC0000"), guide = "none") +
              facet_wrap(~ model) +
              labs(x = sprintf("-log10 p, %s model", primary), y = "-log10 p, other model") +
              theme_bw(base_size = 12),
            width = 4 + 3.5 * length(others), height = 4.5)
} else {
  placeholder_plot(out$model_comparison, "Only one model")
}
# </figure>

# <figure name="candidates">
# Candidate miRNAs of the paper: log2(CPM + 1) per patient at both draws (every miRNA, before filtering)
cand_present <- intersect(candidates, rownames(cpm_all))
if (length(cand_present)) {
  cd <- as.data.frame(t(log2(cpm_all[cand_present, , drop = FALSE] + 1))) |>
    rownames_to_column("sample_id") |>
    pivot_longer(-sample_id, names_to = "mirna", values_to = "log2cpm") |>
    mutate(patient = coldata[sample_id, "patient"], timepoint = coldata[sample_id, "timepoint"],
           mirna = sub("^hsa-", "", mirna))
  save_plot(out$candidates_plot,
            ggplot(cd, aes(x = timepoint, y = log2cpm, group = patient)) +
              geom_line(colour = "grey60") +
              geom_point(aes(colour = timepoint), size = 3) +
              geom_text_repel(data = cd[cd$timepoint == trt_level, ], aes(label = patient),
                              nudge_x = 0.15, size = 3, direction = "y") +
              scale_colour_manual(values = tp_colors) +
              facet_wrap(~ mirna, scales = "free_y") +
              labs(x = NULL, y = "log2 (CPM + 1)", colour = "Timepoint") +
              theme_bw(base_size = 12),
            width = 4 + 3 * length(cand_present), height = 4.5)
} else {
  placeholder_plot(out$candidates_plot, "No candidate miRNA in the count table")
}
# </figure>

message("Done")
