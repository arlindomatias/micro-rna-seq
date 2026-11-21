# ============================================================
# SCRIPT 04: GO BP / KEGG over-representation of the validated targets of the DE miRNAs
# DE miRNAs of the primary model (deg_criterion: FDR, or exploratory p when nothing passes FDR), split by
# direction. Their validated targets (cache/multimir, 03_targets.R) are tested with clusterProfiler
# (enrichGO, enricher on the cached KEGG table) against the targets of every tested miRNA as universe.
# Outputs: results/de_mirna_targets.csv, results/target_enrichment_{GO,KEGG}.csv, plots/target_enrichment.png,
# reports/sessionInfo.txt.
# ============================================================
source(file.path(snakemake@scriptdir, "00_common.R"))
suppressPackageStartupMessages({
  library(clusterProfiler)
  library(enrichplot)
  library(cfg$enrichment$orgdb, character.only = TRUE)
})
out <- snakemake@output
inp <- snakemake@input
ecfg <- cfg$enrichment
orgdb <- get(ecfg$orgdb)

tt <- read_csv(inp$de, show_col_types = FALSE)
targets <- read_tsv(inp$targets, col_types = cols(.default = "c"))
hits <- de_hits(tt)
rule <- deg_criterion(tt)
message("DE miRNAs (", rule$col, " < ", rule$cut, ", |log2FC| > ", cfg$de$lfc, "): ", nrow(hits))

universe <- unique(targets$target_entrez)
message("Universe (targets of tested miRNAs): ", length(universe), " genes")

de_targets <- targets |>
  inner_join(hits |> dplyr::select(mirna, log2FoldChange, pvalue, padj), by = "mirna") |>
  mutate(direction = ifelse(log2FoldChange > 0, "up", "down"))
write_csv(de_targets, out$de_targets)

# Gene sets: targets of the up and of the down miRNAs (a gene counts when >= min_mirnas_per_target hit it).
min_mir <- ecfg$min_mirnas_per_target %||% 1
gene_sets <- lapply(split(de_targets, de_targets$direction), function(d) {
  d |> distinct(mirna, target_entrez) |> count(target_entrez) |> filter(n >= min_mir) |> pull(target_entrez)
})
gene_sets <- gene_sets[lengths(gene_sets) > 0]
message("Target genes per direction: ", paste(names(gene_sets), lengths(gene_sets), collapse = ", "))

go_cols <- c("direction", "ID", "Description", "GeneRatio", "BgRatio", "pvalue", "p.adjust", "qvalue",
             "geneID", "Count")
kegg_link <- read_tsv(inp$kegg_link, col_names = c("pathway", "gene"), col_types = "cc") |>
  mutate(pathway = sub("^path:", "", pathway), gene = sub("^[a-z]+:", "", gene))
kegg_names <- read_tsv(inp$kegg_names, col_names = c("pathway", "name"), col_types = "cc") |>
  mutate(pathway = sub("^path:", "", pathway), name = sub(" - [^-]+\\([^)]*\\)$", "", name))

run_ora <- function(fun) {
  res <- lapply(names(gene_sets), function(dir) {
    r <- tryCatch(fun(gene_sets[[dir]]), error = function(e) {
      message("ORA failed (", dir, "): ", conditionMessage(e)); NULL })
    if (is.null(r) || !nrow(as.data.frame(r))) return(NULL)
    r <- setReadable(r, OrgDb = orgdb, keyType = "ENTREZID")
    as_tibble(as.data.frame(r)) |> mutate(direction = dir, .before = 1)
  })
  bind_rows(res)
}
go_fun <- function(g) enrichGO(g, OrgDb = orgdb, keyType = "ENTREZID", ont = "BP", universe = universe,
                               pvalueCutoff = 1, qvalueCutoff = 1, pAdjustMethod = "BH",
                               minGSSize = ecfg$min_gs_size, maxGSSize = ecfg$max_gs_size)
kegg_fun <- function(g) enricher(g, TERM2GENE = kegg_link[, c("pathway", "gene")],
                                 TERM2NAME = kegg_names, universe = universe, pvalueCutoff = 1,
                                 qvalueCutoff = 1, pAdjustMethod = "BH",
                                 minGSSize = ecfg$min_gs_size, maxGSSize = ecfg$max_gs_size)

if (length(gene_sets)) {
  go <- run_ora(go_fun)
  kegg <- run_ora(kegg_fun)
} else {
  go <- kegg <- tibble()
}
if (nrow(go)) write_csv(go |> arrange(direction, p.adjust), out$go) else write_empty(out$go, go_cols, "GO results")
if (nrow(kegg)) write_csv(kegg |> arrange(direction, p.adjust), out$kegg) else write_empty(out$kegg, go_cols, "KEGG results")
for (nm in c("go", "kegg")) {
  d <- get(nm)
  if (nrow(d)) message(toupper(nm), " terms with padj < ", ecfg$pvalue_cutoff, ": ",
                       paste(names(table(d$direction[d$p.adjust < ecfg$pvalue_cutoff])),
                             table(d$direction[d$p.adjust < ecfg$pvalue_cutoff]), collapse = ", "))
}

# <figure name="target_enrichment">
# Top GO BP and KEGG terms per direction (by adjusted p), gene ratio and adjusted p
top_terms <- bind_rows(
  if (nrow(go)) go |> mutate(source = "GO BP") else NULL,
  if (nrow(kegg)) kegg |> mutate(source = "KEGG") else NULL
)
if (nrow(top_terms)) {
  top_terms <- top_terms |>
    mutate(ratio = sapply(strsplit(GeneRatio, "/"), function(x) as.numeric(x[1]) / as.numeric(x[2]))) |>
    group_by(source, direction) |> slice_min(p.adjust, n = 10, with_ties = FALSE) |> ungroup() |>
    mutate(Description = stringr::str_trunc(Description, 60),
           label = reorder(paste(Description, direction, sep = "___"), -log10(p.adjust)))
  save_plot(out$dotplot,
            ggplot(top_terms, aes(x = ratio, y = label, size = Count, colour = -log10(p.adjust))) +
              geom_point() +
              scale_y_discrete(labels = function(x) sub("___.*$", "", x)) +
              scale_colour_gradient(low = "#3399FF", high = "#CC0000") +
              facet_grid(source ~ direction, scales = "free_y", space = "free_y") +
              labs(x = "Gene ratio", y = NULL, colour = "-log10 padj", size = "Genes") +
              theme_bw(base_size = 10),
            width = 11, height = 2 + 0.22 * nrow(top_terms))
} else {
  placeholder_plot(out$dotplot, "No DE miRNA with validated targets, or no enriched term")
}
# </figure>

writeLines(capture.output(sessionInfo()), out$session)
