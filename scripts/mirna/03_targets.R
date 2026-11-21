# ============================================================
# SCRIPT 03: validated targets of every tested miRNA (multiMiR), cached
# Queries the multiMiR server (validated tables: miRTarBase, TarBase, miRecords) once for all miRNAs that
# passed the filter, so the enrichment universe is "targets of any tested miRNA". Needs internet.
# Output: cache/multimir/validated_targets.tsv (one row per miRNA × target gene × database).
# ============================================================
source(file.path(snakemake@scriptdir, "00_common.R"))
suppressPackageStartupMessages(library(multiMiR))
out <- snakemake@output
tcfg <- cfg$targets

tested <- read_csv(snakemake@input$logcpm, col_types = cols(mirna = "c", .default = "d"))$mirna
message("miRNAs tested: ", length(tested))
dir.create(dirname(out$targets), showWarnings = FALSE, recursive = TRUE)

query <- function(mirs, tries = 3) {
  for (k in seq_len(tries)) {
    res <- tryCatch(
      get_multimir(org = "hsa", mirna = mirs, table = "validated", summary = FALSE, legacy.out = FALSE),
      error = function(e) { message("multiMiR attempt ", k, " failed: ", conditionMessage(e)); NULL })
    if (!is.null(res)) return(res@data)
    Sys.sleep(20)
  }
  stop("multiMiR query failed ", tries, " times")
}

batches <- split(tested, ceiling(seq_along(tested) / 25))
raw <- bind_rows(lapply(seq_along(batches), function(i) {
  message("Batch ", i, "/", length(batches))
  d <- query(batches[[i]])
  if (is.null(d) || !nrow(d)) return(NULL)
  as_tibble(d) |> mutate(across(everything(), as.character))
}))
message("multiMiR rows: ", nrow(raw))

dbs <- tolower(unlist(tcfg$databases))
targets <- raw |>
  filter(tolower(database) %in% dbs, !is.na(target_entrez), target_entrez != "") |>
  filter(!(isTRUE(tcfg$mirtarbase_strong_only) & tolower(database) == "mirtarbase" &
           grepl("weak", support_type, ignore.case = TRUE))) |>
  group_by(mirna = mature_mirna_id, target_symbol, target_entrez, database) |>
  summarise(experiments = paste(sort(unique(na.omit(experiment))), collapse = "; "),
            support_type = paste(sort(unique(na.omit(support_type))), collapse = "; "),
            n_pubmed = n_distinct(pubmed_id), .groups = "drop") |>
  arrange(mirna, target_symbol, database)

write_tsv(targets, out$targets)
message("Validated miRNA-target pairs: ", nrow(distinct(targets, mirna, target_entrez)),
        " | miRNAs with a target: ", n_distinct(targets$mirna), " of ", length(tested),
        " | target genes: ", n_distinct(targets$target_entrez))
print(targets |> count(database))
