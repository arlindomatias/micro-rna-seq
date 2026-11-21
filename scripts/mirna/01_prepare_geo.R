# ============================================================
# SCRIPT 01: GEO per-sample miRNA count files + curated sample sheet -> count matrix and colData
# Downloads the GEO supplementary tar (one count file per GSM) and the family SOFT file into data/ when they
# are missing. Samples and their clinical variables come from input$samples (curated from the paper), not
# from GEO; the GEO characteristics are only compared with it (results/geo_metadata_check.tsv).
# Outputs: data/counts.csv, data/colData.csv, data/dge_raw.rds (edgeR DGEList), results/library_qc.tsv.
# ============================================================
source(file.path(snakemake@scriptdir, "00_common.R"))
suppressPackageStartupMessages(library(edgeR))
out <- snakemake@output
inp_cfg <- cfg$input
scfg <- cfg$samples

### 1. Download (only when missing) -------------------------------------------
raw_dir <- inp_cfg$raw_dir
if (!dir.exists(raw_dir) || !length(list.files(raw_dir, pattern = "\\.txt\\.gz$"))) {
  tar_file <- file.path(dirname(raw_dir), basename(inp_cfg$raw_tar_url))
  message("Downloading ", inp_cfg$raw_tar_url)
  download.file(inp_cfg$raw_tar_url, tar_file, mode = "wb")
  dir.create(raw_dir, showWarnings = FALSE, recursive = TRUE)
  untar(tar_file, exdir = raw_dir)
  unlink(tar_file)
}
if (!file.exists(inp_cfg$series_soft)) {
  message("Downloading ", inp_cfg$series_soft_url)
  download.file(inp_cfg$series_soft_url, inp_cfg$series_soft, mode = "wb")
}

### 2. Count matrix -------------------------------------------------------------
files <- list.files(raw_dir, pattern = inp_cfg$file_pattern, full.names = TRUE)
ids <- sub(inp_cfg$file_pattern, "\\1", basename(files))
gsm <- sub("_.*", "", basename(files))
message("Count files: ", length(files))

tabs <- lapply(files, function(f) {
  x <- read_tsv(f, col_types = cols(.default = "c"), progress = FALSE)
  stopifnot(ncol(x) == 2, !anyDuplicated(x[[1]]))
  setNames(as.numeric(x[[2]]), x[[1]])
})
mirnas <- names(tabs[[1]])
for (i in seq_along(tabs)) {
  if (!setequal(names(tabs[[i]]), mirnas)) stop("Different miRNA lists in ", files[i])
}
counts <- sapply(tabs, function(v) v[mirnas])
dimnames(counts) <- list(mirnas, ids)
if (any(is.na(counts)) || any(counts != round(counts)) || any(counts < 0)) {
  stop("Counts must be non-negative integers")
}
storage.mode(counts) <- "integer"
message("miRNAs: ", nrow(counts), " | samples: ", ncol(counts))

### 3. Sample sheet -> colData --------------------------------------------------
samples <- read_tsv(snakemake@input$samples, col_types = cols(.default = "c"))
id_col <- scfg$id_column
missing_in_counts <- setdiff(samples[[id_col]], ids)
missing_in_table  <- setdiff(ids, samples[[id_col]])
if (length(missing_in_counts) || length(missing_in_table)) {
  stop("Samples differ between samples.tsv and the count files. Only in table: ",
       paste(missing_in_counts, collapse = ", "), " | only in counts: ", paste(missing_in_table, collapse = ", "))
}
if ("gsm" %in% names(samples)) {
  stopifnot(identical(samples$gsm[match(ids, samples[[id_col]])], gsm))
}

labels <- unlist(cfg$groups$labels)
bad <- labels[make.names(labels) != labels]
if (length(bad)) stop("Group labels must be syntactic R names: ", paste(bad, collapse = ", "))
unknown <- setdiff(unique(samples[[scfg$timepoint_column]]), names(labels))
if (length(unknown)) stop("Timepoints without a label in groups$labels: ", paste(unknown, collapse = ", "))

exclude <- unlist(cfg$exclude_samples)
if (length(exclude)) message("Excluding samples: ", paste(exclude, collapse = ", "))

coldata <- samples |>
  filter(!.data[[id_col]] %in% exclude) |>
  mutate(across(any_of(c("folfox_cycles", "ctcae_grade", "oxaliplatin_total_mg", "age")), as.numeric),
         sample_id = .data[[id_col]],
         patient   = factor(.data[[scfg$patient_column]]),
         timepoint = factor(labels[.data[[scfg$timepoint_column]]], levels = c(ref_level, trt_level)),
         grade     = as.numeric(.data[[scfg$grade_column]])) |>
  arrange(timepoint, patient) |>
  group_by(patient) |>
  mutate(delta_grade = grade[timepoint == trt_level][1] - grade[timepoint == ref_level][1]) |>
  ungroup() |>
  mutate(delta_grade_c = delta_grade - mean(delta_grade[timepoint == ref_level], na.rm = TRUE),
         # Interaction term for `~ patient + timepoint + late_delta_grade_c`: the centred grade change on the
         # late draw, 0 on the early one (the main effect of delta_grade is absorbed by the patient term).
         late_delta_grade_c = ifelse(timepoint == trt_level, delta_grade_c, 0)) |>
  as.data.frame()
rownames(coldata) <- coldata$sample_id
message("Samples per timepoint:")
print(table(coldata$timepoint))
message("Patients with both timepoints: ", sum(table(coldata$patient) == 2), " of ", nlevels(coldata$patient))
message("CIPN grade change per patient (late - early):")
print(table(coldata$delta_grade[coldata$timepoint == ref_level]))

counts <- counts[, coldata$sample_id]
write_csv(data.frame(mirna = rownames(counts), counts, check.names = FALSE), out$counts)
write_csv(coldata, out$coldata)

dge <- DGEList(counts = counts, samples = coldata, group = coldata$timepoint)
saveRDS(dge, out$dge)

### 4. Library QC table ---------------------------------------------------------
lib <- colSums(counts)
top_share <- apply(counts, 2, function(v) max(v) / sum(v))
library_qc <- tibble(
  sample_id   = colnames(counts),
  patient     = as.character(coldata$patient),
  timepoint   = as.character(coldata$timepoint),
  reads       = lib,
  mirnas_annotated = nrow(counts),
  detected    = colSums(counts > 0),
  detected_10 = colSums(counts >= 10),
  top_mirna   = rownames(counts)[apply(counts, 2, which.max)],
  top_share   = round(top_share, 4),
  top10_share = round(apply(counts, 2, function(v) sum(sort(v, decreasing = TRUE)[1:10]) / sum(v)), 4)
)
write_tsv(library_qc, out$library_qc)
print(as.data.frame(library_qc))

### 5. GEO characteristics vs the curated sheet --------------------------------
soft <- readLines(gzfile(inp_cfg$series_soft))
starts <- grep("^\\^SAMPLE = ", soft)
ends <- c(starts[-1] - 1, length(soft))
geo <- bind_rows(lapply(seq_along(starts), function(i) {
  block <- soft[starts[i]:ends[i]]
  ch <- sub("^!Sample_characteristics_ch1 = ", "", grep("^!Sample_characteristics_ch1", block, value = TRUE))
  kv <- setNames(trimws(sub("^[^:]+:", "", ch)), trimws(sub(":.*", "", ch)))
  tibble(gsm = sub("^\\^SAMPLE = ", "", block[1]),
         geo_title = sub("^!Sample_title = ", "", grep("^!Sample_title", block, value = TRUE)),
         geo_cycles = as.numeric(sub("\\D*(\\d+).*", "\\1", kv[["cycle"]] %||% NA)),
         geo_grade = as.numeric(sub("\\D*(\\d+).*", "\\1", kv[["grade"]] %||% NA)))
}))
geo_check <- samples |>
  transmute(sample_id = .data[[id_col]], gsm,
            paper_cycles = as.numeric(folfox_cycles), paper_grade = as.numeric(ctcae_grade)) |>
  left_join(geo, by = "gsm") |>
  mutate(cycles_match = paper_cycles == geo_cycles, grade_match = paper_grade == geo_grade)
write_tsv(geo_check, out$geo_check)
message("GEO characteristics that differ from the paper: cycles ", sum(!geo_check$cycles_match),
        ", grade ", sum(!geo_check$grade_match), " of ", nrow(geo_check), " samples")
