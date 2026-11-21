# Circulating miRNAs during FOLFOX chemotherapy (GSE260661)

Re-analysis of the public dataset **GSE260661** (NCBI GEO): plasma small RNA-seq of **8 colorectal cancer
patients** on FOLFOX (5-fluorouracil, leucovorin, oxaliplatin), each sampled **twice**: after 2–3 cycles
(*early*) and after 6–8 cycles (*late*). Chemotherapy-induced peripheral neuropathy (CIPN) was graded with
NCI-CTCAE at each draw.

Paper: Ju Y, Baek DH, Choi GE, Jang A. *Exploring microRNA patterns as biomarkers of FOLFOX
chemotherapy-induced peripheral neuropathy in patients with colorectal cancer.* BBA Mol Basis Dis 1870 (2024)
167209, [doi:10.1016/j.bbadis.2024.167209](https://doi.org/10.1016/j.bbadis.2024.167209).

**Report:** [arlindomatias.github.io/micro-rna-seq/GSE260661/reports/report.html](https://arlindomatias.github.io/micro-rna-seq/GSE260661/reports/report.html)

## What is different from the paper

- **No untreated control**: every patient received FOLFOX, so the analysis tests the *within-patient* change
  between the two draws (`~ patient + timepoint`, edgeR quasi-likelihood), not treated vs control. The
  authors compared the draws as independent groups. That model is kept as `unpaired`, for reference.
- **Corrected sample metadata**: GEO lists almost every late draw as "7 cycles & grade 2". The cycles,
  grades and oxaliplatin doses are taken from Table 1 of the paper (`data/samples.tsv`).
  `results/geo_metadata_check.tsv` lists the differences.
- **Multiple-testing correction** (FDR), and a model that relates the change to the CIPN grade change
  (`cipn_worsening`).

## Repository structure

```
.
├── Snakefile                  # Snakemake workflow (scripts in ../scripts/mirna/)
├── config.yaml                # Inputs, groups, models, cutoffs, candidates, targets, enrichment
├── envs/                      # mirna.yaml (+ pin), report.yaml
├── data/
│   ├── GSE260661_RAW/         # GEO count files (one per sample; downloaded if missing)
│   ├── samples.tsv            # Curated sample sheet (paper Table 1)
│   ├── counts.csv, colData.csv
│   └── *.rds                  # DGEList and fits (generated, not versioned)
├── results/                   # Tables: QC, DE per model, candidates, targets, enrichment
├── plots/                     # Figures
├── reports/report.html        # HTML report of the whole analysis
├── docs/                      # Paper PDF (not versioned)
└── cache/                     # multiMiR and KEGG downloads (not versioned)
```

## Pipeline

| Rule | Script (`../scripts/mirna/`) | Main outputs |
|---|---|---|
| `prepare` | `01_prepare_geo.R` | `data/counts.csv`, `data/colData.csv`, `results/library_qc.tsv`, `results/geo_metadata_check.tsv` |
| `differential` | `02_differential_expression.R` | `results/de_<model>.csv`, `results/de_summary.tsv`, `results/candidates.tsv`, `results/log2_cpm.csv`, QC/DE plots |
| `targets` | `03_targets.R` | `cache/multimir/validated_targets.tsv` (miRTarBase strong evidence + miRecords) |
| `kegg_annotation` | (curl) | `cache/kegg/hsa_*.tsv` |
| `enrichment` | `04_target_enrichment.R` | `results/de_mirna_targets.csv`, `results/target_enrichment_{GO,KEGG}.csv` |
| `report` | `report.py` | `reports/report.html` |

```bash
conda activate genomic-variants
snakemake -n --use-conda
snakemake --use-conda --cores 2
```

## Methods summary

- **Input:** raw read counts of 2,588 mature miRNAs (miRBase) × 16 libraries.
- **Filtering:** ≥ 10 reads in ≥ 8 samples, which leaves 119 miRNAs. Normalized with TMM.
- **Models** (edgeR `glmQLFit(robust = TRUE)` + `glmQLFTest`):
  `paired` `~ patient + timepoint` (primary), `unpaired` `~ timepoint`,
  `cipn_worsening` `~ patient + timepoint + late_delta_grade_c`. Positive log2FC means higher in the late draw.
- **DE miRNAs:** **unadjusted p < 0.01** and |log2FC| > 0.5 in every model (`de.criterion: pvalue`), because
  no miRNA passes FDR < 0.05 in the paired model. The lists are exploratory; FDR is in every table.
- **Targets and enrichment:** validated targets of the DE miRNAs; GO BP and KEGG over-representation
  against the targets of all 119 tested miRNAs.

## Results (2026-11-21)

| Model | FDR < 0.05 | Smallest FDR | p < 0.01 (expected by chance) |
|---|---|---|---|
| `paired` | 0 | 0.135 | 3 (1.2): miR-146a-5p, miR-6836-3p, miR-3960, all up |
| `unpaired` | 1 (let-7a-3p, up) | 0.020 | 2 (1.2): let-7a-3p, miR-6836-3p, both up |
| `cipn_worsening` | 0 | 0.878 | 0 (1.2) |

- The paper's candidates go up in most patients but do not survive correction. **miR-3184-5p**: up in
  6 of 8 patients, paired p = 0.077, FDR 0.76. **miR-222-3p**: up in 7 of 8 patients, paired p = 0.011,
  FDR 0.32. Both behave in an on/off way: near zero in some draws and abundant in others.
- Data quality: the dispersion is very high (common BCV ≈ 1.3), and patient S6 differs from all the
  others at both draws (MDS dimension 1, TMM factor 3.1 at the early draw).
- Enrichment (exploratory, from the 3 paired p < 0.01 miRNAs): 74 GO BP and 48 KEGG terms, with
  Toll-like receptor and NF-κB signalling and cytokine responses at the top. Most of the target genes come
  from miR-146a-5p, so these terms mostly reflect what is already known about that miRNA.
