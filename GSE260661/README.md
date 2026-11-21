# Circulating miRNAs during FOLFOX chemotherapy (GSE260661)

Re-analysis of the public dataset **GSE260661** (NCBI GEO): plasma small RNA-seq of **8 colorectal cancer
patients** on FOLFOX (5-fluorouracil, leucovorin, oxaliplatin), each sampled **twice**: after 2–3 cycles
(*early*) and after 6–8 cycles (*late*). Chemotherapy-induced peripheral neuropathy (CIPN) was graded with
NCI-CTCAE at each draw.

Paper: Ju Y, Baek DH, Choi GE, Jang A. *Exploring microRNA patterns as biomarkers of FOLFOX
chemotherapy-induced peripheral neuropathy in patients with colorectal cancer.* BBA Mol Basis Dis 1870 (2024)
167209, [doi:10.1016/j.bbadis.2024.167209](https://doi.org/10.1016/j.bbadis.2024.167209).

**Report:** [arlindomatias.github.io/micro-rna-seq/GSE260661/reports/report.html](https://arlindomatias.github.io/micro-rna-seq/GSE260661/reports/report.html)

## The paper

The paper PDF is not kept; this section records the information from it that the analysis uses
or that is needed to read the results.

**Cohort and design.** Between June 2020 and February 2022, 43 plasma samples were collected from 35 patients
with colorectal cancer on FOLFOX at Pusan National University Hospital (ethics approval H-1803-014-065,
CRIS KCT0005353). The patients were split by exposure, because the risk of CIPN rises after ≥ 6 cycles:
group 1 had ≤ 3 cycles and group 2 had ≥ 6 cycles. **Eight patients were sampled twice**, once in each
group: these 16 samples are the sequencing (*analysis*) cohort and the GEO series. The other **27 patients**
(12 with ≤ 3 cycles, 15 with ≥ 6) were sampled once and form the qPCR *validation* cohort, which has no
sequencing data. CIPN was graded by clinicians with NCI-CTCAE v4.03. Blood was drawn in the morning from the
median cubital vein, centrifuged at 3,000 g for 10 min at 4 °C, and the plasma was stored at −20 °C.

**Sequencing cohort (paper Table 1, = `data/samples.tsv`).**

| Patient | Age | Sex | 1st draw: cycles, grade | 2nd draw: cycles, grade | Oxaliplatin (mg), 1st → 2nd | T | N | TNM |
|---|---|---|---|---|---|---|---|---|
| S1 | 62 | F | 2, grade 1 | 7, grade 3 | 242 → 847 | T4 | N1 | IIIB |
| S2 | 68 | F | 3, grade 1 | 7, grade 2 | 342 → 798 | T2 | N0 | I |
| S3 | 56 | F | 3, grade 1 | 8, grade 2 | 378 → 1008 | T3 | N0 | IIA |
| S4 | 58 | M | 2, grade 1 | 6, grade 1 | 314 → 942 | T4a | N0 | IIB |
| S5 | 64 | M | 3, grade 2 | 8, grade 3 | 443 → 1176 | T3 | N2b | IIIC |
| S6 | 71 | M | 3, grade 2 | 7, grade 2 | 471 → 1099 | T2 | N1 | IIIA |
| S7 | 59 | F | 2, grade 1 | 6, grade 2 | 270 → 810 | T4a | N1b | IVB |
| S8 | 61 | M | 3, grade 1 | 8, grade 3 | 501 → 1328 | T4a | N0 | IIB |

The figure legends confirm that `G1_S<n>` and `G2_S<n>` are the same patient.

**Validation cohort (paper Table 2; not in GEO).** Group 1 (n = 12): median age 60 (47–72), 6 F / 6 M,
4 patients at cycle 2 and 8 at cycle 3; CTCAE grade 1 in 9 patients, grade 2 in 3. Group 2 (n = 15): median
age 61 (52–71), 7 F / 8 M, 5 / 6 / 4 patients at cycles 6 / 7 / 8; grade 1 in 2, grade 2 in 8, grade 3 in 5
(p < 0.001). Age, sex, T, N and TNM stage did not differ between the groups. Cumulative oxaliplatin dose:
366.5 ± 82.9 mg vs 974.7 ± 174.5 mg (p < 0.0001).

**Laboratory methods.**
- RNA was extracted with TRIzol and quantified with a NanoDrop 2000.
- Libraries: NEBNext Multiplex Small RNA Library Prep, 75 bp single-end, "NextSeq500" in the paper
  (GEO says NextSeq 550). Reads were trimmed and mapped to miRBase with Bowtie2, then normalized with TMM +
  CPM in edgeR (GEO processing notes); 2,588 mature miRNAs.
- qPCR: miRCURY LNA RT and SYBR Green kits (Qiagen) on a QuantStudio 7 Flex, in triplicate, normalized to
  the spike-in cel-miR-39-3p with 2^−ΔCt. Primer sequences: miR-222-3p `AGCUACAUCUGGCUACUGGGU`,
  miR-3184-5p `UGAGGGGCCUCAGACCGAGCUUUU`, cel-miR-39-3p `UCACCGGGUGUAAAUCAGCUUG`.
- Plasma cortisol was measured by ELISA (Thermo Fisher EIAHCOR), and epinephrine, norepinephrine and
  dopamine by Abnova ELISA kits.
- Statistics: Mann–Whitney U in SPSS 26 / GraphPad Prism 10.

**Differential expression.** DE miRNAs were selected with ExDEGA 4.0.3 (ebiogen): unadjusted p < 0.05 and a
fold change described as "log fold-change > 2", although Table 4 reports linear fold changes. The two draws
were compared as independent groups, without multiple-testing correction. Two miRNAs were reported
(Table 4), both higher in the ≥ 6-cycle draw (negative fold change = group 1 relative to group 2):

| miRNA | Fold change | p |
|---|---|---|
| hsa-miR-222-3p | −2.546 | 0.010 |
| hsa-miR-3184-5p | −2.665 | 0.025 |

**Validation and other findings.**
- In the 27-patient qPCR cohort, miR-222-3p did not differ between the groups. miR-3184-5p was lower in
  the ≤ 3-cycle group (p < 0.05). ROC for CIPN: AUC 0.7806, cut-off < 0.5450, specificity 75.00%
  (95% CI 43–95%), sensitivity 73.33% (45–92%).
- Plasma cortisol was higher with ≥ 6 cycles (134.9 ± 20.0 vs 101.1 ± 25.1 ng/ml, p < 0.001).
  Epinephrine (p = 0.98), norepinephrine (p = 0.59) and dopamine (p = 0.32) did not differ.
- Predicted targets of miR-3184-5p (overlap of miRDB, miRWalk and miRNet): SOX12, CBX6, LASP1, NAV1,
  NACC1, SPATA2, WASF2, ASPH, WIPF2, WFIKKN2, ZBTB46. Enrichr GO terms included Fc-gamma receptor
  signalling, Arp2/3 actin nucleation and the actin cytoskeleton. The KEGG terms were Yersinia and
  pathogenic *E. coli* infection, adherens junction, Fc-gamma R phagocytosis and choline metabolism in
  cancer (WASF2, WIPF2, ASPH).
- The authors conclude that miR-3184-5p and cortisol are putative biomarkers of FOLFOX-induced CIPN.
  Their data-availability statement says "on request"; only the sequencing counts are public (this GEO
  series).

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
