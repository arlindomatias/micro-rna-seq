# Plasma microRNA sequencing

Snakemake + R workflows for circulating miRNA (small RNA-seq) of public GEO datasets on chemotherapy side
effects. Same layout and conventions as [bulk-rna-seq](https://github.com/arlindomatias/bulk-rna-seq).

| Project | Dataset | Description |
|---|---|---|
| [`GSE260661/`](GSE260661) | GSE260661 | Plasma of 8 colorectal cancer patients on FOLFOX, sampled at 2–3 and 6–8 cycles; paired edgeR, CIPN grade, validated-target enrichment. [Report](https://arlindomatias.github.io/micro-rna-seq/GSE260661/reports/report.html) |

Shared scripts: `scripts/html_report.py` (report builder, copy of the bulk-rna-seq one) and
`scripts/mirna/` (prepare, edgeR, targets, enrichment, report).
