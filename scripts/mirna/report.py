"""HTML report of a miRNA-seq project: design → input/QC → edgeR models → paper candidates → target enrichment.

Run through Snakemake (`script:` directive, rule `report`). Dataset-specific prose comes from
`snakemake.params`; every number is computed from results/ (never hard-coded). Reads only plots/ and results/.
"""
import os
import sys

import pandas as pd

# Shared report module: <repository>/scripts/html_report.py (this script lives in scripts/mirna/).
sys.path.insert(0, os.path.join(snakemake.scriptdir, ".."))  # noqa: F821
from html_report import Report  # noqa: E402

inp, out, prm = snakemake.input, snakemake.output, snakemake.params  # noqa: F821
cfg = snakemake.config  # noqa: F821
P = "plots"
ref, trt = cfg["groups"]["reference"], cfg["groups"]["treatment"]
desc = cfg["groups"].get("descriptions", {})
padj_cut, lfc_cut = float(cfg["de"]["padj"]), float(cfg["de"]["lfc"])
fallback_cut = float(cfg["de"].get("pvalue_fallback", 0.01))
models = cfg["models"]
primary = cfg["primary_model"]
candidates = list(cfg.get("candidates", []))
ecfg = cfg["enrichment"]


def read(path, sep=","):
    try:
        return pd.read_csv(path, sep=sep)
    except (pd.errors.EmptyDataError, FileNotFoundError):
        return pd.DataFrame()


def short(m):
    return m.replace("hsa-", "")


def hits(tt):
    """Same rule as deg_criterion()/de_hits() in 00_common.R."""
    mode = cfg["de"].get("criterion", "auto")
    fb = mode == "pvalue" or (mode == "auto" and not bool((tt["padj"] < padj_cut).any()))
    col, cut = ("pvalue", fallback_cut) if fb else ("padj", padj_cut)
    return tt[(tt[col] < cut) & (tt["log2FoldChange"].abs() > lfc_cut)], col, cut, fb


coldata = read(inp.coldata)
geo_check = read(inp.geo_check, "\t")
lib = read(inp.library_qc, "\t")
summary = read(inp.summary, "\t")
tables = {m: read(f"results/de_{m}.csv") for m in models}
cand = read(inp.candidates, "\t")
targets = read(inp.targets, "\t")
de_targets = read(inp.de_targets)
go = read(inp.go)
kegg = read(inp.kegg)

n_pat = coldata["patient"].nunique()
tested = len(tables[primary])
prim_hits, prim_col, prim_cut, prim_fb = hits(tables[primary])

r = Report(prm.title, subtitle=prm.subtitle, code_files=list(inp.code))
r.kpis([
    (f"{n_pat} × 2", "patients × draws"),
    (f"{len(coldata)}", "libraries"),
    (f"{tested:,}", "miRNAs tested"),
    (f"{int((prim_hits['log2FoldChange'] > 0).sum())} / {int((prim_hits['log2FoldChange'] < 0).sum())}",
     f"up / down, {primary} model" + (f" (p < {prim_cut:g})" if prim_fb else "")),
    (f"{int((tables[primary]['padj'] < padj_cut).sum())}", f"miRNAs at FDR < {padj_cut:g} ({primary})"),
])

# ---------------------------------------------------------------- design
r.h2("Study design")
r.md(prm.design)
r.md(f"**Source:** {prm.paper}")
r.md(f"* **{ref}**: {desc.get(ref, '')} ({int((coldata['timepoint'] == ref).sum())} samples).\n"
     f"* **{trt}**: {desc.get(trt, '')} ({int((coldata['timepoint'] == trt).sum())} samples).")

per_patient = (coldata.pivot(index="patient", columns="timepoint",
                             values=["folfox_cycles", "ctcae_grade", "oxaliplatin_total_mg"]))
per_patient.columns = [f"{v}_{t}" for v, t in per_patient.columns]
pp = coldata.drop_duplicates("patient").set_index("patient")[["age", "sex", "tnm_stage", "delta_grade"]]
pp = pp.join(per_patient[[f"folfox_cycles_{ref}", f"folfox_cycles_{trt}", f"ctcae_grade_{ref}",
                          f"ctcae_grade_{trt}", f"oxaliplatin_total_mg_{ref}", f"oxaliplatin_total_mg_{trt}"]])
r.table(pp.reset_index(), "Patients (curated from Table 1 of the paper, `data/samples.tsv`): FOLFOX cycles, "
        "NCI-CTCAE CIPN grade and cumulative oxaliplatin dose (mg) at each draw; `delta_grade` = late − early.",
        digits=0)
dg = pp["delta_grade"].value_counts().sort_index()
r.md("CIPN grade change between the draws: " + ", ".join(f"**{int(k):+d}** in {v} patients" for k, v in dg.items())
     + ". Cycles, cumulative dose and grade all rise between the two draws, so they cannot be separated.")

r.h3("GEO metadata")
bad = geo_check[~(geo_check["cycles_match"] & geo_check["grade_match"])]
r.md(f"""
The sample characteristics in GEO do not match the paper for **{len(bad)} of {len(geo_check)} samples**
(cycles differ in {int((~geo_check['cycles_match']).sum())}, grade in {int((~geo_check['grade_match']).sum())}).
GEO lists the late draws almost uniformly as "7 cycles & grade 2". This analysis uses the values of the paper's
Table 1 (`data/samples.tsv`); the comparison is in `results/geo_metadata_check.tsv`.
""")
if len(bad):
    r.table(bad[["sample_id", "gsm", "paper_cycles", "geo_cycles", "paper_grade", "geo_grade"]],
            "Samples whose GEO characteristics differ from the paper.", digits=0)

r.h3("What the authors did")
r.md(prm.paper_analysis)

# ---------------------------------------------------------------- pipeline
if cfg["de"].get("criterion", "auto") == "pvalue":
    de_rule_text = (f"DE miRNAs: **unadjusted `p < {fallback_cut:g}`** and `|log2FC| > {lfc_cut:g}` in every model "
                    f"(no miRNA reaches FDR < {padj_cut:g} in the primary model). The lists are **exploratory**: "
                    "FDR is reported in every table.")
else:
    de_rule_text = (f"DE miRNAs: `FDR < {padj_cut:g}` and `|log2FC| > {lfc_cut:g}`; when no miRNA reaches FDR < "
                    f"{padj_cut:g}, the unadjusted `p < {fallback_cut:g}` is used instead and the list is flagged as "
                    "**exploratory** (same rule as the bulk RNA-seq projects).")
r.h2("Pipeline")
r.md(f"""
1. **prepare** (`01_prepare_geo.R`): reads the {len(coldata)} GEO count files (raw reads per mature miRNA),
   joins the curated sample sheet and checks it against GEO.
2. **differential** (`02_differential_expression.R`): filter (≥ {cfg['filter']['min_count']} reads in
   ≥ {cfg['filter']['min_samples']} samples), TMM normalization, QC plots and one edgeR quasi-likelihood fit
   per model (`glmQLFit(robust = TRUE)`, `glmQLFTest`).
3. **targets** (`03_targets.R`): validated targets (multiMiR: {', '.join(cfg['targets']['databases'])}) of
   every tested miRNA; miRTarBase limited to strong evidence = {cfg['targets'].get('mirtarbase_strong_only')}.
4. **kegg_annotation**: KEGG `{ecfg['organism_kegg']}` gene–pathway table, downloaded once.
5. **enrichment** (`04_target_enrichment.R`): GO BP and KEGG over-representation of the targets of the DE
   miRNAs of the **{primary}** model, against the targets of all tested miRNAs.
6. **report**: this page.

{de_rule_text} Positive log2FC = higher in the **{trt}** draw.
""")

# ---------------------------------------------------------------- QC
r.h2("Input data and quality control")
r.md(f"""
Libraries have {lib['reads'].min() / 1e6:.1f}–{lib['reads'].max() / 1e6:.1f} M reads mapped to miRBase, but
only {lib['detected'].min()}–{lib['detected'].max()} of the {int(lib['mirnas_annotated'].iloc[0]):,} annotated
miRNAs have at least one read per sample, and {lib['detected_10'].min()}–{lib['detected_10'].max()} have
≥ 10 reads. The ten most abundant miRNAs take {100 * lib['top10_share'].min():.0f}–{100 * lib['top10_share'].max():.0f}%
of each library (the most abundant one is {', '.join(sorted(set(lib['top_mirna'].map(short))))}).
After filtering, **{tested} miRNAs** are tested.
""")
r.table(lib, "Reads, detected miRNAs and dominance of the top miRNAs per library (`results/library_qc.tsv`).",
        digits=3, collapsed=True)
r.figure(f"{P}/library_sizes.png", "Reads and detected miRNAs per library.")
r.figure(f"{P}/top_mirnas.png", "Share of reads taken by the 10 most abundant miRNAs.")
r.h3("Hemolysis")
num, den = cfg["qc"]["hemolysis_numerator"], cfg["qc"]["hemolysis_denominator"]
r.md(f"""
Red blood cells release {short(num)} into plasma; the ratio {short(num)} / {short(den)} is a common hemolysis
indicator (in qPCR, ΔCq > 7 ≈ a 128-fold ratio, log2 = 7). Sequencing ratios are not on the same scale, so
use the plot to spot samples that stand out, and pairs where one draw differs from the other.
""")
r.figure(f"{P}/hemolysis.png", f"log2 ratio of CPM {short(num)} / {short(den)} (pseudocount 1) per draw; "
         "crosses mark samples with no read of the denominator.", width="70%")
r.h3("Sample similarity")
r.figure(f"{P}/density_logcpm.png", "Distribution of log2 CPM (TMM) of the tested miRNAs.", width="80%")
r.figure(f"{P}/mds.png", "MDS of log2 CPM; lines join the two draws of the same patient.", width="80%")
r.figure(f"{P}/sample_correlation.png", "Pearson correlation of log2 CPM centred per miRNA.", width="85%")

# ---------------------------------------------------------------- models
r.h2("Differential expression")
r.md("Every model is fitted on the same filtered, TMM-normalized counts.")
show = summary[["model", "formula", "coef", "mirnas_tested", "common_bcv", "prior_df", "padj_only",
                "min_padj", "de_criterion", "de_up", "de_down", "p_below_fallback", "expected_by_chance"]]
r.table(show, f"Models (`results/de_summary.tsv`). `padj_only` = miRNAs at FDR < {padj_cut:g}; "
        f"`p_below_fallback` = miRNAs at p < {fallback_cut:g}, against `expected_by_chance` under the null.",
        digits=3)
r.h3("Model fit")
r.figure(f"{P}/bcv.png", f"Biological coefficient of variation ({primary} model).", width="70%")
r.figure(f"{P}/ql_dispersion.png", f"Quasi-likelihood dispersions ({primary} model).", width="70%")

for m, spec in models.items():
    tt = tables[m]
    h, col, cut, fb = hits(tt)
    s = summary[summary["model"] == m].iloc[0]
    r.h3(f"Model: {m}")
    r.md(f"`{spec['formula']}`, coefficient `{spec['coef']}`. {spec.get('description', '')}")
    if fb:
        r.md(f"""
> **{int(s['padj_only']) if s['padj_only'] else 'No'} miRNA{'' if s['padj_only'] == 1 else 's'} reach{'es' if s['padj_only'] == 1 else ''} FDR < {padj_cut:g}** (smallest FDR: {s['min_padj']:.3g}); the DE list uses the unadjusted p-value. {int(s['p_below_fallback'])}
> miRNAs have p < {cut:g} where {s['expected_by_chance']:g} are expected by chance; with the fold-change cutoff,
> **{len(h)}** remain ({int((h['log2FoldChange'] > 0).sum())} up, {int((h['log2FoldChange'] < 0).sum())} down).
> Treat them as exploratory.
""")
    else:
        r.md(f"**{len(h)} miRNAs** pass FDR < {padj_cut:g} and |log2FC| > {lfc_cut:g}: "
             f"{int((h['log2FoldChange'] > 0).sum())} up and {int((h['log2FoldChange'] < 0).sum())} down.")
    r.figure(f"{P}/volcano_{m}.png", f"Volcano plot, {m} model (top 10 and the paper's candidates labelled).",
             code=["volcano"])
    r.figure(f"{P}/pvalue_histogram_{m}.png", f"p-value histogram, {m} model. A flat histogram means "
             "no detectable signal; a peak near 0 means some miRNAs change.", width="60%",
             code=["pvalue_histogram"])
    cols = ["mirna", "logCPM", "log2FoldChange", "pvalue", "padj", f"mean_cpm_{ref}", f"mean_cpm_{trt}", "detected"]
    r.table(tt[cols].head(10), f"Top 10 miRNAs by p-value, {m} model (`results/de_{m}.csv`).", digits=3)

r.h3("Top miRNAs within patients")
r.figure(f"{P}/heatmap.png", f"Top {min(tested, cfg['plots'].get('heatmap_max_mirnas', 40))} miRNAs of the "
         f"{primary} model: log2 CPM centred within each patient (each patient's two draws sum to zero), "
         "so the colours show the within-patient change.")
r.h3("Paired vs unpaired")
r.md(f"""
The {primary} model removes the between-patient differences from the residual variance. Where it gains
power, miRNAs lie above the diagonal on the x-axis side; where the patient effect is weak, the two models agree.
""")
r.figure(f"{P}/model_comparison.png", f"-log10 p of each miRNA in the {primary} model vs the other models; "
         "the paper's candidates in red.")

# ---------------------------------------------------------------- candidates
r.h2("The paper's candidate miRNAs")
lines = []
for _, c in cand.iterrows():
    nm = short(c["mirna"])
    if not c["in_data"]:
        lines.append(f"* **{nm}** is not in the count table.")
        continue
    base = (f"* **{nm}**: {int(c['reads_total']):,} reads in total, detected in {int(c['detected'])} of "
            f"{len(coldata)} libraries (≥ {cfg['filter']['min_count']} reads in {int(c['samples_ge_min_count'])}); "
            f"higher in the {trt} draw in {int(c['patients_up'])} of {n_pat} patients "
            f"(median paired change {c['median_paired_log2_change']:+.2f} log2 CPM, paired Wilcoxon "
            f"p = {c['paired_wilcoxon_p']:.3g}).")
    if c["tested"]:
        base += " edgeR: " + "; ".join(
            f"{m} log2FC {c[f'{m}_log2FC']:+.2f}, p = {c[f'{m}_pvalue']:.3g}, FDR = {c[f'{m}_padj']:.3g}"
            for m in models) + "."
    else:
        base += " Filtered out before edgeR (too few reads)."
    lines.append(base)
r.md("\n".join(lines))
r.figure(f"{P}/candidates.png", "Candidate miRNAs per patient at both draws, log2(CPM + 1) before filtering.")
r.table(cand, "Candidate miRNAs in every model (`results/candidates.tsv`).", digits=3, collapsed=True)

# ---------------------------------------------------------------- targets
r.h2("Validated targets and enrichment")
n_univ = targets["target_entrez"].nunique() if len(targets) else 0
n_mir_t = targets["mirna"].nunique() if len(targets) else 0
r.md(f"""
{n_mir_t} of the {tested} tested miRNAs have validated targets in multiMiR, covering **{n_univ:,} genes**
(the enrichment universe). The DE miRNAs of the **{primary}** model ({len(prim_hits)}
{"exploratory, p < " + format(prim_cut, "g") if prim_fb else "at FDR < " + format(prim_cut, "g")}) have
{de_targets['target_entrez'].nunique() if len(de_targets) else 0:,} validated targets
({", ".join(f"{d}: {g['target_entrez'].nunique():,}" for d, g in de_targets.groupby('direction')) if len(de_targets) else "none"}).
Over-representation uses the targets of all tested miRNAs as background, so it asks what is special about the
targets of the changed miRNAs among the targets of plasma miRNAs in general. Validated-target databases are
biased towards well-studied miRNAs and genes, and a target list says nothing about whether the miRNA acts on
it in the tissue that matters (here, peripheral nerve).
""")
if len(de_targets):
    per_mir = (de_targets.groupby(["mirna", "direction"])["target_entrez"].nunique()
               .rename("validated_targets").reset_index().sort_values("validated_targets", ascending=False))
    top = per_mir.iloc[0]
    total = de_targets["target_entrez"].nunique()
    r.table(per_mir, "Validated targets per DE miRNA.", digits=0)
    r.md(f"{short(top['mirna'])} alone contributes {int(top['validated_targets'])} of the {total} target genes "
         f"({100 * top['validated_targets'] / total:.0f}%), so the enrichment mostly describes its known biology.")
for name, d in (("GO BP", go), ("KEGG", kegg)):
    if len(d):
        sig = d[d["p.adjust"] < float(ecfg["pvalue_cutoff"])]
        r.md(f"**{name}:** {len(sig)} terms at adjusted p < {ecfg['pvalue_cutoff']} "
             + (f"({', '.join(f'{k}: {v}' for k, v in sig['direction'].value_counts().items())})." if len(sig) else "."))
        r.table(d.sort_values("p.adjust")[["direction", "ID", "Description", "GeneRatio", "BgRatio", "pvalue",
                                           "p.adjust", "Count"]].head(15),
                f"Top {name} terms (`results/target_enrichment_{'GO' if name == 'GO BP' else 'KEGG'}.csv`).", digits=3)
    else:
        r.md(f"**{name}:** no result (no DE miRNA with validated targets, or no term tested).")
r.figure(f"{P}/target_enrichment.png", "Top GO BP and KEGG terms per direction of the DE miRNAs.")

# ---------------------------------------------------------------- interpretation
r.h2("Interpretation and limitations")
cand_txt = []
for _, c in cand.iterrows():
    if c["in_data"] and c["tested"]:
        cand_txt.append(f"{short(c['mirna'])} (paired p = {c[f'{primary}_pvalue']:.3g}, "
                        f"FDR = {c[f'{primary}_padj']:.3g})")
    elif c["in_data"]:
        cand_txt.append(f"{short(c['mirna'])} (filtered out: ≥ {cfg['filter']['min_count']} reads in only "
                        f"{int(c['samples_ge_min_count'])} samples)")
r.md(f"""
* **No control group.** Every change is relative to the same patient's earlier draw. It reflects
  cumulative FOLFOX exposure, but also time, disease course, surgery recovery and supportive medication, and
  it cannot be attributed to oxaliplatin or to neuropathy specifically.
* **Small and sparse.** {n_pat} patients, {tested} miRNAs with enough reads; the {primary} model finds
  {int((tables[primary]['padj'] < padj_cut).sum())} miRNAs at FDR < {padj_cut:g}.
* **The paper's candidates**: {"; ".join(cand_txt) if cand_txt else "not tested"}.
  The authors used unadjusted p < 0.05 and unpaired groups; with ~2,600 miRNAs some hits are expected by chance.
* **CIPN.** The `cipn_worsening` model asks whether the change is larger in patients whose grade rose more;
  with Δgrade of 0, +1 or +2 in only {n_pat} patients, it can only flag large effects.
""")

# ---------------------------------------------------------------- appendix
outputs = [
    ("data/samples.tsv", "Curated sample sheet (paper Table 1)"),
    ("data/counts.csv", "Raw counts, all miRNAs × samples"),
    ("data/colData.csv", "Sample metadata used in the models"),
    ("results/geo_metadata_check.tsv", "GEO characteristics vs the paper"),
    ("results/library_qc.tsv", "Reads and detection per library"),
    ("results/log2_cpm.csv", "log2 CPM (TMM) of the tested miRNAs"),
    ("results/de_summary.tsv", "One row per model"),
] + [(f"results/de_{m}.csv", f"edgeR QL results, {m} model") for m in models] + [
    ("results/candidates.tsv", "The paper's candidate miRNAs in every model"),
    ("results/de_mirna_targets.csv", f"Validated targets of the DE miRNAs ({primary})"),
    ("results/target_enrichment_GO.csv", "GO BP over-representation of the targets"),
    ("results/target_enrichment_KEGG.csv", "KEGG over-representation of the targets"),
]
r.appendix(outputs, env_files=list(inp.envs), config_files=list(inp.config), session_file=inp.session)
r.write(out.html)
