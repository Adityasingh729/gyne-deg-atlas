# Gyne-DEG Atlas

Tissue-stratified differential-expression meta-analysis across **PCOS, endometriosis and
adenomyosis**, built from uniformly pre-quantified public datasets. Designed to run on
**8 GB RAM** machines (no genome alignment; pre-computed counts only).

## What this repo does
1. Pulls analysis-ready counts (recount3 / ARCHS4) or GEO processed matrices.
2. Runs per-dataset DEG (DESeq2 for RNA-seq, limma for microarray).
3. Meta-analyses DEGs **within each (condition x tissue) stratum** (random-effects + rank aggregation).
4. Functional enrichment + PPI hub genes.
5. A cross-dataset-validated ML gene panel.
6. Feeds robust hub targets into network pharmacology / docking (separate repos).

## Golden rules (do not break)
- **Never pool across tissues.** One analysis per (condition x tissue) stratum.
- **One primary tissue per dataset.** Auto tissue tags overlap; fix in the sample sheet.
- **RNA-seq and microarray are analysed separately**, then combined at the results level.
- **Human only.** Drop non-human / mixed-species series.
- See `docs/tissue_rules.md` before adding any dataset.

## Quick start
```bash
mamba env create -f environment.yml          # base (snakemake)
conda activate gyne-deg
# per-step envs are created automatically by snakemake --use-conda
snakemake -n                                  # dry run (shows plan)
snakemake --use-conda --cores 2 -p            # 8 GB: keep cores low
```

## Layout
```
config/      sample sheet, strata table, parameters
workflow/    Snakefile + rules/ + scripts/
docs/        DEG schema, tissue rules, contributing
data/        inputs (gitignored)
results/     outputs (gitignored)
```

## 8 GB survival notes
- Use 1-2 cores; counts matrices are small, but R can spike. Call `gc()` often.
- For WGCNA, filter to the top ~5000 variable genes (see `scripts/wgcna.R` TODO).
- Do NOT download FASTQ on these machines. Alignment (if ever needed) runs only on the 64 GB workstation.
