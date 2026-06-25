# Combine per-dataset DEGs WITHIN a (condition x tissue) stratum.
# Robust Rank Aggregation (default) or random-effects (metafor).

def get_stratum_datasets(wildcards):
    datasets_in_stratum = []
    for ds, row in samples.iterrows():
        cond = str(row["condition"]).strip()
        tissue = str(row["tissue_primary"]).strip()
        s_id = f"{cond}__{tissue}".replace(" ", "_").replace("/", "-")
        if s_id == wildcards.stratum:
            datasets_in_stratum.append(ds)
    return [f"results/deg/{ds}.deg.symbols.tsv" for ds in datasets_in_stratum]

rule meta:
    input: get_stratum_datasets
    output: "results/meta/{stratum}.meta.tsv"
    params: method=config["meta"]["method"], minds=config["meta"]["min_datasets"]
    conda: "../../envs/meta.yaml"
    resources: mem_mb=3000
    log: "results/logs/meta_{stratum}.log"
    shell:
        "Rscript workflow/scripts/meta_combine.R --stratum {wildcards.stratum} --method {params.method} --min_datasets {params.minds} --deg_dir results/deg --samplesheet {config[samplesheet]} --out {output} > {log} 2>&1"


