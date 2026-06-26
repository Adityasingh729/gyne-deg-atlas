# limma on GEO processed expression (microarray strata).
MICROARRAY_DATASETS = [d for d in DATASETS if samples.loc[d, "technique"] == "Microarray"]

rule deg_microarray:
    input: "data/expr/{dataset}.expr.rds"
    output: "results/deg/{dataset}.deg.tsv"
    wildcard_constraints:
        dataset="|".join(MICROARRAY_DATASETS) if MICROARRAY_DATASETS else "NONE"
    params:
        group=lambda wc: samples.loc[wc.dataset, "group_column"],
        case=lambda wc: samples.loc[wc.dataset, "case_label"],
        ctrl=lambda wc: samples.loc[wc.dataset, "control_label"],
    conda: "../../envs/limma.yaml"
    resources: mem_mb=3000
    log: "results/logs/deg_micro_{dataset}.log"
    shell:
        "Rscript workflow/scripts/deg_limma.R --expr {input} --group \"{params.group}\" --case \"{params.case}\" --control \"{params.ctrl}\" --out {output} > {log} 2>&1"


 