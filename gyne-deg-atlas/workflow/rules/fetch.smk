# Pull analysis-ready expression for one dataset.
# RNA-seq -> recount3/ARCHS4 counts ; Microarray -> GEO series matrix.
# NO FASTQ, NO ALIGNMENT here (8 GB constraint).

def _src(ds): return samples.loc[ds, "data_source"]

rule fetch_counts:
    output: "data/expr/{dataset}.expr.rds"
    params:
        technique=lambda wc: samples.loc[wc.dataset, "technique"],
        source=lambda wc: samples.loc[wc.dataset, "data_source"],
    conda: "../../envs/deseq2.yaml"
    resources: mem_mb=3000
    log: "results/logs/fetch_{dataset}.log"
    shell:
        r"""
        Rscript workflow/scripts/fetch_expr.R \
            --dataset {wildcards.dataset} \
            --technique "{params.technique}" \
            --source "{params.source}" \
            --out {output} > {log} 2>&1
        """
