rule make_variant_calls_probeset:
    input:
         vcf = lambda wildcards: data.xs((wildcards.sample_id, wildcards.coverage, wildcards.tool))["vcf"],
         vcf_ref = lambda wildcards: data.xs((wildcards.sample_id, wildcards.coverage, wildcards.tool))["vcf_reference"]
    output:
          probeset = "analysis/variant_calls_probesets/{sample_id}/{coverage}/{tool}.variant_calls_probeset.fa"
    params:
          flank_length = config["variant_calls_probe_length"]
    threads: 1
    resources:
        mem_mb = lambda wildcards, attempt: 1000 * attempt
    log:
        "logs/make_variant_calls_probeset/{sample_id}_{coverage}_{tool}.log"
    script:
        "../scripts/make_variant_calls_probeset.py"


rule merge_precision_and_recall_dfs:
    input:
         precision_file = precision_file,
         recall_file = recall_file
    output:
         error_rate_and_recall_file = "analysis/plot/error_rate_and_recall_gt_min_{min_gt}_step_{step_gt}_max_{max_gt}.tsv"
    threads: 1
    resources:
        mem_mb = lambda wildcards, attempt: 2000 * attempt
    log:
        "logs/merge_precision_and_recall_dfs/gt_min_{min_gt}_step_{step_gt}_max_{max_gt}.log"
    script:
        "../scripts/merge_precision_and_recall_dfs.py"
