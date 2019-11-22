rule merge_precision_and_recall_dfs:
    input:
         all_precision_files = all_precision_files,
         all_recall_files = all_recall_files
    output:
         error_rate_and_recall_file = "analysis/plot/error_rate_and_recall_gt_min_{min_gt}_step_{step_gt}_max_{max_gt}.tsv"
    threads: 1
    resources:
        mem_mb = lambda wildcards, attempt: 2000 * attempt
    log:
        "logs/merge_precision_and_recall_dfs/gt_min_{min_gt}_step_{step_gt}_max_{max_gt}.log"
    script:
        "../scripts/merge_precision_and_recall_dfs.py"


rule generate_ROC_curve:
    input:
        error_rate_and_recall_file = rules.merge_precision_and_recall_dfs.output.error_rate_and_recall_file
    output:
        plot = "analysis/plot/error_rate_and_recall_gt_min_{min_gt}_step_{step_gt}_max_{max_gt}.pdf"
    threads: 1
    resources:
        mem_mb = lambda wildcards, attempt: 2000 * attempt
    log:
        "logs/generate_ROC_curve/gt_min_{min_gt}_step_{step_gt}_max_{max_gt}.log"
    script:
        "../scripts/generate_ROC_curve.py"
