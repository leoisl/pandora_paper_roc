rule merge_precision_and_recall_dfs:
    input:
         precision_file = rules.calculate_precision.output.precision_file,
         recall_file = rules.calculate_recall.output.recall_file
    output:
         error_rate_and_recall_file = "analysis/plot/error_rate_and_recall_gt_min_{min_gt}_step_{step_gt}_max_{max_gt}.tsv"
    threads: 1
    resources:
        mem_mb = lambda wildcards, attempt: 2000 * attempt
    log:
        "logs/merge_precision_and_recall_dfs/gt_min_{min_gt}_step_{step_gt}_max_{max_gt}.log"
    script:
        "../scripts/merge_precision_and_recall_dfs.py"

#
# rule generate_ROC_curve:
#     input:
