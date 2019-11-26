rule merge_precision_and_recall_dfs:
    input:
         precision_file = output_folder + "/precision/precision_files/{coverage}/{tool}/coverage_filter_{coverage_threshold}/strand_bias_filter_{strand_bias_threshold}/gaps_filter_{gaps_threshold}/precision_gt_min_{min_gt}_step_{step_gt}_max_{max_gt}.tsv",
         recall_file = output_folder + "/recall/recall_files/{coverage}/{tool}/coverage_filter_{coverage_threshold}/strand_bias_filter_{strand_bias_threshold}/gaps_filter_{gaps_threshold}/recall_gt_min_{min_gt}_step_{step_gt}_max_{max_gt}.tsv"
    output:
         error_rate_and_recall_file = output_folder + "/plot_data/{coverage}/{tool}/coverage_filter_{coverage_threshold}/strand_bias_filter_{strand_bias_threshold}/gaps_filter_{gaps_threshold}/ROC_data_gt_min_{min_gt}_step_{step_gt}_max_{max_gt}.tsv"
    threads: 1
    resources:
        mem_mb = lambda wildcards, attempt: 2000 * attempt
    log:
        "logs/merge_precision_and_recall_dfs/{coverage}/{tool}/coverage_filter_{coverage_threshold}/strand_bias_filter_{strand_bias_threshold}/gaps_filter_{gaps_threshold}/ROC_data_gt_min_{min_gt}_step_{step_gt}_max_{max_gt}.log"
    script:
        "../scripts/merge_precision_and_recall_dfs.py"


rule merge_all_plot_data:
    input:
         all_plot_data_intermediate_files = all_plot_data_intermediate_files
    output:
         final_plot_data_file =output_folder + "/plot_data/ROC_data_gt_min_{min_gt}_step_{step_gt}_max_{max_gt}.tsv"
    threads: 1
    resources:
        mem_mb = lambda wildcards, attempt: 2000 * attempt
    log:
        "logs/merge_all_plot_data/ROC_data_gt_min_{min_gt}_step_{step_gt}_max_{max_gt}.log"
    script:
        "../scripts/merge_all_plot_data.py"


# TODO : no need for this anymore, we are building interactive plots based on the ROC data
# rule generate_ROC_curve:
#     input:
#         error_rate_and_recall_file = rules.merge_precision_and_recall_dfs.output.error_rate_and_recall_file
#     output:
#         plot = output_folder + "/plot/error_rate_and_recall_gt_min_{min_gt}_step_{step_gt}_max_{max_gt}.pdf"
#     threads: 1
#     resources:
#         mem_mb = lambda wildcards, attempt: 2000 * attempt
#     log:
#         "logs/generate_ROC_curve/gt_min_{min_gt}_step_{step_gt}_max_{max_gt}.log"
#     script:
#         "../scripts/generate_ROC_curve.py"
