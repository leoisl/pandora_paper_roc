rule merge_precision_and_recall_dfs:
    input:
         precision_file = output_folder + "/precision/precision_files/{coverage}/{tool}/coverage_filter_{coverage_threshold}/strand_bias_filter_{strand_bias_threshold}/gaps_filter_{gaps_threshold}/precision.tsv",
         recall_file = output_folder + "/recall/recall_files/{coverage}/{tool}/coverage_filter_{coverage_threshold}/strand_bias_filter_{strand_bias_threshold}/gaps_filter_{gaps_threshold}/recall.tsv"
    output:
         error_rate_and_recall_file = output_folder + "/plot_data/{coverage}/{tool}/coverage_filter_{coverage_threshold}/strand_bias_filter_{strand_bias_threshold}/gaps_filter_{gaps_threshold}/ROC_data.tsv"
    threads: 1
    resources:
        mem_mb = lambda wildcards, attempt: 4000 * attempt
    log:
        "logs/merge_precision_and_recall_dfs/{coverage}/{tool}/coverage_filter_{coverage_threshold}/strand_bias_filter_{strand_bias_threshold}/gaps_filter_{gaps_threshold}/ROC_data.log"
    script:
        "../scripts/merge_precision_and_recall_dfs.py"


rule concat_all_plot_data:
    input:
         all_plot_data_intermediate_files = all_plot_data_intermediate_files
    output:
         final_plot_data_file = output_folder + "/plot_data/ROC_data.tsv"
    threads: 1
    resources:
        mem_mb = lambda wildcards, attempt: 4000 * attempt
    log:
        "logs/concat_all_plot_data/ROC_data.log"
    script:
        "../scripts/concat_all_plot_data.py"


rule concat_all_nb_of_records_removed_with_mapq_sam_records_filter_files_for_precision:
    input:
         all_nb_of_records_removed_with_mapq_sam_records_filter_files_for_precision = all_nb_of_records_removed_with_mapq_sam_records_filter_files_for_precision
    output:
         nb_of_records_removed_with_mapq_sam_records_filter_for_precision_filepath = output_folder + "/plot_data/nb_of_records_removed_with_mapq_sam_records_filter_for_precision.csv"
    threads: 1
    resources:
        mem_mb = lambda wildcards, attempt: 4000 * attempt
    log:
        "logs/concat_all_nb_of_records_removed_with_mapq_sam_records_filter_files_for_precision/nb_of_records_removed_with_mapq_sam_records_filter_for_precision.log"
    script:
        "../scripts/concat_all_nb_of_records_removed_with_mapq_sam_records_filter_files.py"


rule concat_all_recall_per_sample_no_gt_conf_filter:
    input:
         all_recall_per_sample_no_gt_conf_filter = all_recall_per_sample_no_gt_conf_filter
    output:
         recall_per_sample = output_folder + "/plot_data/recall_per_sample/recall_per_sample.tsv"
    threads: 1
    resources:
        mem_mb = lambda wildcards, attempt: 4000 * attempt
    log:
        "logs/concat_all_recall_per_sample_no_gt_conf_filter/recall_per_sample.log"
    script:
        "../scripts/concat_all_recall_per_sample_no_gt_conf_filter.py"


rule concat_all_recall_per_sample_per_nb_of_samples:
    input:
         all_recalls_per_sample_per_nb_of_samples = cov_tool_and_filters_recall_per_sample_per_number_of_samples.values()
    output:
         aggregated_recall_per_sample_per_nb_of_samples = output_folder + "/plot_data/recall_per_sample_per_number_of_samples.tsv"
    threads: 1
    resources:
        mem_mb = lambda wildcards, attempt: 4000 * attempt
    log:
        "logs/concat_all_recall_per_sample_per_nb_of_samples.log"
    run:
        import pandas as pd
        aggregated_df = pd.DataFrame(columns=["sample", "coverage", "tool", "coverage_threshold", "strand_bias_threshold",
                                   "gaps_threshold", "nb_of_samples", "recalls_wrt_truth_probes",
                                   "nbs_of_truth_probes_found", "nbs_of_truth_probes_in_total"])

        for (sample, coverage, tool, coverage_threshold, strand_bias_threshold, gaps_threshold), df_filepath \
            in cov_tool_and_filters_recall_per_sample_per_number_of_samples.items():
            df = pd.read_csv(df_filepath)
            aggregated_df = pd.concat([aggregated_df, df], ignore_index=True)

        aggregated_df.to_csv(output.aggregated_recall_per_sample_per_nb_of_samples, index=False, sep="\t")


rule aggregate_recall_per_number_of_samples:
    input:
         all_recalls_per_number_of_samples = cov_tool_and_filters_recall_per_number_of_samples.values()
    output:
         aggregated_recall_per_number_of_samples = output_folder + "/plot_data/recall_per_number_of_samples.tsv"
    threads: 1
    resources:
        mem_mb = lambda wildcards, attempt: 4000 * attempt
    log:
        "logs/aggregate_recall_per_number_of_samples.log"
    run:
        import pandas as pd

        aggregated_df = pd.DataFrame(columns=["coverage", "tool", "coverage_threshold", "strand_bias_threshold",
                                   "gaps_threshold", "NB_OF_SAMPLES", "recall"])

        for (coverage, tool, coverage_threshold, strand_bias_threshold, gaps_threshold), df_filepath \
            in cov_tool_and_filters_recall_per_number_of_samples.items():
            df = pd.read_csv(df_filepath)
            df["coverage"] = coverage
            df["tool"] = tool
            df["coverage_threshold"] = coverage_threshold
            df["strand_bias_threshold"] = strand_bias_threshold
            df["gaps_threshold"] = gaps_threshold
            aggregated_df = pd.concat([aggregated_df, df], ignore_index=True)

        aggregated_df.to_csv(output.aggregated_recall_per_number_of_samples, index=False, sep="\t")


rule concat_all_precision_per_sample_no_gt_conf_filter:
    input:
         all_precision_per_sample_no_gt_conf_filter = all_precision_per_sample_no_gt_conf_filter
    output:
         precision_per_sample = output_folder + "/plot_data/precision_per_sample.tsv"
    threads: 1
    resources:
        mem_mb = lambda wildcards, attempt: 4000 * attempt
    log:
        "logs/concat_all_precision_per_sample_no_gt_conf_filter.log"
    script:
        "../scripts/concat_all_precision_per_sample_no_gt_conf_filter.py"


rule make_enrichment_of_FPs_per_sample_plot:
    input:
        precision_per_sample = rules.concat_all_precision_per_sample_no_gt_conf_filter.output.precision_per_sample
    output:
        csv_data = output_folder + "/plot_data/enrichment_of_FPs/enrichment_of_FPs.csv",
        plot = output_folder + "/plot_data/enrichment_of_FPs/enrichment_of_FPs.png",
    threads: 1
    resources:
        mem_mb=lambda wildcards, attempt: 4000 * attempt
    notebook:
        "../../eda/enrichment_of_FPs_per_sample/enrichment_of_FPs_per_sample.ipynb"



rule make_precision_per_ref_per_clade_csv:
    input:
        precision_per_sample = rules.concat_all_precision_per_sample_no_gt_conf_filter.output.precision_per_sample
    output:
        csv_data = output_folder + "/plot_data/precision_per_ref_per_clade/precision_per_ref_per_clade_{tools_to_keep}.csv"
    threads: 1
    resources:
        mem_mb=lambda wildcards, attempt: 8000 * attempt
    log:
        notebook="logs/make_precision_per_ref_per_clade_csv/{tools_to_keep}.ipynb"
    notebook:
        "../../eda/precision_per_ref_per_clade/precision_per_ref_per_clade.ipynb"

rule make_precision_per_ref_per_clade_plot:
    input:
        csv_data = rules.make_precision_per_ref_per_clade_csv.output.csv_data
    output:
        plot = output_folder + "/plot_data/precision_per_ref_per_clade/precision_per_ref_per_clade_{tools_to_keep}.png"
    threads: 1
    resources:
        mem_mb=lambda wildcards, attempt: 4000 * attempt
    singularity:
        "docker://leandroishilima/pandora1_paper_r:pandora_paper_tag1"
    shell:
        "Rscript eda/precision_per_ref_per_clade/clade_plots.R {input.csv_data} {wildcards.tools_to_keep} {output.plot}"


rule make_recall_per_ref_per_clade_csv:
    input:
        recall_per_sample = rules.concat_all_recall_per_sample_no_gt_conf_filter.output.recall_per_sample
    output:
        csv_data = output_folder + "/plot_data/recall_per_ref_per_clade/recall_per_ref_per_clade_{tools_to_keep}.csv"
    threads: 1
    resources:
        mem_mb=lambda wildcards, attempt: 8000 * attempt
    log:
        notebook="logs/make_recall_per_ref_per_clade_csv/{tools_to_keep}.ipynb"
    notebook:
        "../../eda/recall_per_ref_per_clade/recall_per_ref_per_clade.ipynb"

rule make_recall_per_ref_per_clade_plot:
    input:
        csv_data = rules.make_recall_per_ref_per_clade_csv.output.csv_data
    output:
        plot = output_folder + "/plot_data/recall_per_ref_per_clade/recall_per_ref_per_clade_{tools_to_keep}.png"
    threads: 1
    resources:
        mem_mb=lambda wildcards, attempt: 4000 * attempt
    singularity:
        "docker://leandroishilima/pandora1_paper_r:pandora_paper_tag1"
    shell:
        "Rscript eda/recall_per_ref_per_clade/clade_plots.R {input.csv_data} 0.75 {wildcards.tools_to_keep} {output.plot}"


rule make_recall_per_ref_per_nb_of_samples_per_clade_csv:
    input:
        aggregated_recall_per_sample_per_nb_of_samples = rules.concat_all_recall_per_sample_per_nb_of_samples.output.aggregated_recall_per_sample_per_nb_of_samples
    output:
        csv_data = expand(output_folder + "/plot_data/recall_per_ref_per_nb_of_samples_per_clade/recall_per_ref_per_nb_of_samples_per_clade.{{tools_to_keep}}.nb_of_samples_{nb_of_samples}.csv", nb_of_samples = list_with_number_of_samples)
    params:
        list_with_number_of_samples = list_with_number_of_samples,
        output_file_string_format = lambda wildcards: output_folder + f"/plot_data/recall_per_ref_per_nb_of_samples_per_clade/recall_per_ref_per_nb_of_samples_per_clade.{wildcards.tools_to_keep}.nb_of_samples_{{nb_of_samples}}.csv"
    threads: 1
    resources:
        mem_mb=lambda wildcards, attempt: 8000 * attempt
    log:
        notebook="logs/make_recall_per_ref_per_nb_of_samples_per_clade_csv/{tools_to_keep}.ipynb"
    notebook:
        "../../eda/recall_per_ref_per_nb_of_samples_per_clade/recall_per_ref_per_nb_of_samples_per_clade.ipynb"

rule make_recall_per_ref_per_nb_of_samples_per_clade_plot:
    input:
        csv_data = rules.make_recall_per_ref_per_nb_of_samples_per_clade_csv.output.csv_data
    output:
        plots = expand(output_folder + "/plot_data/recall_per_ref_per_nb_of_samples_per_clade/recall_per_ref_per_nb_of_samples_per_clade.{{tools_to_keep}}.nb_of_samples_{nb_of_samples}.png", nb_of_samples = list_with_number_of_samples),
        gif = output_folder + "/plot_data/recall_per_ref_per_nb_of_samples_per_clade/recall_per_ref_per_nb_of_samples_per_clade.{tools_to_keep}.gif"
    params:
        list_with_number_of_samples = list_with_number_of_samples,
    threads: 1
    resources:
        mem_mb=lambda wildcards, attempt: 4000 * attempt
    singularity:
        "docker://leandroishilima/pandora1_paper_r:pandora_paper_tag1"
    script:
        "../scripts/make_recall_per_ref_per_nb_of_samples_per_clade_plot.py"


rule make_recall_per_sample_plot:
    input:
         recall_per_sample = rules.concat_all_recall_per_sample_no_gt_conf_filter.output.recall_per_sample
    output:
         lineplot = output_folder + "/plot_data/recall_per_sample/recall_per_sample_{tools_to_keep}.lineplot.png",
         boxplot = output_folder + "/plot_data/recall_per_sample/recall_per_sample_{tools_to_keep}.boxplot.png",
    threads: 1
    resources:
        mem_mb=lambda wildcards, attempt: 4000 * attempt
    log:
        notebook="logs/make_recall_per_sample_plot/{tools_to_keep}.ipynb"
    notebook:
        "../../eda/recall_per_sample/recall_per_sample.ipynb"
