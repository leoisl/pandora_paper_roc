from utils import get_concatenated_df

# setup
all_plot_data_intermediate_files = snakemake.input.all_plot_data_intermediate_files
output_file = snakemake.output.final_plot_data_file

# read
all_dfs = get_concatenated_df(all_plot_data_intermediate_files,separator="\t",
                              fields_to_keep=["tool", "coverage", "coverage_threshold", "strand_bias_threshold",
                                              "gaps_threshold", "step_GT", "error_rate", "nb_of_correct_calls",
                                              "nb_of_total_calls", "recall", "nb_of_truth_probes_found",
                                              "nb_of_truth_probes_in_total"])


# some sanity checks
all_nb_of_truth_probes_in_total = all_dfs.nb_of_truth_probes_in_total.unique()
assert len(all_nb_of_truth_probes_in_total) == 1

# output
all_dfs.to_csv(output_file, sep="\t")
