from utils import get_concatenated_df

# setup
all_plot_data_intermediate_files = snakemake.input.all_plot_data_intermediate_files
output_file = snakemake.output.final_plot_data_file

# read
all_dfs = get_concatenated_df(all_plot_data_intermediate_files,separator="\t", fields_to_keep=["tool", "coverage", "coverage_threshold", "strand_bias_threshold", "gaps_threshold", "step_GT", "error_rate", "recall"])

# output
all_dfs.to_csv(output_file, sep="\t")
