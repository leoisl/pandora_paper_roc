import pandas as pd

def get_concatenated_df(files):
    dfs = [pd.read_csv(file, sep="\t") for file in files]
    concatenated_df = pd.concat(dfs, ignore_index=True)[
       ["tool", "coverage", "coverage_threshold", "strand_bias_threshold", "gaps_threshold", "GT", "error_rate", "recall"]
    ]
    return concatenated_df

# setup
all_plot_data_intermediate_files = snakemake.input.all_plot_data_intermediate_files
output_file = snakemake.output.final_plot_data_file

# read
all_dfs = get_concatenated_df(all_plot_data_intermediate_files)

# output
all_dfs.to_csv(output_file, sep="\t")
