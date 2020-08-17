from utils import get_concatenated_df

# setup
all_precision_per_sample_no_gt_conf_filter = snakemake.input.all_precision_per_sample_no_gt_conf_filter
precision_per_sample = snakemake.output.precision_per_sample

# read
all_dfs = get_concatenated_df(all_precision_per_sample_no_gt_conf_filter, separator="\t")

# output
all_dfs.to_csv(precision_per_sample, sep="\t", index=False)
