from utils import get_concatenated_df

# setup
all_recall_per_sample_no_gt_conf_filter = snakemake.input.all_recall_per_sample_no_gt_conf_filter
recall_per_sample = snakemake.output.recall_per_sample

# read
all_dfs = get_concatenated_df(all_recall_per_sample_no_gt_conf_filter, separator="\t")

# output
all_dfs.to_csv(recall_per_sample, sep="\t", index=False)
