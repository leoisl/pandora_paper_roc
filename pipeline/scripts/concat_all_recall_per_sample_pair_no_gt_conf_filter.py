from utils import get_concatenated_df

# setup
all_recall_per_sample_pair_no_gt_conf_filter = snakemake.input.all_recall_per_sample_pair_no_gt_conf_filter
recall_per_sample_pair = snakemake.output.recall_per_sample_pair

# read
all_dfs = get_concatenated_df(all_recall_per_sample_pair_no_gt_conf_filter, separator="\t")

# output
all_dfs.to_csv(recall_per_sample_pair, sep="\t", index=False)
