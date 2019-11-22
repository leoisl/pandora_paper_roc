import pandas as pd


def get_concatenated_df(files):
    dfs = [pd.read_csv(file, sep="\t") for file in files]
    concatenated_df = pd.concat(dfs, ignore_index=True)
    return concatenated_df


all_precision_files = snakemake.input.all_precision_files
all_recall_files = snakemake.input.all_recall_files

precision_df = get_concatenated_df(all_precision_files)
recall_df = get_concatenated_df(all_recall_files)

error_rate_and_recall_df = pd.merge(precision_df, recall_df, on=["GT", "label"])[
    ["GT", "label", "error_rate", "recall"]
]
error_rate_and_recall_df.to_csv(snakemake.output.error_rate_and_recall_file, sep="\t")
