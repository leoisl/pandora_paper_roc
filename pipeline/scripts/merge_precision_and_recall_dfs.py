import pandas as pd

# setup
precision_file = snakemake.input.precision_file
recall_file = snakemake.input.recall_file
output_file = snakemake.output.error_rate_and_recall_file

# read
precision_df = pd.read_csv(precision_file, sep="\t")
recall_df = pd.read_csv(recall_file, sep="\t")

# merge
error_rate_and_recall_df = pd.merge(precision_df, recall_df, on=["GT", "label"])[
    ["GT", "label", "error_rate", "recall"]
]

# output
error_rate_and_recall_df.to_csv(output_file, sep="\t")
