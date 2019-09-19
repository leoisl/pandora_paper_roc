import pandas as pd

precision_df = pd.read_csv(snakemake.input.precision_file, sep="\t")
recall_df = pd.read_csv(snakemake.input.recall_file, sep="\t")

error_rate_and_recall_df = pd.merge(precision_df, recall_df, on="GT")[
    ["GT", "error_rate", "recall"]
]
error_rate_and_recall_df.to_csv(snakemake.output.error_rate_and_recall_file, sep="\t")
