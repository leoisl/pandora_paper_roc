from pathlib import Path
import sys
sys.path.append(str(Path().absolute()))
import logging
log_level = "INFO"
logging.basicConfig(
    filename=str(snakemake.log),
    filemode="w",
    level=log_level,
    format="[%(asctime)s]:%(levelname)s: %(message)s",
    datefmt="%d/%m/%Y %I:%M:%S %p",
)
import pickle
import re

# inputs
snps_df_filenames = snakemake.input.snps_dfs
deduplicated_snps_df_filenames = snakemake.input.deduplicated_snps_dfs

logging.info("Reading show-snps dataframes")
for snps_df_filename, deduplicated_snps_df_filename in zip(snps_df_filenames, deduplicated_snps_df_filenames):
    with open(snps_df_filename, "rb") as snps_df_fh, open(deduplicated_snps_df_filename, "wb") as deduplicated_snps_df_fh:
