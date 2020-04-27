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
from evaluate.deduplicate_pairwise_snps import DeduplicationGraph, ConsistentPangenomeVariations
import pickle


# setup
snps_dfs_filepaths = snakemake.input.snps_dfs
deduplicated_snps_dfs_filepaths = snakemake.output.deduplicated_snps_dfs
deduplicated_snps_dfs_text_filepaths = snakemake.output.deduplicated_snps_dfs_text



# API usage
# create the deduplication graph
deduplication_graph = DeduplicationGraph(number_of_positions_in_each_index_bucket=100)
logging.info("Building nodes...")
for snps_df in snps_dfs_filepaths:
    deduplication_graph.add_variants_from_ShowSNPsDataframe_filepath(snps_df)

logging.info("Building edges...")
deduplication_graph.build_edges()

# create the consistent pangenome variations
logging.info("Getting pangenome variations...")
pangenome_variations = deduplication_graph.get_pangenome_variations()
logging.info("Getting consistent pangenome variations...")
consistent_pangenome_variations = ConsistentPangenomeVariations(pangenome_variations)

# write the enriched and filtered deduplicated variations
logging.info("Outputting...")
for snps_df_filepath, deduplicated_snps_df_filepath, deduplicated_snps_df_text_filepath \
        in zip(snps_dfs_filepaths, deduplicated_snps_dfs_filepaths, deduplicated_snps_dfs_text_filepaths):
    deduplicated_snps_df = consistent_pangenome_variations.build_DeduplicatedVariationsDataframe_from_ShowSNPsDataframe(snps_df_filepath)
    with open(deduplicated_snps_df_filepath, "wb") as deduplicated_snps_df_fh:
        pickle.dump(deduplicated_snps_df, file=deduplicated_snps_df_fh)
    deduplicated_snps_df.to_csv(deduplicated_snps_df_text_filepath, index=False)