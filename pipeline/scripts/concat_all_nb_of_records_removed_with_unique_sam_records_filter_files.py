from utils import get_concatenated_df

# setup
all_nb_of_records_removed_with_unique_sam_records_filter_files = snakemake.input
output_file = str(snakemake.output)

# read
all_dfs = get_concatenated_df(all_nb_of_records_removed_with_unique_sam_records_filter_files,separator=",",
                              fields_to_keep=["tool", "nb_of_records_before_mapq_sam_records_filter",
                                              "nb_of_records_after_mapq_sam_records_filter",
                                              "nb_of_records_removed_with_mapq_sam_records_filter",
                                              "nb_of_records_removed_with_mapq_sam_records_filter_proportion"])

# output
all_dfs.to_csv(output_file)
