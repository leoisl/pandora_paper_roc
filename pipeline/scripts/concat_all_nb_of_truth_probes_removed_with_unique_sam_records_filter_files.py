from utils import get_concatenated_df

# setup
all_nb_of_truth_probes_removed_with_unique_sam_records_filter_files_for_recall = snakemake.input.all_nb_of_truth_probes_removed_with_unique_sam_records_filter_files_for_recall
nb_of_truth_probes_removed_with_unique_sam_records_filter_for_recall_filepath = snakemake.output.nb_of_truth_probes_removed_with_unique_sam_records_filter_for_recall_filepath

# read
all_dfs = get_concatenated_df(all_nb_of_truth_probes_removed_with_unique_sam_records_filter_files_for_recall,separator=",",
                              fields_to_keep=["sample_id", "sample_pair", "nb_of_truth_probes_before_unique_sam_records_filter",
                                              "nb_of_truth_probes_after_unique_sam_records_filter",
                                              "nb_of_truth_probes_removed_with_unique_sam_records_filter",
                                              "nb_of_truth_probes_removed_with_unique_sam_records_filter_proportion"])

# output
all_dfs.to_csv(nb_of_truth_probes_removed_with_unique_sam_records_filter_for_recall_filepath)
