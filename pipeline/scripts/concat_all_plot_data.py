from utils import get_concatenated_df

# setup
all_plot_data_intermediate_files = snakemake.input.all_plot_data_intermediate_files
output_file = snakemake.output.final_plot_data_file

# read
all_dfs = get_concatenated_df(all_plot_data_intermediate_files,separator="\t",
                              fields_to_keep=["tool", "coverage", "coverage_threshold", "strand_bias_threshold",
                                              "gaps_threshold", "step_GT",
                                              "error_rate", "nb_of_correct_calls", "nb_of_total_calls",
                                              "recalls_wrt_truth_probes", "nbs_of_truth_probes_found", "nbs_of_truth_probes_in_total",
                                              "recalls_wrt_variants_where_all_allele_seqs_were_found", "recalls_wrt_variants_found_wrt_alleles",
                                              "nbs_variants_where_all_allele_seqs_were_found", "nbs_variants_found_wrt_alleles", "nbs_variants_total"])


# TODO: check this
# some sanity checks
# all_nbs_of_truth_probes_in_total = all_dfs.nbs_of_truth_probes_in_total.unique()
# assert len(all_nbs_of_truth_probes_in_total) == 1

# output
all_dfs.to_csv(output_file, sep="\t")
