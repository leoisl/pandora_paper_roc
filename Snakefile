from pathlib import Path
import itertools
from snakemake.utils import min_version, validate
import pandas as pd
from collections import defaultdict
from pipeline.scripts.utils import get_sample_pairs_containing_given_sample

min_version("5.1.0")

# ======================================================
# Config files
# ======================================================
configfile: "config.yaml"

validate(config, "pipeline/schemas/config.schema.yaml")
samples = pd.read_csv(config["samples"])
validate(samples, "pipeline/schemas/samples.schema.yaml")
samples.rename(columns={"reference": "reference_assembly"}, inplace=True)
samples["reference_assembly"] += ".chrom_only"

variant_calls = pd.read_csv(config["variant_calls"])
validate(variant_calls, "pipeline/schemas/variant_calls.schema.yaml")
variant_calls.rename(columns={"reference": "vcf_reference"}, inplace=True)
variant_calls["vcf"] += ".~~sample~~names~~fixed.vcf"

max_gt_conf_percentile = int(config['max_gt_conf_percentile'])

# ======================================================
# Global variables
# ======================================================
output_folder = config['output_folder']
data: pd.DataFrame = pd.merge(variant_calls, samples, on="sample_id")
data = data.set_index(["sample_id", "coverage", "tool"], drop=False)
samples = samples.set_index(["sample_id"], drop=False)
sample_pairs = [(sample1, sample2) for sample1, sample2 in itertools.combinations(sorted(samples["sample_id"]), r=2)]
gt_conf_percentiles = list(range(0, max_gt_conf_percentile))


# ======================================================
# Helper functions
# ======================================================
def get_coverage_filters(tool):
    return [str(elem) for elem in config['coverage_filters']]

def get_strand_bias_filters(tool):
    if tool.startswith("pandora"):
        return [str(elem) for elem in config['strand_bias_filters']]
    else:
        return ["Not_App"]

def get_gaps_filters(tool):
    if tool.startswith("pandora"):
        return [str(elem) for elem in config['gaps_filters']]
    else:
        return ["Not_App"]


# ======================================================
# Pipeline files
# ======================================================
files = []

# Precision files
all_precision_files=[]

cov_tool_and_filters_to_precision_report_files = defaultdict(list)
all_nb_of_records_removed_with_mapq_sam_records_filter_files_for_precision = []
for index, row in data.iterrows():
    sample_id, coverage, tool = row["sample_id"], row["coverage"], row["tool"]
    for coverage_threshold in get_coverage_filters(tool):
        for strand_bias_threshold in get_strand_bias_filters(tool):
            for gaps_threshold in get_gaps_filters(tool):
                report_file = f"{output_folder}/precision/reports_from_probe_mappings/{sample_id}/{coverage}/{tool}/coverage_filter_{coverage_threshold}/strand_bias_filter_{strand_bias_threshold}/gaps_filter_{gaps_threshold}/variant_calls_probeset_report.tsv"
                cov_tool_and_filters_to_precision_report_files[(coverage, tool, coverage_threshold, strand_bias_threshold, gaps_threshold)].append(report_file)

                nb_of_records_removed_with_mapq_sam_records_filter_file = f"{output_folder}/precision/reports_from_probe_mappings/{sample_id}/{coverage}/{tool}/coverage_filter_{coverage_threshold}/strand_bias_filter_{strand_bias_threshold}/gaps_filter_{gaps_threshold}/nb_of_records_removed_with_mapq_sam_records_filter.csv"
                all_nb_of_records_removed_with_mapq_sam_records_filter_files_for_precision.append(nb_of_records_removed_with_mapq_sam_records_filter_file)

for coverage, tool, coverage_threshold, strand_bias_threshold, gaps_threshold in cov_tool_and_filters_to_precision_report_files:
    all_precision_files.append(f"{output_folder}/precision/precision_files/{coverage}/{tool}/coverage_filter_{coverage_threshold}/strand_bias_filter_{strand_bias_threshold}/gaps_filter_{gaps_threshold}/precision.tsv")

files.extend(all_precision_files)



# Recall files
all_recall_files=[]

for index, row in data.iterrows():
    sample_id, coverage, tool = row["sample_id"], row["coverage"], row["tool"]
    for filename_prefix in get_sample_pairs_containing_given_sample(sample_pairs, sample_id):
        files_with_filters = expand(f"{output_folder}/recall/map_probes/{sample_id}/{coverage}/{tool}/coverage_filter_{{coverage_threshold}}/strand_bias_filter_{{strand_bias_threshold}}/gaps_filter_{{gaps_threshold}}/gt_conf_percentile_{{gt_conf_percentile}}/{filename_prefix}.sam", coverage_threshold = get_coverage_filters(tool), strand_bias_threshold = get_strand_bias_filters(tool), gaps_threshold = get_gaps_filters(tool), gt_conf_percentile=gt_conf_percentiles)

cov_tool_and_filters_to_recall_report_files = defaultdict(list)
all_recall_per_sample_no_gt_conf_filter = set()
all_recall_per_sample_pair_no_gt_conf_filter = set()
for index, row in data.iterrows():
    sample_id, coverage, tool = row["sample_id"], row["coverage"], row["tool"]
    for filename_prefix in get_sample_pairs_containing_given_sample(sample_pairs, sample_id):
        for coverage_threshold, strand_bias_threshold, gaps_threshold, gt_conf_percentile in \
            itertools.product(get_coverage_filters(tool), get_strand_bias_filters(tool), get_gaps_filters(tool), gt_conf_percentiles):
            report_file = f"{output_folder}/recall/reports/{sample_id}/{coverage}/{tool}/coverage_filter_{coverage_threshold}/strand_bias_filter_{strand_bias_threshold}/gaps_filter_{gaps_threshold}/gt_conf_percentile_{gt_conf_percentile}/{filename_prefix}.report.tsv"
            cov_tool_and_filters_to_recall_report_files[(coverage, tool, str(coverage_threshold), str(strand_bias_threshold), str(gaps_threshold))].append(report_file)
            all_recall_per_sample_no_gt_conf_filter.add(f"{output_folder}/recall/recall_files_per_sample/{sample_id}/{coverage}/{tool}/coverage_filter_{coverage_threshold}/strand_bias_filter_{strand_bias_threshold}/gaps_filter_{gaps_threshold}/recall.tsv")
            all_recall_per_sample_pair_no_gt_conf_filter.add(f"{output_folder}/recall/recall_files_per_sample_pair/{sample_id}/{coverage}/{tool}/coverage_filter_{coverage_threshold}/strand_bias_filter_{strand_bias_threshold}/gaps_filter_{gaps_threshold}/{filename_prefix}.recall.tsv")
all_recall_per_sample_no_gt_conf_filter = list(all_recall_per_sample_no_gt_conf_filter)
all_recall_per_sample_pair_no_gt_conf_filter = list(all_recall_per_sample_pair_no_gt_conf_filter)


for coverage, tool, coverage_threshold, strand_bias_threshold, gaps_threshold in cov_tool_and_filters_to_recall_report_files:
    all_recall_files.append(f"{output_folder}/recall/recall_files/{coverage}/{tool}/coverage_filter_{coverage_threshold}/strand_bias_filter_{strand_bias_threshold}/gaps_filter_{gaps_threshold}/recall.tsv")

files.extend(all_recall_files)



# Plot files
all_plot_data_intermediate_files = []
for coverage, tool, coverage_threshold, strand_bias_threshold, gaps_threshold in cov_tool_and_filters_to_recall_report_files:
    all_plot_data_intermediate_files.append(f"{output_folder}/plot_data/{coverage}/{tool}/coverage_filter_{coverage_threshold}/strand_bias_filter_{strand_bias_threshold}/gaps_filter_{gaps_threshold}/ROC_data.tsv")
final_plot_data_file = f"{output_folder}/plot_data/ROC_data.tsv"
final_all_nb_of_records_removed_with_mapq_sam_records_filter_file = f"{output_folder}/plot_data/nb_of_records_removed_with_mapq_sam_records_filter_for_precision.csv"
recall_per_sample_file = f"{output_folder}/plot_data/recall_per_sample.tsv"
recall_per_sample_pair_file = f"{output_folder}/plot_data/recall_per_sample_pair.tsv"

files.extend(all_plot_data_intermediate_files)
files.append(final_plot_data_file)
files.append(final_all_nb_of_records_removed_with_mapq_sam_records_filter_file)
files.append(recall_per_sample_file)
files.append(recall_per_sample_pair_file)



# ======================================================
# Rules
# ======================================================
rule all:
    input: files

rules_dir = Path("pipeline/rules/")
include: str(rules_dir / "common.smk")
include: str(rules_dir / "recall.smk")
include: str(rules_dir / "precision.smk")
include: str(rules_dir / "plot.smk")

localrules: fix_snippy_sample_name, make_empty_depth_file