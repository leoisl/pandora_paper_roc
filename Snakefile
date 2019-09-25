from pathlib import Path
import itertools
from snakemake.utils import min_version, validate
import pandas as pd
from collections import defaultdict

min_version("5.1.0")

# ======================================================
# Config files
# ======================================================
configfile: "config.yaml"

validate(config, "analysis/schemas/config.schema.yaml")
samples = pd.read_csv(config["samples"])
validate(samples, "analysis/schemas/samples.schema.yaml")
samples.rename(columns={"reference": "reference_assembly"}, inplace=True)

variant_calls = pd.read_csv(config["variant_calls"])
validate(variant_calls, "analysis/schemas/variant_calls.schema.yaml")
variant_calls.rename(columns={"reference": "vcf_reference"}, inplace=True)


# ======================================================
# Global variables
# ======================================================
data: pd.DataFrame = pd.merge(variant_calls, samples, on="sample_id")
data = data.set_index(["sample_id", "coverage", "tool"], drop=False)
samples = samples.set_index(["sample_id"], drop=False)
sample_pairs = [(sample1, sample2) for sample1, sample2 in itertools.combinations(sorted(samples["sample_id"]), r=2)]
min_gt = config['genotype_confidence_min']
step_gt = config['genotype_confidence_step']
max_gt = config['genotype_confidence_max']
coverage_filters = config['coverage_filters']
strand_bias_filters = config['strand_bias_filters']
gaps_filters = config['gaps_filters']

files = []

# Common files
for index, row in data.iterrows():
    sample_id, coverage, tool = row["sample_id"], row["coverage"], row["tool"]
    files_with_filters = expand(f"analysis/variant_calls_probesets/{sample_id}/{coverage}/{tool}/coverage_filter_{{coverage_threshold}}/strand_bias_filter_{{strand_bias_threshold}}/gaps_filter_{{gaps_threshold}}/variant_calls_probeset.fa", coverage_threshold = coverage_filters, strand_bias_threshold = strand_bias_filters, gaps_threshold = gaps_filters)
    files.extend(files_with_filters)

all_precision_files=[]
all_recall_files=[]



# Precision files
for index, row in data.iterrows():
    sample_id, coverage, tool = row["sample_id"], row["coverage"], row["tool"]
    files.extend(
        [
            f"analysis/precision/variant_calls_probesets_mapped_to_refs/{sample_id}/{coverage}/{tool}/variant_calls_probeset_mapped.sam",
        ]
    )

tool_and_coverage_to_precision_report_files = defaultdict(list)
for index, row in data.iterrows():
    sample_id, coverage, tool = row["sample_id"], row["coverage"], row["tool"]
    tool_and_coverage_to_precision_report_files[f"{tool}_{coverage}"].append(f"analysis/precision/reports_from_probe_mappings/{sample_id}/{coverage}/{tool}/variant_calls_probeset_report.tsv")
for precision_report_files in tool_and_coverage_to_precision_report_files.values():
    files.extend(precision_report_files)

all_precision_files = [f"analysis/precision/precision_{tool_and_coverage}_gt_min_{min_gt}_step_{step_gt}_max_{max_gt}.tsv"
                       for tool_and_coverage in tool_and_coverage_to_precision_report_files]




# Recall files
for sample1, sample2 in sample_pairs:
    files.extend(
        [
            f"analysis/recall/truth_probesets/{sample1}/{sample1}_and_{sample2}.truth_probeset.fa",
            f"analysis/recall/truth_probesets/{sample2}/{sample1}_and_{sample2}.truth_probeset.fa",
        ]
    )

tool_and_coverage_to_recall_report_files = defaultdict(list)
for index, row in data.iterrows():
    sample_id, coverage, tool = row["sample_id"], row["coverage"], row["tool"]
    for sample1, sample2 in [pair for pair in sample_pairs if sample_id in pair]:
        filename_prefix = f"{sample1}_and_{sample2}"
        tool_and_coverage_to_recall_report_files[f"{tool}_{coverage}"].append(f"analysis/recall/reports/{sample_id}/{coverage}/{tool}/{filename_prefix}.report.tsv")
for recall_report_files in tool_and_coverage_to_recall_report_files.values():
    files.extend(recall_report_files)

all_recall_files = [f"analysis/recall/recall_{tool_and_coverage}_gt_min_{min_gt}_step_{step_gt}_max_{max_gt}.tsv"
                       for tool_and_coverage in tool_and_coverage_to_precision_report_files]




# Plot files
files.append(f"analysis/plot/error_rate_and_recall_gt_min_{config['genotype_confidence_min']}_step_{config['genotype_confidence_step']}_max_{config['genotype_confidence_max']}.tsv")
files.append(f"analysis/plot/error_rate_and_recall_gt_min_{config['genotype_confidence_min']}_step_{config['genotype_confidence_step']}_max_{config['genotype_confidence_max']}.pdf")




# ======================================================
# Rules
# ======================================================
rule all:
    input: files

rules_dir = Path("analysis/rules/")
include: str(rules_dir / "common.smk")
include: str(rules_dir / "recall.smk")
include: str(rules_dir / "precision.smk")
include: str(rules_dir / "plot.smk")
