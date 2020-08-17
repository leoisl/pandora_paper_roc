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
from evaluate.calculator import RecallCalculator
from evaluate.report import RecallReport
import pandas as pd


# setup
all_recall_reports_for_one_sample = snakemake.input.all_recall_reports_for_one_sample
sample = snakemake.wildcards.sample_id
tool = snakemake.wildcards.tool
coverage = snakemake.wildcards.coverage
coverage_threshold = snakemake.wildcards.coverage_threshold
strand_bias_threshold = snakemake.wildcards.strand_bias_threshold
gaps_threshold = snakemake.wildcards.gaps_threshold
list_with_number_of_samples = snakemake.params.list_with_number_of_samples
recall_file_for_one_sample_vs_nb_samples_filename = Path(snakemake.output.recall_file_for_one_sample_vs_nb_samples)


# API usage
logging.info(f"Loading report")
recall_report = RecallReport.from_files(all_recall_reports_for_one_sample,
                                        concatenate_dfs_one_by_one_keeping_only_best_mappings=True)

logging.info(f"Creating calculator")
recall_calculator = RecallCalculator(recall_report)

logging.info(f"Calculating recall")
recall_df = recall_calculator.get_recall_report_wrt_truth_probes_for_those_present_in_a_given_nb_of_samples(list_with_number_of_samples)

metadata_df = pd.DataFrame(
    data={
        "tool": [tool] * len(recall_df),
        "coverage": [coverage] * len(recall_df),
        "coverage_threshold": [coverage_threshold] * len(recall_df),
        "strand_bias_threshold": [strand_bias_threshold] * len(recall_df),
        "gaps_threshold": [gaps_threshold] * len(recall_df),
        "sample": [sample] * len(recall_df)
    }
)
output_df = pd.concat([recall_df, metadata_df], axis=1)


# output
logging.info(f"Outputting recall file")
output_df.to_csv(recall_file_for_one_sample_vs_nb_samples_filename, index=False)
logging.info(f"Done")
