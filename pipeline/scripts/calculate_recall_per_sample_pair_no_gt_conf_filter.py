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



from evaluate.calculator import RecallCalculator, EmptyReportError
from evaluate.report import RecallReport
import pandas as pd


# setup
all_recall_reports_for_one_sample_pair_with_no_gt_conf_filter = (
    snakemake.input.all_recall_reports_for_one_sample_pair_with_no_gt_conf_filter
)
sample = snakemake.wildcards.sample_id
sample_pair = snakemake.wildcards.sample_pair
tool = snakemake.wildcards.tool
coverage = snakemake.wildcards.coverage
coverage_threshold = snakemake.wildcards.coverage_threshold
strand_bias_threshold = snakemake.wildcards.strand_bias_threshold
gaps_threshold = snakemake.wildcards.gaps_threshold
aligned_bases_percentage_file = Path(snakemake.input.aligned_bases_percentage)
aligned_bases_percentage = float(aligned_bases_percentage_file.read_text())
recall_file_for_one_sample_pair_with_no_gt_conf_filter = Path(snakemake.output.recall_file_for_one_sample_pair_with_no_gt_conf_filter)


# API usage
logging.info(f"Loading report")
recall_report = RecallReport.from_files([all_recall_reports_for_one_sample_pair_with_no_gt_conf_filter])

logging.info(f"Creating calculator")
recall_calculator = RecallCalculator(recall_report)

logging.info(f"Calculating recall")
recall_df = recall_calculator.get_recall_report([0])

metadata_df = pd.DataFrame(
    data={
        "tool": [tool] * len(recall_df),
        "coverage": [coverage] * len(recall_df),
        "coverage_threshold": [coverage_threshold] * len(recall_df),
        "strand_bias_threshold": [strand_bias_threshold] * len(recall_df),
        "gaps_threshold": [gaps_threshold] * len(recall_df),
        "sample": [sample] * len(recall_df),
        "sample_pair": [sample_pair] * len(recall_df),
        "aligned_bases_percentage": [aligned_bases_percentage] * len(recall_df)
    }
)
output_df = pd.concat([recall_df, metadata_df], axis=1)


# output
logging.info(f"Outputting recall file")
output_df.to_csv(recall_file_for_one_sample_pair_with_no_gt_conf_filter, sep="\t", index=False)
logging.info(f"Done")
