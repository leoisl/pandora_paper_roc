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
recall_report_files_for_all_samples_and_all_gt_conf_percentile = (
    snakemake.input.recall_report_files_for_all_samples_and_all_gt_conf_percentile
)
output = Path(snakemake.output.recall_file_for_all_samples_and_all_gt_conf_percentile)

number_of_points_in_ROC_curve = int(snakemake.params.number_of_points_in_ROC_curve)

tool = snakemake.wildcards.tool
coverage = snakemake.wildcards.coverage
coverage_threshold = snakemake.wildcards.coverage_threshold
strand_bias_threshold = snakemake.wildcards.strand_bias_threshold
gaps_threshold = snakemake.wildcards.gaps_threshold


# API usage
logging.info(f"Loading report")
recall_report = RecallReport.from_files(recall_report_files_for_all_samples_and_all_gt_conf_percentile)

logging.info(f"Creating calculator")
recall_calculator = RecallCalculator(recall_report)

logging.info(f"Calculating recall")
recall_df = recall_calculator.get_recall_report(number_of_points_in_ROC_curve)

metadata_df = pd.DataFrame(
    data={
        "tool": [tool] * len(recall_df),
        "coverage": [coverage] * len(recall_df),
        "coverage_threshold": [coverage_threshold] * len(recall_df),
        "strand_bias_threshold": [strand_bias_threshold] * len(recall_df),
        "gaps_threshold": [gaps_threshold] * len(recall_df),
    }
)
output_df = pd.concat([recall_df, recall_df], axis=1)


# output
logging.info(f"Outputting recall file")
output_df.to_csv(output, sep="\t")

logging.info(f"Done")
