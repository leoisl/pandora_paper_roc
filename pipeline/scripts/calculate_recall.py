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
gt_conf_percentiles = snakemake.params.gt_conf_percentiles
tool = snakemake.wildcards.tool
coverage = snakemake.wildcards.coverage
coverage_threshold = snakemake.wildcards.coverage_threshold
strand_bias_threshold = snakemake.wildcards.strand_bias_threshold
gaps_threshold = snakemake.wildcards.gaps_threshold

recall_file_for_all_samples_and_all_gt_conf_percentile = Path(snakemake.output.recall_file_for_all_samples_and_all_gt_conf_percentile)
recall_final_report = snakemake.output.recall_final_report


# API usage
logging.info(f"Loading report")
recall_report = RecallReport.from_files(recall_report_files_for_all_samples_and_all_gt_conf_percentile)

logging.info(f"Creating calculator")
recall_calculator = RecallCalculator(recall_report)

logging.info(f"Calculating recall")
recall_df = recall_calculator.get_recall_report(gt_conf_percentiles)

metadata_df = pd.DataFrame(
    data={
        "tool": [tool] * len(recall_df),
        "coverage": [coverage] * len(recall_df),
        "coverage_threshold": [coverage_threshold] * len(recall_df),
        "strand_bias_threshold": [strand_bias_threshold] * len(recall_df),
        "gaps_threshold": [gaps_threshold] * len(recall_df),
    }
)
output_df = pd.concat([recall_df, metadata_df], axis=1)


# output
logging.info(f"Outputting recall file")
output_df.to_csv(recall_file_for_all_samples_and_all_gt_conf_percentile, sep="\t")
with open(recall_final_report, "w") as recall_report_filehandler:
    recall_report.save_report(recall_report_filehandler)

logging.info(f"Done")
