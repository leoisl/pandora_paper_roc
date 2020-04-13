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
from evaluate.report import RecallReport


# setup
recall_report_files_for_one_sample_and_all_gt_conf_percentiles = (
    snakemake.input.recall_report_files_for_one_sample_and_all_gt_conf_percentiles
)
recall_report_per_sample_for_calculator = snakemake.output.recall_report_per_sample_for_calculator


# API usage
logging.info(f"Loading report")
recall_report = RecallReport.from_files(recall_report_files_for_one_sample_and_all_gt_conf_percentiles,
                                        concatenate_dfs_one_by_one_keeping_only_best_mappings=True)

with open(recall_report_per_sample_for_calculator, "w") as recall_report_per_sample_for_calculator_filehandler:
    recall_report.save_report(recall_report_per_sample_for_calculator_filehandler)

logging.info(f"Done")
