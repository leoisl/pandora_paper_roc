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
from evaluate.calculator import RecallCalculator


# setup
all_recall_reports_with_no_gt_conf_filter = snakemake.input.all_recall_reports_with_no_gt_conf_filter
recall_per_number_of_samples_filename = snakemake.output.recall_per_number_of_samples
list_with_number_of_samples = snakemake.params.list_with_number_of_samples


# API usage
logging.info(f"Loading report")
recall_report = RecallReport.from_files(all_recall_reports_with_no_gt_conf_filter,
                                        concatenate_dfs_one_by_one_keeping_only_best_mappings=True)

logging.info(f"Creating calculator")
recall_calculator = RecallCalculator(recall_report)

logging.info(f"Calculating recall")
recall_per_number_of_samples = recall_calculator.get_recall_alleles_vs_nb_of_samples_report(list_with_number_of_samples)


# output
recall_per_number_of_samples.to_csv(recall_per_number_of_samples_filename, index=False)