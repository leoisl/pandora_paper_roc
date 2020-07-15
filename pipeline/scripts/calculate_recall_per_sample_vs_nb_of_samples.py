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


# setup
all_recall_reports_for_one_sample_until_nb_samples = (
    snakemake.input.all_recall_reports_for_one_sample_until_nb_samples
)
sample = snakemake.wildcards.sample_id
tool = snakemake.wildcards.tool
coverage = snakemake.wildcards.coverage
coverage_threshold = snakemake.wildcards.coverage_threshold
strand_bias_threshold = snakemake.wildcards.strand_bias_threshold
gaps_threshold = snakemake.wildcards.gaps_threshold
list_with_number_of_samples = snakemake.params.list_with_number_of_samples
recall_file_for_one_sample_until_nb_samples_filename = Path(snakemake.output.recall_file_for_one_sample_until_nb_samples)


# API usage
logging.info(f"Loading report")
recall_report = RecallReport.from_files(all_recall_reports_for_one_sample_until_nb_samples,
                                        concatenate_dfs_one_by_one_keeping_only_best_mappings=True)

logging.info(f"Creating calculator")
recall_calculator = RecallCalculator(recall_report)

logging.info(f"Calculating recall")

recall_file_for_one_sample_until_nb_samples = recall_calculator.get_recall_vs_nb_of_samples_report(list_with_number_of_samples)


# output
recall_file_for_one_sample_until_nb_samples.to_csv(recall_file_for_one_sample_until_nb_samples_filename, index=False)