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


# setup
all_recall_reports_with_no_gt_conf_filter = snakemake.input.all_recall_reports_with_no_gt_conf_filter
recalls_per_number_of_samples = snakemake.output.recalls_per_number_of_samples

for recall_per_number_of_samples in recalls_per_number_of_samples:
    Path(recall_per_number_of_samples).write_text("asd")
