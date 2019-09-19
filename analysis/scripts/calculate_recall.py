from pathlib import Path
import sys

sys.path.append(str(Path().absolute()))
import logging

log_level = "DEBUG"
logging.basicConfig(
    filename=str(snakemake.log),
    filemode="w",
    level=log_level,
    format="[%(asctime)s]:%(levelname)s: %(message)s",
    datefmt="%d/%m/%Y %I:%M:%S %p",
)

from evaluate.calculator import RecallCalculator, EmptyReportError
import pandas as pd
import numpy as np

all_recall_report_files = snakemake.input.all_recall_report_files
recall_calculator = RecallCalculator.from_files(all_recall_report_files)

output = Path(snakemake.output.recall_file)

min_gt = float(snakemake.wildcards.min_gt)
step_gt = float(snakemake.wildcards.step_gt)
max_gt = float(snakemake.wildcards.max_gt)
max_gt = min(max_gt, recall_calculator.get_maximum_gt_conf())

logging.info(
    f"Generating recall file with min_gt = {min_gt}, step_gt = {step_gt}, and max_gt = {max_gt}"
)

gts = []
recalls = []
all_gts = list(np.arange(min_gt, max_gt, step_gt)) + [max_gt]
for gt in all_gts:
    try:
        recall = recall_calculator.calculate_recall(gt)
        gts.append(gt)
        recalls.append(recall)
    except EmptyReportError:
        pass

recall_df = pd.DataFrame(data={"GT": gts, "recall": recalls})
recall_df.to_csv(output, sep="\t")
