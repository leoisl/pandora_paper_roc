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

recall_report_files_for_tool_and_coverage = (
    snakemake.input.recall_report_files_for_tool_and_coverage
)
recall_calculator = RecallCalculator.from_files(
    recall_report_files_for_tool_and_coverage
)

output = Path(snakemake.output.recall_file_for_tool_and_coverage)

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

labels = [snakemake.wildcards.tool_and_coverage] * len(gts)
recall_df = pd.DataFrame(data={"GT": gts, "recall": recalls, "label": labels})
recall_df.to_csv(output, sep="\t")
