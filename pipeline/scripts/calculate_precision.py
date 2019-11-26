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



from evaluate.calculator import PrecisionCalculator, EmptyReportError
import pandas as pd
import numpy as np


# setup
precision_report_files_for_all_samples = (
    snakemake.input.precision_report_files_for_all_samples
)
output = Path(snakemake.output.precision_file_for_all_samples)

min_gt = float(snakemake.wildcards.min_gt)
step_gt = float(snakemake.wildcards.step_gt)
max_gt = float(snakemake.wildcards.max_gt)

tool = snakemake.wildcards.tool
coverage = snakemake.wildcards.coverage
coverage_threshold = snakemake.wildcards.coverage_threshold
strand_bias_threshold = snakemake.wildcards.strand_bias_threshold
gaps_threshold = snakemake.wildcards.gaps_threshold


# API usage
logging.info(f"Creating calculator")
precision_calculator = PrecisionCalculator.from_files(
    precision_report_files_for_all_samples
)


max_gt = min(max_gt, precision_calculator.get_maximum_gt_conf())
logging.info(
    f"Calculating precision with min_gt = {min_gt}, step_gt = {step_gt}, and max_gt = {max_gt}"
)

gts = []
precisions = []
error_rates = []
all_gts = list(np.arange(min_gt, max_gt, step_gt)) + [max_gt]
for gt in all_gts:
    try:
        precision = precision_calculator.calculate_precision(gt)
        gts.append(gt)
        precisions.append(precision)
        error_rates.append(1 - precision)
    except EmptyReportError:
        pass


precision_df = pd.DataFrame(
    data={
        "tool": [tool] * len(gts),
        "coverage": [coverage] * len(gts),
        "coverage_threshold": [coverage_threshold] * len(gts),
        "strand_bias_threshold": [strand_bias_threshold] * len(gts),
        "gaps_threshold": [gaps_threshold] * len(gts),
        "GT": gts,
        "precision": precisions,
        "error_rate": error_rates
    }
)


# output
logging.info(f"Outputting precision file")
precision_df.to_csv(output, sep="\t")


logging.info(f"Done")