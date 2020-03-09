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



import pysam
from evaluate.masker import RecallMasker
from evaluate.classifier import RecallClassifier
from evaluate.reporter import RecallReporter


# setup
mask_filepath = snakemake.input.mask
sams_filepath = snakemake.input.sams
sample_id = snakemake.wildcards.sample_id
variant_call_recall_reports = snakemake.output.reports
gt_conf_percentiles = snakemake.params.gt_conf_percentiles


# API usage
logging.info(f"Creating masker from {mask_filepath}")
with open(mask_filepath) as bed:
    masker = RecallMasker.from_bed(bed)

for sam_filepath, variant_call_recall_report, gt_conf_percentile in zip(sams_filepath, variant_call_recall_reports, gt_conf_percentiles):
    logging.info(f"Masking SAM records")
    with pysam.AlignmentFile(sam_filepath) as sam:
        records = masker.filter_records(sam)

    logging.info("Creating classifier")
    classifier = RecallClassifier(sam=records, name=sample_id)

    logging.info("Creating reporter")
    reporter = RecallReporter(classifiers=[classifier])

    logging.info("Generating report")

    # TODO: we are passing gt_conf_percentile (values in [0, 100, 1]) as gt_conf
    # TODO: fix this? It does not really matter as we use step gt (which is gt_conf_percentile) anyway later
    report = reporter.generate_report(gt_conf_percentile)


    # output
    logging.info("Saving report")
    with open(variant_call_recall_report, "w") as output:
        reporter.save_report(report, output)

logging.info("Done")