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

logging.info(f"Creating masker from {snakemake.input.mask}")
with open(snakemake.input.mask) as bed:
    mask = RecallMasker.from_bed(bed)

with pysam.AlignmentFile(snakemake.input.sam) as sam:
    logging.info(f"Masking SAM records")
    records = mask.filter_records(sam)

logging.info("Creating classifier")
classifier = RecallClassifier(sam=records, name=snakemake.wildcards.sample_id)
logging.info("Creating reporter")
reporter = RecallReporter(classifiers=[classifier])

with open(snakemake.output.report, "w") as output:
    logging.info("Generating and saving report")
    reporter.save(output)
