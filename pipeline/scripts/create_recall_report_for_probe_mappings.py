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
sam_filepath = snakemake.input.sam
sample_id = snakemake.wildcards.sample_id
output = snakemake.output.report


# API usage
logging.info(f"Creating masker from {mask_filepath}")
with open(mask_filepath) as bed:
    masker = RecallMasker.from_bed(bed)

logging.info(f"Masking SAM records")
with pysam.AlignmentFile(sam_filepath) as sam:
    records = masker.filter_records(sam)

logging.info("Creating classifier")
classifier = RecallClassifier(sam=records, name=sample_id)

logging.info("Creating reporter")
reporter = RecallReporter(classifiers=[classifier])


# output
logging.info("Generating and saving report")
with open(output, "w") as output:
    reporter.save(output)

logging.info("Done")