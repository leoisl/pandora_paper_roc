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
from evaluate.classifier import PrecisionClassifier
from evaluate.masker import PrecisionMasker
from evaluate.reporter import PrecisionReporter


# setup
sam_filepath = snakemake.input.variant_call_probeset_mapped_to_ref
sample_id = snakemake.wildcards.sample_id
mask_filepath = snakemake.input.mask
variant_call_precision_report = snakemake.output.variant_call_precision_report


# API usage
logging.info(f"Creating masker from {mask_filepath}")
with open(mask_filepath) as bed:
    masker = PrecisionMasker.from_bed(bed)

logging.info(f"Masking SAM records")
with pysam.AlignmentFile(sam_filepath) as sam:
    records = masker.filter_records(sam)

logging.info("Creating classifier")
classifier = PrecisionClassifier(sam=records, name=sample_id)

logging.info("Creating reporter")
reporter = PrecisionReporter(classifiers=[classifier])


# output
logging.info("Generating and saving report")
with open(variant_call_precision_report, "w") as output:
    reporter.save(output)

logging.info("Done")