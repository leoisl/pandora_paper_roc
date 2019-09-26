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
sam_file = snakemake.input.variant_call_probeset_mapped_to_ref
sample_id = snakemake.wildcards.sample_id
mask = snakemake.input.mask
variant_call_precision_report = snakemake.output.variant_call_precision_report


# API usage
with pysam.AlignmentFile(sam_file) as variant_call_probeset_mapped_to_ref:
    logging.info(f"Creating mask from {mask}")
    with open(mask) as mask_file:
        precision_masker = PrecisionMasker.from_bed(mask_file)

    logging.info(f"Filtering sam records")
    filtered_sam_records = precision_masker.filter_records(
        variant_call_probeset_mapped_to_ref
    )

    logging.info(f"Creating precision report")
    precision_classifier = PrecisionClassifier(filtered_sam_records, sample_id)
    precision_reporter = PrecisionReporter([precision_classifier])


    # output
    logging.info(f"Outputting precision report")
    precision_reporter.save(Path(variant_call_precision_report))

logging.info(f"Done")