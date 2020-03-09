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
from evaluate.mapq_sam_records_filter import MAPQSamRecordsFilter
import pandas as pd


# setup
sam_filepath = snakemake.input.variant_call_probeset_mapped_to_ref
sample_id = snakemake.wildcards.sample_id
tool = snakemake.wildcards.tool
mask_filepath = snakemake.input.mask
variant_call_precision_report = snakemake.output.variant_call_precision_report
nb_of_records_removed_with_mapq_sam_records_filter_filepath = snakemake.output.nb_of_records_removed_with_mapq_sam_records_filter_filepath


# API usage
with pysam.AlignmentFile(sam_filepath) as sam:
    records = [record for record in sam]

logging.info(f"Applying MAPQ SAM records filter")
nb_of_records_before_mapq_sam_records_filter = len(records)
mapq_sam_records_filter = MAPQSamRecordsFilter(records)
records = mapq_sam_records_filter.filter_records(records)
nb_of_records_after_mapq_sam_records_filter = len(records)
nb_of_records_removed_with_mapq_sam_records_filter = nb_of_records_before_mapq_sam_records_filter - nb_of_records_after_mapq_sam_records_filter
nb_of_records_removed_with_mapq_sam_records_filter_proportion = nb_of_records_removed_with_mapq_sam_records_filter/nb_of_records_before_mapq_sam_records_filter if nb_of_records_before_mapq_sam_records_filter>0 else 0

nb_of_records_removed_with_mapq_sam_records_filter_df = pd.DataFrame({
    "tool": [tool],
    "nb_of_records_before_mapq_sam_records_filter": [nb_of_records_before_mapq_sam_records_filter],
    "nb_of_records_after_mapq_sam_records_filter": [nb_of_records_after_mapq_sam_records_filter],
    "nb_of_records_removed_with_mapq_sam_records_filter": [nb_of_records_removed_with_mapq_sam_records_filter],
    "nb_of_records_removed_with_mapq_sam_records_filter_proportion": [nb_of_records_removed_with_mapq_sam_records_filter_proportion]
})
nb_of_records_removed_with_mapq_sam_records_filter_df.to_csv(nb_of_records_removed_with_mapq_sam_records_filter_filepath, index=False)

logging.info(f"Masking SAM records")
with open(mask_filepath) as bed:
    masker = PrecisionMasker.from_bed(bed)
records = masker.filter_records(records)

logging.info("Creating classifier")
classifier = PrecisionClassifier(sam=records, name=sample_id)

logging.info("Creating reporter")
reporter = PrecisionReporter(classifiers=[classifier])

logging.info("Generating report")
report = reporter.generate_report()


# output
logging.info("Saving report")
with open(variant_call_precision_report, "w") as output:
    reporter.save_report(report, output)

logging.info("Done")