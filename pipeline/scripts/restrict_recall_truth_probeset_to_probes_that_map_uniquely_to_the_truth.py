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
from evaluate.bwa import BWA
import pysam
from evaluate.unique_sam_records_filter import UniqueSamRecordsFilter
import pandas as pd

# setup
sample_id = Path(snakemake.wildcards.sample_id)
sample_pair = Path(snakemake.wildcards.sample_pair)

truth_filepath = Path(snakemake.input.truth)
unrestricted_probeset_filepath = Path(snakemake.input.unrestricted_probeset)

unrestricted_probeset_mapped_to_truth_sam_filepath = Path(snakemake.output.unrestricted_probeset_mapped_to_truth_sam)
nb_of_truth_probes_removed_with_unique_sam_records_filter_filepath = Path(snakemake.output.nb_of_truth_probes_removed_with_unique_sam_records_filter_filepath)
restricted_probeset_filepath = Path(snakemake.output.restricted_probeset)

threads = int(snakemake.threads)


# API usage

# just in order to keep the exact same sequence - if we get the sequence from the SAM file, it can be rev complemented
unrestricted_probeset_header_to_sequence = {}
with pysam.FastxFile(unrestricted_probeset_filepath) as unrestricted_probeset_file:
    for record in unrestricted_probeset_file:
        unrestricted_probeset_header_to_sequence[record.name] = record.sequence


logging.info(f"Mapping {unrestricted_probeset_filepath} to {truth_filepath}")
BWA.map_query_to_ref(
    query=unrestricted_probeset_filepath,
    ref=truth_filepath,
    output=unrestricted_probeset_mapped_to_truth_sam_filepath,
    threads=threads,
)

# get all mappings
logging.info("Filtering truth probes to unique truth probes only")
with pysam.AlignmentFile(unrestricted_probeset_mapped_to_truth_sam_filepath) as sam:
    records = [record for record in sam]

# create the UniqueSamRecordsFilter
unique_sam_records_filter = UniqueSamRecordsFilter(records)

# output only the truth probes that are uniquely mapped
nb_of_truth_probes_before_unique_sam_records_filter = len(records)
nb_of_truth_probes_after_unique_sam_records_filter = 0
with open(restricted_probeset_filepath, "w") as restricted_probeset_file:
    for record in records:
        if not unique_sam_records_filter.record_should_be_filtered_out(record):
            print(f">{record.query_name}\n{unrestricted_probeset_header_to_sequence[record.query_name]}", file=restricted_probeset_file)
            nb_of_truth_probes_after_unique_sam_records_filter += 1

nb_of_truth_probes_removed_with_unique_sam_records_filter = nb_of_truth_probes_before_unique_sam_records_filter - nb_of_truth_probes_after_unique_sam_records_filter
nb_of_truth_probes_removed_with_unique_sam_records_filter_df = pd.DataFrame({
    "sample_id": [sample_id],
    "sample_pair": [sample_pair],
    "nb_of_truth_probes_before_unique_sam_records_filter": [nb_of_truth_probes_before_unique_sam_records_filter],
    "nb_of_truth_probes_after_unique_sam_records_filter": [nb_of_truth_probes_after_unique_sam_records_filter],
    "nb_of_truth_probes_removed_with_unique_sam_records_filter": [nb_of_truth_probes_removed_with_unique_sam_records_filter],
    "nb_of_truth_probes_removed_with_unique_sam_records_filter_proportion": [nb_of_truth_probes_removed_with_unique_sam_records_filter/nb_of_truth_probes_before_unique_sam_records_filter]
})
nb_of_truth_probes_removed_with_unique_sam_records_filter_df.to_csv(nb_of_truth_probes_removed_with_unique_sam_records_filter_filepath, index=False)
