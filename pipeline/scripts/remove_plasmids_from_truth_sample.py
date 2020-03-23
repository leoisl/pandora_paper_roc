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


# setup
truth_sample = snakemake.input.truth_sample
truth_sample_with_chrom_only = snakemake.output.truth_sample_with_chrom_only
threshold = 0.1


with pysam.FastxFile(truth_sample) as truth_sample_filehandler:
    entries = [entry for entry in truth_sample_filehandler]

sequences_lengths = [len(entry.sequence) for entry in entries]
sequences_lengths_max_to_min = sorted(sequences_lengths, reverse=True)
longest_sequence_length = sequences_lengths_max_to_min[0]
second_longest_sequence_length = sequences_lengths_max_to_min[1] if len(sequences_lengths_max_to_min) > 1 else 0

longest_sequence_is_way_longer_than_second_longest_sequence = (0.1 * longest_sequence_length) > second_longest_sequence_length
assert longest_sequence_is_way_longer_than_second_longest_sequence,\
    f"In sample {truth_sample}, the longest sequence has {longest_sequence_length} bps and " \
    f"the second longest sequences has {second_longest_sequence_length}. " \
    f"The longest sequence is assumed to be the bacterial chromosome, and the rest to be plasmids. " \
    f"Plasmids are supposed to be at most {threshold*100}% the length of the chromosome (in this pipeline). " \
    f"Are you sure your genome is well assembled?"

longest_entry = [entry for entry in entries if len(entry.sequence)==longest_sequence_length][0]
with open(truth_sample_with_chrom_only, "w") as truth_sample_with_chrom_only_filehandler:
    print(longest_entry, file=truth_sample_with_chrom_only_filehandler)