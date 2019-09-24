import logging
from contextlib import ExitStack
from pathlib import Path
from typing import Tuple, List, Dict

import pysam

from .probe import ProbeHeader, Probe, ProbeInterval
from .vcf import VCF
from .vcf_file import VCFFile

class OverlappingRecordsError(Exception):
    pass


class Query:
    def __init__(
        self, vcf_filepath: Path, vcf_ref: Path, samples: List[str], flank_width: int = 0
    ):
        self.vcf_filepath = vcf_filepath
        self.genes = vcf_ref
        self._probe_names = set()
        self.flank_width = flank_width
        self._entry_number = 0
        self.samples = samples

    def make_probes(self) -> Dict[str, str]:
        query_probes: Dict[str, str] = {s: "" for s in self.samples}
        with pysam.FastxFile(str(self.genes)) as genes_fasta:
            for gene in genes_fasta:
                for sample in self.samples:
                    try:
                        vcfs = self.vcf_file.get_VCF_records_given_sample_and_gene(sample, gene.name)
                    except ValueError as error:
                        if str(error).startswith("invalid contig"):
                            continue
                        else:
                            raise error

                    probes_for_gene: Dict[str, str] = self._create_probes_for_gene_variants(
                        gene, vcfs
                    )
                    for sample, probe in probes_for_gene.items():
                        query_probes[sample] += probe

        return query_probes

    # TODO : tagged for refactoring - this function does a lot of things
    def _create_probes_for_gene_variants(
        self, gene: pysam.FastxRecord, vcfs: List[VCF]
    ) -> Dict[str, str]:
        """Note: An assumption is made with this function that the variants you pass in
        are from the gene passed with them."""

        sample_to_probes: Dict[str, str] = {s: "" for s in self.samples}
        sample_to_intervals_to_probes: Dict[str, Dict[ProbeInterval, Probe]] = {
            s: {} for s in self.samples
        }

        for vcf in vcfs:
            sample = vcf.sample
            if vcf.is_invalid_vcf_entry:
                continue
            interval = self.calculate_probe_boundaries_for_entry(vcf)

            if (
                interval in sample_to_intervals_to_probes[sample]
                and sample_to_intervals_to_probes[sample][interval].gt_conf
                > vcf.genotype_confidence
            ):
                continue

            mutated_consensus = ""
            consensus = gene.sequence[slice(*interval)]
            last_idx = 0

            start_idx_of_variant_on_consensus = vcf.start - interval.start
            mutated_consensus += consensus[
                last_idx:start_idx_of_variant_on_consensus
            ]
            mutated_consensus += vcf.variant_sequence
            last_idx = start_idx_of_variant_on_consensus + vcf.rlen
            mutated_consensus += consensus[last_idx:]
            probe_header = self._create_probe_header(sample, vcf, interval)
            probe = Probe(header=probe_header, full_sequence=mutated_consensus)
            if sample not in sample_to_intervals_to_probes:
                sample_to_intervals_to_probes[sample] = {interval: probe}
            else:
                sample_to_intervals_to_probes[sample][interval] = probe

        for sample in sample_to_intervals_to_probes:
            for interval, probe in sample_to_intervals_to_probes[sample].items():
                sample_to_probes[sample] += str(probe) + "\n"

        return sample_to_probes

    def calculate_probe_boundaries_for_entry(
        self, vcf: VCF
    ) -> ProbeInterval:
        probe_start = max(0, vcf.start - self.flank_width)
        probe_stop = vcf.stop + self.flank_width

        return ProbeInterval(probe_start, probe_stop)

    @staticmethod
    def _create_probe_header(
        sample: str, vcf: VCF, interval: ProbeInterval
    ) -> ProbeHeader:
        call_start_idx = max(0, vcf.start - interval[0])
        call_end_idx = call_start_idx + vcf.variant_length
        call_interval = ProbeInterval(call_start_idx, call_end_idx)
        return ProbeHeader(
            chrom=vcf.chrom,
            sample=sample,
            pos=vcf.pos,
            interval=call_interval,
            svtype=vcf.svtype,
            mean_fwd_covg=vcf.mean_coverage_forward,
            mean_rev_covg=vcf.mean_coverage_reverse,
            gt_conf=vcf.genotype_confidence,
        )


# TODO: refactor all these functions into an Intervals class
def merge_overlap_intervals(intervals: List[List[int]]) -> List[Tuple[int, ...]]:
    """Checks consecutive intervals and if they overlap it merges them into a
    single interval.
    Args:
        intervals: A list of intervals where each interval is a List with two
        elements corresponding to the start and end of the interval
        respectively.
    Returns:
        A new intervals list where any intervals that overlapped have been
        merged into a single interval.
    Example:
        >>> intervals = [[1, 4], [3, 7], [10, 14]]
        >>> merge_overlap_intervals(intervals)
        [(1, 7), (10, 14)]
    """
    merged_intervals = []
    cached_interval = None

    for interval in intervals:
        if cached_interval is None:
            cached_interval = interval
            continue

        if outside_interval(cached_interval, interval):
            merged_intervals.append(tuple(cached_interval))
            cached_interval = interval
        else:
            cached_interval = extend_interval(cached_interval, interval)

    if cached_interval is not None:
        merged_intervals.append(tuple(cached_interval))

    return merged_intervals


def outside_interval(first_interval: List[int], second_interval: List[int]) -> bool:
    """Determines whether two intervals overlap.
    Args:
        first_interval: The interval with the lower start index.
        second_interval: The interval with the higher start index.
    Returns:
        Whether the start of the second interval is less than the end of the
        first interval. i.e do they overlap?
    Notes:
        If the end index of the first interval is equal to the start of the
        second interval than they are deemed to NOT be overlapping.
    Example:
        >>> first_interval = [0, 4]
        >>> second_interval = [3, 7]
        >>> outside_interval(first_interval, second_interval)
        False
    """
    return second_interval[0] >= first_interval[1]


def extend_interval(interval_to_extend: List[int], interval: List[int]) -> List[int]:
    """Extends an interval to encompass another.
    Args:
        interval_to_extend: The interval to extend.
        interval: The interval to extend by.
    Returns:
        A new interval with the same start as interval_to_extend and the same
        end as interval.
    """
    interval_to_extend[1] = interval[1]

    return interval_to_extend


def find_index_in_intervals(intervals: List[Tuple[int, int]], query: int) -> int:
    """Return the index of the interval that the query lies within"""
    for i, (start, end) in enumerate(intervals):
        if start <= query <= end:
            return i
    return -1


def write_vcf_probes_to_file(
    vcf_probes: Dict[str, str], query_name: str, tempdir: Path
) -> Path:
    query_vcf_probes = vcf_probes[query_name]
    query_vcf_probes_path: Path = tempdir / f"{query_name}.query_probes.fa"
    query_vcf_probes_path.write_text(query_vcf_probes)
    logging.info(f"VCF probes written to file: {query_vcf_probes_path}")
    return query_vcf_probes_path
