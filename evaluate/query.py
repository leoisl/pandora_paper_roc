import pysam
from pathlib import Path
from typing import Tuple, List
from contextlib import ExitStack


class OverlappingRecordsError(Exception):
    pass


class Query:
    def __init__(self, vcf: Path, vcf_ref: Path, flank_width: int = 0):
        self.vcf = Path(
            pysam.tabix_index(str(vcf), preset="vcf", keep_original=True, force=True)
        )
        self.genes = vcf_ref
        self._probe_names = set()
        self.flank_width = flank_width
        self._entry_number = 0

    def make_probes(self) -> str:
        query_probes = ""
        with ExitStack() as stack:
            vcf = stack.enter_context(pysam.VariantFile(self.vcf))
            genes_fasta = stack.enter_context(pysam.FastxFile(str(self.genes)))

            for gene in genes_fasta:
                try:
                    entries = vcf.fetch(contig=gene.name)
                except ValueError as error:
                    if str(error).startswith("invalid contig"):
                        continue
                    else:
                        raise error

                probes_for_gene = self._create_probes_for_gene_variants(gene, entries)
                query_probes += probes_for_gene

        return query_probes

    def _create_probes_for_gene_variants(
        self, gene: pysam.FastxRecord, variants: pysam.tabix_iterator
    ) -> str:
        """Note: An assumption is made with this function that the variants you pass in
        are from the gene passed with them."""
        probes = ""
        variants = [entry for entry in variants if not is_invalid_vcf_entry(entry)]
        intervals_to_probes = dict()

        for variant in variants:
            interval = self.calculate_probe_boundaries_for_entry(variant)
            if interval in intervals_to_probes and float(
                intervals_to_probes[interval].name.split("=")[-1]
            ) > get_genotype_confidence(variant):
                continue

            mutated_consensus = ""
            consensus = gene.sequence[slice(*interval)]
            last_idx = 0

            start_idx_of_variant_on_consensus = variant.start - interval[0]
            mutated_consensus += consensus[last_idx:start_idx_of_variant_on_consensus]
            mutated_consensus += get_variant_sequence(variant)
            last_idx = start_idx_of_variant_on_consensus + variant.rlen
            mutated_consensus += consensus[last_idx:]
            probe = pysam.FastxRecord()
            probe.set_name(
                f"{variant.chrom}_POS={variant.pos}_interval={interval}_GT_CONF={get_genotype_confidence(variant)}".replace(
                    " ", ""
                )
            )
            probe.set_sequence(mutated_consensus)
            intervals_to_probes[interval] = probe

        for probe in intervals_to_probes.values():
            probes += str(probe) + "\n"

        return probes

    def calculate_probe_boundaries_for_entry(
        self, entry: pysam.VariantRecord
    ) -> Tuple[int, int]:
        probe_start = max(0, entry.start - self.flank_width)
        probe_stop = entry.stop + self.flank_width

        return probe_start, probe_stop


def merge_overlap_intervals(intervals: List[List[int]]) -> List[Tuple[int, int]]:
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


def is_invalid_vcf_entry(entry: pysam.VariantRecord) -> bool:
    genotype = get_genotype(entry)

    return genotype is None


def get_genotype_confidence(variant: pysam.VariantRecord) -> float:
    return float(variant.samples["sample"].get("GT_CONF", 0))


def get_genotype(variant: pysam.VariantRecord) -> int:
    samples = list(variant.samples.keys())
    assert len(samples) == 1
    return variant.samples[samples[0]]["GT"][0]


def get_variant_sequence(variant: pysam.VariantRecord) -> str:
    genotype = get_genotype(variant)

    if genotype is None:
        return variant.ref
    else:
        return variant.alleles[genotype]


def get_variant_length(variant: pysam.VariantRecord) -> int:
    return len(get_variant_sequence(variant))
