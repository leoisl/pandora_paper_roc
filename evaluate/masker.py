from intervaltree import IntervalTree, Interval
from evaluate.classification import Classification
from typing import TextIO, Iterable, Type, Optional


class Masker:
    def __init__(self, tree: IntervalTree = None):
        if tree is None:
            tree = IntervalTree()
        self.tree = tree

    def __eq__(self, other: "Masker") -> bool:
        return self.tree == other.tree

    @staticmethod
    def from_bed(bed: TextIO) -> "Masker":
        tree = IntervalTree()
        for region in bed:
            chrom, start, end = region.strip().split("\t")
            tree.addi(int(start), int(end), chrom)
        return Masker(tree=tree)

    def filter_records(
        self, classifications: Iterable[Type[Classification]]
    ) -> Iterable[Type[Classification]]:
        return [
            classification
            for classification in classifications
            if not self.record_overlaps_mask(classification)
        ]

    def record_overlaps_mask(self, record: Type[Classification]) -> bool:
        interval = self.get_interval_where_probe_aligns_to_truth(record)
        if interval is None:
            return False

        overlaps = self.tree.overlap(interval)
        return any(interval.data == iv.data for iv in overlaps)

    @staticmethod
    def get_interval_where_probe_aligns_to_truth(record: Classification) -> Interval:
        raise NotImplementedError()


class PrecisionMasker(Masker):
    @staticmethod
    def get_interval_where_probe_aligns_to_truth(
        record: Classification
    ) -> Optional[Interval]:
        if record.is_unmapped:
            return None

        probe_aligned_pairs = record.get_probe_aligned_pairs()
        _, ref_positions, _ = zip(*probe_aligned_pairs)
        ref_positions = sorted(pos for pos in ref_positions if pos)
        # todo: what to do if ref_positions is all None?
        # could get aligned pairs and keep track of ref positions up until the first
        # query positions from the probe aligned pairs and then get the next non-None
        # ref position after it
        chromosome = record.ref_probe.chrom

        return Interval(ref_positions[0], ref_positions[-1] + 1, chromosome)
