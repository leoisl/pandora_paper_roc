from intervaltree import IntervalTree
from typing import TextIO, Iterable
import pysam


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
        self, records: Iterable[pysam.AlignedSegment]
    ) -> Iterable[pysam.AlignedSegment]:
        return [record for record in records if not self.record_overlaps_mask(record)]

    def record_overlaps_mask(self, record: pysam.AlignedSegment) -> bool:
        raise NotImplementedError("We're all adults here. Deal with the consequences.")
