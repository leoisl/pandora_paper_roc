from intervaltree import IntervalTree
from typing import TextIO


class Masker:
    def __init__(self, tree: IntervalTree = None):
        self.tree = tree

    def __eq__(self, other: "Masker") -> bool:
        return self.tree == other.tree

    @staticmethod
    def from_bed(bed: TextIO) -> "Masker":
        tree = IntervalTree()
        for region in bed:
            chrom, start, end = region.strip().split("\t")
            tree.addi(int(start), int(end), chrom)
        return Masker()
