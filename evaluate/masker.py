from intervaltree import IntervalTree, Interval
import math
from evaluate.classification import Classification
from typing import TextIO, Type, Optional
import pysam
from evaluate.filter import Filter


class Masker(Filter):
    def __init__(self, tree: IntervalTree = None):
        if tree is None:
            tree = IntervalTree()
        self.tree = tree

    def __eq__(self, other: "Masker") -> bool:
        return self.tree == other.tree

    @classmethod
    def from_bed(cls: Type, bed: TextIO) -> "Type[Masker]":
        tree = IntervalTree()
        for region in bed:
            chrom, start, end = region.strip().split("\t")
            tree.addi(int(start), int(end), chrom)
        return cls(tree=tree)

    def record_should_be_filtered_out(self, record: pysam.AlignedSegment) -> bool:
        return self.record_overlaps_mask(record)

    def record_overlaps_mask(self, record: pysam.AlignedSegment) -> bool:
        classification = Classification(record)
        interval = self.get_interval_where_probe_aligns_to_truth(classification)
        if interval is None:
            return False

        overlaps = self.tree.overlap(interval)
        return any(interval.data == iv.data for iv in overlaps)

    @staticmethod
    def get_interval_where_probe_aligns_to_truth(
        record: Classification
    ) -> Optional[Interval]:
        raise NotImplementedError()


class PrecisionMasker(Masker):
    @staticmethod
    def get_interval_where_probe_aligns_to_truth(
        record: Classification
    ) -> Optional[Interval]:
        if record.is_unmapped:
            return None

        aligned_pairs = record.get_aligned_pairs(with_seq=True)
        query_interval = aligned_pairs.get_index_of_query_interval(
            Interval(*record.query_probe.get_interval_or_default_interval_if_none())
        )
        ref_positions = aligned_pairs.get_ref_positions(
            transform_Nones_into_halfway_positions=True
        )
        ref_positions_query_aligns_to = ref_positions[slice(*query_interval)]
        if len(ref_positions_query_aligns_to)==0:
            return None

        ref_start, ref_end = (
            math.floor(ref_positions_query_aligns_to[0]),
            math.ceil(ref_positions_query_aligns_to[-1]),
        )
        chromosome = record.ref_probe.chrom

        return Interval(max(0, ref_start), ref_end + 1, chromosome)


class RecallMasker(Masker):
    @staticmethod
    def get_interval_where_probe_aligns_to_truth(record: Classification) -> Interval:
        begin = record.query_probe.pos
        end = begin + len(record.query_probe.get_interval(False))
        chrom = record.query_probe.chrom
        return Interval(begin, end, data=chrom)
