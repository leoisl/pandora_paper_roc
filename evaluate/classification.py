from collections import defaultdict
from enum import Enum
from typing import List

import pysam

from .probe import Probe, ProbeHeader, DELIM


class Classification:
    def __init__(self, record: pysam.AlignedSegment = None):
        self.record = record

        if self.record is not None:
            self.query_probe = Probe(
                header=ProbeHeader.from_string(self.record.query_name),
                full_sequence=self.record.query_sequence,
            )
            reference_name = f"CHROM={self.record.reference_name or ''}{DELIM}"
            self.ref_probe = Probe(header=ProbeHeader.from_string(reference_name))
        else:
            self.query_probe = Probe()
            self.ref_probe = Probe()

    def __eq__(self, other: "Classification") -> bool:
        return self.query_probe == other.query_probe

    @property
    def is_unmapped(self) -> bool:
        return self.record.is_unmapped

    @property
    def is_secondary(self) -> bool:
        return self.record.is_secondary

    @property
    def is_supplementary(self) -> bool:
        return self.record.is_supplementary

    def get_aligned_pairs(
        self, matches_only: bool = False, with_seq: bool = False
    ) -> List[tuple]:
        return self.record.get_aligned_pairs(
            matches_only=matches_only, with_seq=with_seq
        )

    @property
    def query_alignment_start(self) -> int:
        return self.record.query_alignment_start

    @property
    def query_alignment_end(self) -> int:
        return self.record.query_alignment_end

    def _whole_probe_maps(self) -> bool:
        if self.is_unmapped:
            return False

        truth_interval = self.query_probe.interval
        truth_starts_before_alignment = (
            truth_interval.start < self.query_alignment_start
        )
        truth_ends_after_alignment = truth_interval.end > self.query_alignment_end

        if truth_starts_before_alignment or truth_ends_after_alignment:
            return False

        return True

    class AlignmentType(Enum):
        INSERTION = 0
        DELETION = 1
        MISMATCH = 2
        MATCH = 3

    @staticmethod
    def __get_alignment_type(
        query_pos: int, query_base: str, ref_pos: int, ref_base: str
    ) -> AlignmentType:
        # infers which alignments are present
        alignment_type_to_present_flag = {
            Classification.AlignmentType.INSERTION: ref_pos is None
            and ref_base is None,
            Classification.AlignmentType.DELETION: query_pos is None
            and query_base is None,
            Classification.AlignmentType.MISMATCH: query_base is not None
            and ref_base is not None
            and query_base != ref_base,
            Classification.AlignmentType.MATCH: query_base is not None
            and ref_base is not None
            and query_base == ref_base,
        }
        only_one_flag_is_set = (
            sum([int(flag) for flag in alignment_type_to_present_flag.values()]) == 1
        )
        assert only_one_flag_is_set, "Error, inferring more than one alignment type"

        # get the only alignment type set
        alignment_type = [
            alignment_type
            for alignment_type, present_flag in alignment_type_to_present_flag.items()
            if present_flag
        ][0]
        return alignment_type

    def get_query_probe_mapping_score(self) -> float:
        if not self.query_probe.is_deletion:
            query_start = self.query_probe.interval.start
            query_stop = self.query_probe.interval.end - 1
            query_sequence = self.query_probe.core_sequence
        else:
            query_start = max(0, self.query_probe.interval.start - 1)
            query_stop = self.query_probe.interval.end
            query_sequence = self.query_probe.full_sequence[
                query_start : query_stop + 1
            ]

        within_probe = False
        alignment_type_to_count = defaultdict(int)
        for query_pos, ref_pos, ref_base in self.get_aligned_pairs(with_seq=True):
            arrived_in_the_probe = query_pos is not None and query_pos == query_start
            if arrived_in_the_probe:
                within_probe = True

            if within_probe:
                query_base = (
                    None
                    if query_pos is None
                    else query_sequence[query_pos - query_start]
                )
                alignment_type = Classification.__get_alignment_type(
                    query_pos, query_base, ref_pos, ref_base
                )
                alignment_type_to_count[alignment_type] += 1

            exiting_probe = query_pos is not None and query_pos == query_stop
            if exiting_probe:
                break

        total_nb_of_alignments_checked = sum(alignment_type_to_count.values())
        query_probe_mapping_score = (
            alignment_type_to_count[Classification.AlignmentType.MATCH]
            / total_nb_of_alignments_checked
        )
        return query_probe_mapping_score

    def assessment(self) -> str:
        raise NotImplementedError()


class RecallClassification(Classification):
    def __init__(self, record: pysam.AlignedSegment = None):
        super().__init__(record)

    def is_correct(self) -> bool:
        return self.get_query_probe_mapping_score() == 1.0

    def assessment(self) -> str:
        if self.is_unmapped:
            assessment = "unmapped"
        elif not self._whole_probe_maps():
            assessment = "partially_mapped"
        else:
            is_correct = self.is_correct()
            if self.is_secondary:
                assessment = (
                    "secondary_correct" if is_correct else "secondary_incorrect"
                )
            elif self.is_supplementary:
                assessment = (
                    "supplementary_correct" if is_correct else "supplementary_incorrect"
                )
            else:
                assessment = "correct" if is_correct else "incorrect"

        return assessment
