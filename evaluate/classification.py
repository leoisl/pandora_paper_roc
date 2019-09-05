from collections import Counter
from enum import Enum
from typing import List
import pysam

from .probe import Probe, ProbeHeader, DELIM


class AlignmentType(Enum):
    INSERTION = 0
    DELETION = 1
    MISMATCH = 2
    MATCH = 3


class Classification:
    def __init__(self, record: pysam.AlignedSegment = None):
        self.record = record

        if self.record is not None:
            self.query_probe = Probe(
                header=ProbeHeader.from_string(self.record.query_name),
                full_sequence=self.record.query_sequence,
            )

            reference_name = self.record.reference_name or ""
            if reference_name.startswith(">"):
                reference_name = self.record.reference_name
            else:
                reference_name = f"CHROM={reference_name}{DELIM}"
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

    def _whole_query_probe_maps(self) -> bool:
        if self.is_unmapped:
            return False

        query_interval = self.query_probe.interval
        query_starts_before_alignment = (
            query_interval.start < self.query_alignment_start
        )
        query_ends_after_alignment = query_interval.end > self.query_alignment_end

        if query_starts_before_alignment or query_ends_after_alignment:
            return False

        return True

    @staticmethod
    def __get_alignment_type(query_pos: int, ref_base: str) -> AlignmentType:
        if ref_base is None:
            return AlignmentType.DELETION
        elif query_pos is None:
            return AlignmentType.INSERTION
        elif ref_base.islower():
            return AlignmentType.MISMATCH
        else:
            return AlignmentType.MATCH

    def get_query_probe_mapping_score(self) -> float:
        if not self.query_probe.is_deletion:
            query_start = self.query_probe.interval.start
            query_stop = self.query_probe.interval.end - 1
        else:
            query_start = max(0, self.query_probe.interval.start - 1)
            query_stop = self.query_probe.interval.end

        within_probe = False
        alignment_type_count = Counter()
        for query_pos, ref_pos, ref_base in self.get_aligned_pairs(with_seq=True):
            arrived_in_the_probe = query_pos is not None and query_pos == query_start
            if arrived_in_the_probe:
                within_probe = True

            if within_probe:
                alignment_type = self.__get_alignment_type(query_pos, ref_base)
                alignment_type_count[alignment_type] += 1

            exiting_probe = query_pos is not None and query_pos == query_stop
            if exiting_probe:
                break

        total_nb_of_alignments_checked = sum(alignment_type_count.values())
        query_probe_mapping_score = (
            alignment_type_count[AlignmentType.MATCH] / total_nb_of_alignments_checked
        )
        return query_probe_mapping_score

    def assessment(self) -> str:
        raise NotImplementedError()


class AlignmentAssessment(Enum):
    UNMAPPED = "unmapped"
    PARTIALLY_MAPPED = "partially_mapped"
    PRIMARY_CORRECT = "primary_correct"
    PRIMARY_INCORRECT = "primary_incorrect"
    SECONDARY_CORRECT = "secondary_correct"
    SECONDARY_INCORRECT = "secondary_incorrect"
    SUPPLEMENTARY_CORRECT = "supplementary_correct"
    SUPPLEMENTARY_INCORRECT = "supplementary_incorrect"


class RecallClassification(Classification):
    def is_correct(self) -> bool:
        return self.get_query_probe_mapping_score() == 1.0

    def assessment(self) -> AlignmentAssessment:
        if self.is_unmapped:
            assessment = AlignmentAssessment.UNMAPPED
        elif not self._whole_query_probe_maps():
            assessment = AlignmentAssessment.PARTIALLY_MAPPED
        else:
            is_correct = self.is_correct()
            if self.is_secondary:
                assessment = (
                    AlignmentAssessment.SECONDARY_CORRECT
                    if is_correct
                    else AlignmentAssessment.SECONDARY_INCORRECT
                )
            elif self.is_supplementary:
                assessment = (
                    AlignmentAssessment.SUPPLEMENTARY_CORRECT
                    if is_correct
                    else AlignmentAssessment.SUPPLEMENTARY_INCORRECT
                )
            else:
                assessment = (
                    AlignmentAssessment.PRIMARY_CORRECT
                    if is_correct
                    else AlignmentAssessment.PRIMARY_INCORRECT
                )

        return assessment


class PrecisionClassification(Classification):
    def assessment(self) -> float:
        return self.get_query_probe_mapping_score()
