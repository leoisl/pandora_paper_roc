from collections import Counter
from enum import Enum

import pysam
from intervaltree import Interval

from evaluate.aligned_pairs import AlignmentType, AlignedPairs
from .probe import Probe, ProbeHeader, DELIM


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
    ) -> AlignedPairs:
        return AlignedPairs(
            self.record.get_aligned_pairs(matches_only=matches_only, with_seq=with_seq)
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

    def _get_query_probe_mapping_score(self) -> float:
        probe_aligned_pairs = self.get_probe_aligned_pairs()
        alignment_types = probe_aligned_pairs.get_alignment_types()
        alignment_type_count = Counter(alignment_types)

        total_nb_of_alignments_checked = len(probe_aligned_pairs)
        query_probe_mapping_score = (
            alignment_type_count[AlignmentType.MATCH] / total_nb_of_alignments_checked
        )

        assert 0.0 <= query_probe_mapping_score <= 1.0

        return query_probe_mapping_score

    def get_probe_aligned_pairs(self) -> AlignedPairs:
        if not self.query_probe.is_deletion:
            query_start = self.query_probe.interval.start
            query_stop = self.query_probe.interval.end
        else:
            query_start = max(0, self.query_probe.interval.start - 1)
            query_stop = self.query_probe.interval.end + 1

        probe_aligned_pairs = AlignedPairs(self.get_aligned_pairs(with_seq=True))
        return probe_aligned_pairs.get_pairs_in_query_interval(
            Interval(query_start, query_stop)
        )

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

    def __str__(self):
        return self.value

class RecallClassification(Classification):
    def is_correct(self) -> bool:
        return self._get_query_probe_mapping_score() == 1.0

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
        query_probe_does_not_map_completely = (
            self.is_unmapped or not self._whole_query_probe_maps()
        )

        if query_probe_does_not_map_completely:
            assessment = 0.0
        else:
            assessment = self._get_query_probe_mapping_score()

        return assessment
