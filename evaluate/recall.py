from typing import List, Tuple
import pysam

from evaluate.probe import ProbeHeader, Probe
from .probe import Probe, Interval, ProbeHeader


class RecallClassification:
    def __init__(self, truth_probe: Probe = None, record: pysam.AlignedSegment = None):
        if truth_probe is None:
            truth_probe = Probe()

        self.truth_probe = truth_probe
        self.record = record

    def __eq__(self, other: "RecallClassification") -> bool:
        return self.truth_probe == other.truth_probe

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

        truth_interval = self.truth_probe.interval
        truth_starts_before_alignment = (
            truth_interval.start < self.query_alignment_start
        )
        truth_ends_after_alignment = truth_interval.end > self.query_alignment_end

        if truth_starts_before_alignment or truth_ends_after_alignment:
            return False

        return True

    def is_correct(self) -> bool:
        within_probe = False
        ref_seq = ""
        # todo: do this more succinctly
        if len(self.truth_probe.interval) > 0:
            truth_start = self.truth_probe.interval.start
            truth_stop = self.truth_probe.interval.end - 1
            truth = self.truth_probe.core_sequence
        else:
            truth_start = max(0, self.truth_probe.interval.start - 1)
            truth_stop = self.truth_probe.interval.end
            truth = self.truth_probe.full_sequence[truth_start : truth_stop + 1]

        for query_pos, ref_pos, ref_base in self.get_aligned_pairs(with_seq=True):
            if query_pos is not None and query_pos == truth_start:
                within_probe = True
            if query_pos is not None and query_pos == truth_stop:
                if ref_base is None:
                    return False
                else:
                    ref_seq += ref_base
                break

            if within_probe:
                if ref_base is None:
                    return False
                else:
                    ref_seq += ref_base

        return ref_seq == truth

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


class RecallClassifier:
    def __init__(self, sam: pysam.AlignmentFile):
        self.sam = sam

    def classify(self) -> List[RecallClassification]:
        classifications = []
        for record in self.sam:
            truth_probe = Probe(
                header=ProbeHeader.from_string(record.query_name),
                full_sequence=record.query_sequence,
            )
            classification = RecallClassification(
                truth_probe=truth_probe, record=record
            )
            classifications.append(classification)

        return classifications
