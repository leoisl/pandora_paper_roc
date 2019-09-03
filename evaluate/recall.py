from typing import List, TextIO, Iterable
from enum import Enum

import pandas as pd
import pysam

from .probe import Probe, ProbeHeader


class RecallClassification:
    def __init__(self, record: pysam.AlignedSegment = None):
        self.record = record

        if self.record is not None:
            self.truth_probe = Probe(
                header=ProbeHeader.from_string(self.record.query_name),
                full_sequence=self.record.query_sequence,
            )
            reference_name = self.record.reference_name or ""
            self.vcf_probe = Probe(header=ProbeHeader.from_string(reference_name))
        else:
            self.truth_probe = Probe()
            self.vcf_probe = Probe()

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
    def __init__(self, sam: Iterable[pysam.AlignedSegment] = None, name: str = ""):
        if sam is None:
            sam = []
        self.sam = sam
        self.name = name

    def classify(self) -> List[RecallClassification]:
        classifications = []
        for record in self.sam:
            classification = RecallClassification(record=record)
            classifications.append(classification)

        return classifications


class RecallReporter:
    def __init__(self, classifiers: Iterable[RecallClassifier], delim: str = "\t"):
        self.classifiers = classifiers
        self.delim = delim
        self.columns = [
            "sample",
            "truth_probe_header",
            "vcf_probe_header",
            "classification",
        ]

    def generate_report(self) -> pd.DataFrame:
        report_entries = []
        for classifier in self.classifiers:
            classifications = classifier.classify()
            for classification in classifications:
                assessment = classification.assessment()
                truth_probe_header = str(classification.truth_probe.header)
                vcf_probe_header = str(classification.vcf_probe.header)
                report_entries.append(
                    [classifier.name, truth_probe_header, vcf_probe_header, assessment]
                )

        return pd.DataFrame(data=report_entries, columns=self.columns)

    def save(self, file_handle: TextIO) -> None:
        report = self.generate_report()
        report.to_csv(file_handle, sep=self.delim, header=True, index=False)


class StatisticalClassification(Enum):
    FALSE_NEGATIVE = "fn"
    FALSE_POSITIVE = "fp"
    TRUE_POSITIVE = "tp"
    TRUE_NEGATIVE = "tn"


class RecallCalculator:
    @staticmethod
    def statistical_classification(
        row: pd.Series, conf_threshold: float = 0
    ) -> StatisticalClassification:
        gt_conf = ProbeHeader.from_string(row.vcf_probe_header).gt_conf
        if gt_conf < conf_threshold:
            return StatisticalClassification.FALSE_NEGATIVE
        else:
            return {
                "unmapped": StatisticalClassification.FALSE_NEGATIVE,
                "partially_mapped": StatisticalClassification.FALSE_NEGATIVE,
                "correct": StatisticalClassification.TRUE_POSITIVE,
                "incorrect": StatisticalClassification.FALSE_POSITIVE,
                "secondary_correct": StatisticalClassification.TRUE_POSITIVE,
                "secondary_incorrect": StatisticalClassification.FALSE_POSITIVE,
                "supplementary_correct": StatisticalClassification.TRUE_POSITIVE,
                "supplementary_incorrect": StatisticalClassification.FALSE_POSITIVE,
            }[row.classification]
