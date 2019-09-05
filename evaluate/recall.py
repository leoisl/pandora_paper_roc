from typing import TextIO, Iterable
import logging

import pandas as pd

from evaluate.classifier import RecallClassifier
from .classification import *


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
                truth_probe_header = str(classification.query_probe.header)
                vcf_probe_header = str(classification.ref_probe.header)
                report_entries.append(
                    [classifier.name, truth_probe_header, vcf_probe_header, assessment]
                )

        return pd.DataFrame(data=report_entries, columns=self.columns)

    def save(self, file_handle: TextIO) -> pd.DataFrame:
        report = self.generate_report()
        report.to_csv(file_handle, sep=self.delim, header=True, index=False)
        return report


class StatisticalClassification(Enum):
    FALSE_NEGATIVE = "fn"
    FALSE_POSITIVE = "fp"
    TRUE_POSITIVE = "tp"
    TRUE_NEGATIVE = "tn"


class RecallCalculator:
    def __init__(self, reports: Iterable[pd.DataFrame]):
        self.report = pd.concat(reports)

    def calculate_recall(self, conf_threshold: float = 0) -> float:
        counter = Counter()
        for index, row in self.report.iterrows():
            classification = self.statistical_classification(row, conf_threshold)
            counter[classification] += 1

        true_positives = counter[StatisticalClassification.TRUE_POSITIVE]
        false_negatives = counter[StatisticalClassification.FALSE_NEGATIVE]

        logging.info(
            (
                f"Got {true_positives} true positives and {false_negatives}"
                " false negatives when calculating recall."
            )
        )

        try:
            return true_positives / (true_positives + false_negatives)
        except ZeroDivisionError:
            return 0

    @staticmethod
    def statistical_classification(
        row: pd.Series, conf_threshold: float = 0
    ) -> StatisticalClassification:
        gt_conf = ProbeHeader.from_string(row.vcf_probe_header).gt_conf
        if gt_conf < conf_threshold:
            return StatisticalClassification.FALSE_NEGATIVE
        else:
            return {
                AlignmentAssessment.UNMAPPED: StatisticalClassification.FALSE_NEGATIVE,
                AlignmentAssessment.PARTIALLY_MAPPED: StatisticalClassification.FALSE_NEGATIVE,
                AlignmentAssessment.PRIMARY_CORRECT: StatisticalClassification.TRUE_POSITIVE,
                AlignmentAssessment.PRIMARY_INCORRECT: StatisticalClassification.FALSE_POSITIVE,
                AlignmentAssessment.SECONDARY_CORRECT: StatisticalClassification.TRUE_POSITIVE,
                AlignmentAssessment.SECONDARY_INCORRECT: StatisticalClassification.FALSE_POSITIVE,
                AlignmentAssessment.SUPPLEMENTARY_CORRECT: StatisticalClassification.TRUE_POSITIVE,
                AlignmentAssessment.SUPPLEMENTARY_INCORRECT: StatisticalClassification.FALSE_POSITIVE,
            }[AlignmentAssessment(row.classification)]
