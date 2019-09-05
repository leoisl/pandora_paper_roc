from typing import Iterable
import logging

import pandas as pd

from .classification import *


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
