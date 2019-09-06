from typing import Iterable
import logging

import pandas as pd

from .classification import *


class StatisticalClassification(Enum):
    FALSE_NEGATIVE = "fn"
    FALSE_POSITIVE = "fp"
    TRUE_POSITIVE = "tp"
    TRUE_NEGATIVE = "tn"


class Calculator:
    def create_gt_conf_column_from(self, probe_header: str) -> None:
        self.report["gt_conf"] = self.report[probe_header].apply(
            lambda column_name: ProbeHeader.from_string(column_name).gt_conf
        )

    def __init__(self, reports: Iterable[pd.DataFrame]):
        self.report = pd.concat(reports)


class RecallCalculator(Calculator):
    def __init__(self, reports: Iterable[pd.DataFrame]):
        super().__init__(reports)
        self.create_gt_conf_column_from("ref_probe_header")

    @staticmethod
    def statistical_classification(classification: str) -> StatisticalClassification:
        return {
            AlignmentAssessment.UNMAPPED: StatisticalClassification.FALSE_NEGATIVE,
            AlignmentAssessment.PARTIALLY_MAPPED: StatisticalClassification.FALSE_NEGATIVE,
            AlignmentAssessment.PRIMARY_CORRECT: StatisticalClassification.TRUE_POSITIVE,
            AlignmentAssessment.PRIMARY_INCORRECT: StatisticalClassification.FALSE_POSITIVE,
            AlignmentAssessment.SECONDARY_CORRECT: StatisticalClassification.TRUE_POSITIVE,
            AlignmentAssessment.SECONDARY_INCORRECT: StatisticalClassification.FALSE_POSITIVE,
            AlignmentAssessment.SUPPLEMENTARY_CORRECT: StatisticalClassification.TRUE_POSITIVE,
            AlignmentAssessment.SUPPLEMENTARY_INCORRECT: StatisticalClassification.FALSE_POSITIVE,
        }[AlignmentAssessment(classification)]

    def calculate_recall(self, conf_threshold: float = 0) -> float:
        confident_classifications = self.report.query(
            "gt_conf >= @conf_threshold"
        ).classification.to_list()
        counter = Counter(
            [
                RecallCalculator.statistical_classification(classification)
                for classification in confident_classifications
            ]
        )
        true_positives = counter[StatisticalClassification.TRUE_POSITIVE]

        num_unconfident_classifications = len(self.report) - len(
            confident_classifications
        )
        false_negatives = (
            counter[StatisticalClassification.FALSE_NEGATIVE]
            + num_unconfident_classifications
        )

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
