from typing import Iterable, Type, List
import logging

import pandas as pd

from .classification import *

from pathlib import Path


class StatisticalClassification(Enum):
    FALSE_NEGATIVE = "fn"
    FALSE_POSITIVE = "fp"
    TRUE_POSITIVE = "tp"
    TRUE_NEGATIVE = "tn"


class EmptyReportError(Exception):
    pass


class Calculator:
    def get_confident_classifications(
        self, conf_threshold: float
    ) -> List[str or float]:
        confident_classifications = self.report.query(
            "gt_conf >= @conf_threshold"
        ).classification.to_list()
        return confident_classifications

    def create_gt_conf_column_from(self, probe_header: str) -> None:
        self.report["gt_conf"] = self.report[probe_header].apply(
            lambda column_name: ProbeHeader.from_string(column_name).gt_conf
        )

    def get_maximum_gt_conf(self) -> float:
        return self.report["gt_conf"].max()

    def get_minimum_gt_conf(self) -> float:
        return self.report["gt_conf"].min()

    def __init__(self, reports: Iterable[pd.DataFrame]):
        self.report = pd.concat(reports)

    @classmethod
    def from_files(cls, paths: List[Path]) -> Type["Calculator"]:
        reports = [pd.read_csv(path, sep="\t", keep_default_na=False) for path in paths]
        return cls(reports)

    def __eq__(self, other: "Calculator"):
        return self.report.equals(other.report)


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
        confident_classifications = self.get_confident_classifications(conf_threshold)
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
            raise EmptyReportError(
                "There are not classifications to compute recall on (no true_positives or false_negatives)"
            )


class PrecisionCalculator(Calculator):
    def __init__(self, reports: Iterable[pd.DataFrame]):
        super().__init__(reports)
        self.create_gt_conf_column_from("query_probe_header")

    def calculate_precision(self, conf_threshold: float = 0.0) -> float:
        confident_classifications = self.get_confident_classifications(conf_threshold)
        true_positives = sum(confident_classifications)
        number_of_positives = len(confident_classifications)
        try:
            return true_positives / number_of_positives
        except ZeroDivisionError:
            raise EmptyReportError(
                "There are not classifications to compute precision on"
            )
