from .report import Report, PrecisionReport, RecallReport
from collections import Counter
from enum import Enum
from .classification import AlignmentAssessment


class EmptyReportError(Exception):
    pass
class CalculatorInfo:
    def __init__(self, true_positives: float, total: float):
        self.true_positives = float(true_positives)
        self.total = float(total)

        try:
            self.ratio = true_positives / total
        except ZeroDivisionError:
            raise EmptyReportError(
                "There are not classifications to compute recall/precision on."
            )

class PrecisionInfo(CalculatorInfo):
    def __init__(self, true_positives: float, number_of_calls: float):
        super().__init__(true_positives, number_of_calls)
        self.precision = self.ratio
class RecallInfo(CalculatorInfo):
    def __init__(self, true_positives: float, number_of_truth_probes: float):
        super().__init__(true_positives, number_of_truth_probes)
        self.recall = self.ratio



class Calculator:
    def __init__(self, report: Report):
        self.report = report

class PrecisionCalculator(Calculator):
    def __init__(self, report: PrecisionReport):
        super().__init__(report)

    def calculate_precision(self, conf_threshold: float = 0.0) -> PrecisionInfo:
        confident_classifications = self.report.get_confident_classifications(conf_threshold)
        true_positives = sum(confident_classifications)
        number_of_calls = len(confident_classifications)
        return PrecisionInfo(true_positives, number_of_calls)

class RecallCalculator(Calculator):
    def __init__(self, report: RecallReport):
        super().__init__(report)

    def calculate_recall(self, conf_threshold: float = 0) -> RecallInfo:
        confident_classifications = self.report.get_confident_classifications(conf_threshold)
        counter = Counter(
            [
                RecallCalculator._statistical_classification(classification)
                for classification in confident_classifications
            ]
        )
        true_positives = counter[RecallCalculator.StatisticalClassification.TRUE_POSITIVE]
        return RecallInfo(true_positives, self.report.get_number_of_truth_probes())

    class StatisticalClassification(Enum):
        FALSE_NEGATIVE = "fn"
        FALSE_POSITIVE = "fp"
        TRUE_POSITIVE = "tp"
        TRUE_NEGATIVE = "tn"

    @staticmethod
    def _statistical_classification(classification: str) -> StatisticalClassification:
        return {
            AlignmentAssessment.UNMAPPED: RecallCalculator.StatisticalClassification.FALSE_NEGATIVE,
            AlignmentAssessment.PARTIALLY_MAPPED: RecallCalculator.StatisticalClassification.FALSE_NEGATIVE,
            AlignmentAssessment.PRIMARY_CORRECT: RecallCalculator.StatisticalClassification.TRUE_POSITIVE,
            AlignmentAssessment.PRIMARY_INCORRECT: RecallCalculator.StatisticalClassification.FALSE_POSITIVE,
            AlignmentAssessment.SECONDARY_CORRECT: RecallCalculator.StatisticalClassification.TRUE_POSITIVE,
            AlignmentAssessment.SECONDARY_INCORRECT: RecallCalculator.StatisticalClassification.FALSE_POSITIVE,
            AlignmentAssessment.SUPPLEMENTARY_CORRECT: RecallCalculator.StatisticalClassification.TRUE_POSITIVE,
            AlignmentAssessment.SUPPLEMENTARY_INCORRECT: RecallCalculator.StatisticalClassification.FALSE_POSITIVE,
        }[AlignmentAssessment(classification)]
