from .report import Report, PrecisionReport, RecallReport
from collections import Counter
from enum import Enum
from .classification import AlignmentAssessment
import numpy as np
import pandas as pd


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

    def _get_all_genotype_points(self, number_of_datapoints):
        min_gt = self.report.get_minimum_gt_conf()
        max_gt = self.report.get_maximum_gt_conf()
        all_gts = np.linspace(min_gt, max_gt, number_of_datapoints)
        return all_gts


class PrecisionCalculator(Calculator):
    def __init__(self, report: PrecisionReport):
        super().__init__(report)

    def get_precision_report(self, number_of_datapoints) -> pd.DataFrame:
        all_gts = self._get_all_genotype_points(number_of_datapoints)

        gts = []
        precisions = []
        error_rates = []
        nb_of_correct_calls = []
        nb_of_total_calls = []
        for gt in all_gts:
            try:
                precision_info = self._calculate_precision_for_a_given_confidence(gt)
                gts.append(gt)
                precisions.append(precision_info.precision)
                error_rates.append(1 - precision_info.precision)
                nb_of_correct_calls.append(precision_info.true_positives)
                nb_of_total_calls.append(precision_info.total)
            except EmptyReportError:
                pass

        precision_df = pd.DataFrame(
            data={
                "GT": gts,
                "step_GT": list(range(len(gts))),
                "precision": precisions,
                "error_rate": error_rates,
                "nb_of_correct_calls": nb_of_correct_calls,
                "nb_of_total_calls": nb_of_total_calls
            }
        )

        return precision_df

    def _calculate_precision_for_a_given_confidence(self, conf_threshold: float = 0.0) -> PrecisionInfo:
        confident_classifications = self.report.get_confident_classifications(conf_threshold)
        true_positives = sum(confident_classifications)
        number_of_calls = len(confident_classifications)
        return PrecisionInfo(true_positives, number_of_calls)


class RecallCalculator(Calculator):
    def __init__(self, report: RecallReport):
        super().__init__(report)

    def get_recall_report(self, number_of_datapoints) -> pd.DataFrame:
        all_gts = self._get_all_genotype_points(number_of_datapoints)

        gts = []
        recalls = []
        nb_of_truth_probes_found = []
        nb_of_truth_probes_in_total = []
        for gt in all_gts:
            try:
                recall_info = self._calculate_recall_for_a_given_confidence(gt)
                gts.append(gt)
                recalls.append(recall_info.recall)
                nb_of_truth_probes_found.append(recall_info.true_positives)
                nb_of_truth_probes_in_total.append(recall_info.total)
            except EmptyReportError:
                pass

        recall_df = pd.DataFrame(
            data={
                "GT": gts,
                "step_GT": list(range(len(gts))),
                "recall": recalls,
                "nb_of_truth_probes_found": nb_of_truth_probes_found,
                "nb_of_truth_probes_in_total": nb_of_truth_probes_in_total
            }
        )

        return recall_df

    def _calculate_recall_for_a_given_confidence(self, conf_threshold: float = 0) -> RecallInfo:
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
