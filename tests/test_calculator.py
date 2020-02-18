import pandas as pd
from evaluate.calculator import (
    Calculator,
    RecallCalculator,
    PrecisionCalculator,
    EmptyReportError,
)
from evaluate.classification import AlignmentAssessment
import pytest
from unittest.mock import patch
from evaluate.report import (
    Report,
    PrecisionReport,
    RecallReport
)
from numpy.testing import assert_allclose
from tests.common import create_recall_report_row, create_precision_report_row

class TestCalculator:
    @patch.object(Report, Report.get_minimum_gt_conf.__name__, return_value = 2.5)
    @patch.object(Report, Report.get_maximum_gt_conf.__name__, return_value = 33.3)
    def test____get_all_genotype_points___no_datapoints(self, *patches):
        report_mock = Report([pd.DataFrame()])
        calculator = Calculator(report_mock)
        actual = calculator._get_all_genotype_points(0)
        expected = []
        assert_allclose(actual, expected)

    @patch.object(Report, Report.get_minimum_gt_conf.__name__, return_value = 2.5)
    @patch.object(Report, Report.get_maximum_gt_conf.__name__, return_value = 33.3)
    def test____get_all_genotype_points___one_datapoint(self, *patches):
        report_mock = Report([pd.DataFrame()])
        calculator = Calculator(report_mock)
        actual = calculator._get_all_genotype_points(1)
        expected = [2.5]
        assert_allclose(actual, expected)

    @patch.object(Report, Report.get_minimum_gt_conf.__name__, return_value = 2.5)
    @patch.object(Report, Report.get_maximum_gt_conf.__name__, return_value = 33.3)
    def test____get_all_genotype_points___two_datapoints(self, *patches):
        report_mock = Report([pd.DataFrame()])
        calculator = Calculator(report_mock)
        actual = calculator._get_all_genotype_points(2)
        expected = [2.5, 33.3]
        assert_allclose(actual, expected)

    @patch.object(Report, Report.get_minimum_gt_conf.__name__, return_value = 2.5)
    @patch.object(Report, Report.get_maximum_gt_conf.__name__, return_value = 33.3)
    def test____get_all_genotype_points___three_datapoints(self, *patches):
        report_mock = Report([pd.DataFrame()])
        calculator = Calculator(report_mock)
        actual = calculator._get_all_genotype_points(3)
        expected = [2.5, 17.9 ,33.3]
        assert_allclose(actual, expected)

    @patch.object(Report, Report.get_minimum_gt_conf.__name__, return_value = 2.5)
    @patch.object(Report, Report.get_maximum_gt_conf.__name__, return_value = 33.3)
    def test____get_all_genotype_points___five_datapoints(self, *patches):
        report_mock = Report([pd.DataFrame()])
        calculator = Calculator(report_mock)
        actual = calculator._get_all_genotype_points(5)
        expected = [2.5, 10.2, 17.9, 25.6, 33.3]
        assert_allclose(actual, expected)


class TestPrecisionCalculator:
    def test_calculatePrecision_NoReportsRaisesEmptyReportError(self):
        columns = ["sample", "query_probe_header", "ref_probe_header", "classification"]
        df = pd.DataFrame(columns=columns)
        report = PrecisionReport([df])
        calculator = PrecisionCalculator(report)

        with pytest.raises(EmptyReportError):
            calculator._calculate_precision_for_a_given_confidence()

    def test_calculatePrecision_OneReportWithOneRowCompletelyCorrectReturnsOne(self):
        columns = ["sample", "query_probe_header", "ref_probe_header", "classification"]
        df = pd.DataFrame(
            data=[create_precision_report_row(1.0, gt_conf=100)], columns=columns
        )
        report = PrecisionReport([df])
        calculator = PrecisionCalculator(report)

        actual = calculator._calculate_precision_for_a_given_confidence()

        assert actual.precision == 1.0
        assert actual.true_positives == 1.0
        assert actual.total == 1.0

    def test_calculatePrecision_OneReportWithOneRowCompletelyIncorrectReturnsZero(self):
        columns = ["sample", "query_probe_header", "ref_probe_header", "classification"]
        df = pd.DataFrame(
            data=[create_precision_report_row(0.0, gt_conf=100)], columns=columns
        )
        report = PrecisionReport([df])
        calculator = PrecisionCalculator(report)

        actual = calculator._calculate_precision_for_a_given_confidence()

        assert actual.precision == 0.0
        assert actual.true_positives == 0.0
        assert actual.total == 1.0

    def test_calculatePrecision_OneReportWithOneRowCompletelyCorrectBelowConfThreasholdRaisesEmptyReportError(
        self
    ):
        columns = ["sample", "query_probe_header", "ref_probe_header", "classification"]
        df = pd.DataFrame(
            data=[create_precision_report_row(1.0, gt_conf=10)], columns=columns
        )
        report = PrecisionReport([df])
        calculator = PrecisionCalculator(report)

        confidence_threshold = 60

        with pytest.raises(EmptyReportError):
            calculator._calculate_precision_for_a_given_confidence(confidence_threshold)

    def test_calculatePrecision_OneReportWithOneRowCompletelyCorrectEqualConfThreasholdReturnsOne(
        self
    ):
        columns = ["sample", "query_probe_header", "ref_probe_header", "classification"]
        df = pd.DataFrame(
            data=[create_precision_report_row(1.0, gt_conf=60)], columns=columns
        )
        report = PrecisionReport([df])
        calculator = PrecisionCalculator(report)

        confidence_threshold = 60

        actual = calculator._calculate_precision_for_a_given_confidence(confidence_threshold)

        assert actual.precision == 1.0
        assert actual.true_positives == 1.0
        assert actual.total == 1.0

    def test_calculatePrecision_OneReportWithTwoRowsPartiallyCorrect(self):
        columns = ["sample", "query_probe_header", "ref_probe_header", "classification"]
        df = pd.DataFrame(
            data=[
                create_precision_report_row(0.5, gt_conf=100),
                create_precision_report_row(0.7, gt_conf=100),
            ],
            columns=columns,
        )
        report = PrecisionReport([df])
        calculator = PrecisionCalculator(report)

        actual = calculator._calculate_precision_for_a_given_confidence()

        assert actual.precision == 1.2/2
        assert actual.true_positives == 1.2
        assert actual.total == 2.0


    def test_calculatePrecision_OneReportWithThreeRowsTwoPartiallyCorrectOneBelowThreshold(
        self
    ):
        columns = ["sample", "query_probe_header", "ref_probe_header", "classification"]
        df = pd.DataFrame(
            data=[
                create_precision_report_row(0.4, gt_conf=100),
                create_precision_report_row(0.8, gt_conf=20),
                create_precision_report_row(0.3, gt_conf=100),
            ],
            columns=columns,
        )
        report = PrecisionReport([df])
        calculator = PrecisionCalculator(report)

        confidence_threshold = 80

        actual = calculator._calculate_precision_for_a_given_confidence(confidence_threshold)

        assert actual.precision == 0.7/2.0
        assert actual.true_positives == 0.7
        assert actual.total == 2.0





class TestRecallCalculator:
    def test_calculateRecall_noReportsRaisesEmptyReportError(self):
        columns = ["sample", "query_probe_header", "ref_probe_header", "classification"]
        df = pd.DataFrame(columns=columns)
        report = RecallReport([df])
        calculator = RecallCalculator(report)
        threshold = 0

        with pytest.raises(EmptyReportError):
            calculator._calculate_recall_for_a_given_confidence(conf_threshold=threshold)

    def test_calculateRecall_oneReportNoTruePositivesReturnsZero(self):
        columns = ["sample", "query_probe_header", "ref_probe_header", "classification"]
        df = pd.DataFrame(
            data=[
                create_recall_report_row("truth_probe_1", AlignmentAssessment.UNMAPPED, gt_conf=100),
                create_recall_report_row("truth_probe_2", AlignmentAssessment.UNMAPPED, gt_conf=100),
                create_recall_report_row(
                    "truth_probe_3", AlignmentAssessment.PRIMARY_CORRECT, gt_conf=10
                ),
                create_recall_report_row(
                    "truth_probe_4", AlignmentAssessment.PRIMARY_INCORRECT, gt_conf=100
                ),
            ],
            columns=columns,
        )
        report = RecallReport([df])
        calculator = RecallCalculator(report)

        threshold = 60

        actual = calculator._calculate_recall_for_a_given_confidence(conf_threshold=threshold)

        assert actual.recall == 0.0
        assert actual.true_positives == 0.0
        assert actual.total == 4.0

    def test_calculateRecall_oneReportNoTruePositivesTwoTruthProbesReturnsZero(self):
        columns = ["sample", "query_probe_header", "ref_probe_header", "classification"]
        df = pd.DataFrame(
            data=[
                create_recall_report_row("truth_probe_1", AlignmentAssessment.UNMAPPED, gt_conf=100),
                create_recall_report_row("truth_probe_2", AlignmentAssessment.UNMAPPED, gt_conf=100),
                create_recall_report_row(
                    "truth_probe_1", AlignmentAssessment.PRIMARY_CORRECT, gt_conf=10
                ),
                create_recall_report_row(
                    "truth_probe_2", AlignmentAssessment.PRIMARY_INCORRECT, gt_conf=100
                ),
                create_recall_report_row(
                    "truth_probe_1", AlignmentAssessment.PARTIALLY_MAPPED, gt_conf=10
                ),
                create_recall_report_row(
                    "truth_probe_2", AlignmentAssessment.SECONDARY_INCORRECT, gt_conf=100
                ),
            ],
            columns=columns,
        )
        report = RecallReport([df])
        calculator = RecallCalculator(report)

        threshold = 60

        actual = calculator._calculate_recall_for_a_given_confidence(conf_threshold=threshold)

        assert actual.recall == 0.0
        assert actual.true_positives == 0.0
        assert actual.total == 2.0

    def test_calculateRecall_oneReportNoFalseNegativesReturnsOne(self):
        columns = ["sample", "query_probe_header", "ref_probe_header", "classification"]
        df = pd.DataFrame(
            data=[
                create_recall_report_row(
                    "truth_probe_1", AlignmentAssessment.PRIMARY_CORRECT, gt_conf=100
                ),
                create_recall_report_row(
                    "truth_probe_2", AlignmentAssessment.SECONDARY_CORRECT, gt_conf=100
                ),
                create_recall_report_row(
                    "truth_probe_3", AlignmentAssessment.SUPPLEMENTARY_CORRECT, gt_conf=100
                ),
                create_recall_report_row(
                    "truth_probe_4", AlignmentAssessment.SUPPLEMENTARY_CORRECT, gt_conf=100
                ),
            ],
            columns=columns,
        )
        report = RecallReport([df])
        calculator = RecallCalculator(report)

        threshold = 60

        actual = calculator._calculate_recall_for_a_given_confidence(conf_threshold=threshold)

        assert actual.recall == 1.0
        assert actual.true_positives == 4.0
        assert actual.total == 4.0


    def test_calculateRecall_oneReportNoFalseNegativesTwoProbesReturnsOne(self):
        columns = ["sample", "query_probe_header", "ref_probe_header", "classification"]
        df = pd.DataFrame(
            data=[
                create_recall_report_row(
                    "truth_probe_1", AlignmentAssessment.PRIMARY_CORRECT, gt_conf=100
                ),
                create_recall_report_row(
                    "truth_probe_2", AlignmentAssessment.SECONDARY_CORRECT, gt_conf=100
                ),
                create_recall_report_row(
                    "truth_probe_1", AlignmentAssessment.SUPPLEMENTARY_CORRECT, gt_conf=100
                ),
                create_recall_report_row(
                    "truth_probe_2", AlignmentAssessment.PRIMARY_CORRECT, gt_conf=100
                ),
                create_recall_report_row(
                    "truth_probe_1", AlignmentAssessment.SUPPLEMENTARY_INCORRECT, gt_conf=100
                ),
                create_recall_report_row(
                    "truth_probe_2", AlignmentAssessment.PRIMARY_INCORRECT, gt_conf=100
                ),
                create_recall_report_row(
                    "truth_probe_1", AlignmentAssessment.UNMAPPED, gt_conf=100
                ),
                create_recall_report_row(
                    "truth_probe_2", AlignmentAssessment.PARTIALLY_MAPPED, gt_conf=100
                ),
            ],
            columns=columns,
        )
        report = RecallReport([df])
        calculator = RecallCalculator(report)

        threshold = 60

        actual = calculator._calculate_recall_for_a_given_confidence(conf_threshold=threshold)

        assert actual.recall == 1.0
        assert actual.true_positives == 2.0
        assert actual.total == 2.0


    def test_calculateRecall_oneReportHalfTruePositiveHalfFalseNegativeReturnsFifty(
        self
    ):
        columns = ["sample", "query_probe_header", "ref_probe_header", "classification"]
        df = pd.DataFrame(
            data=[
                create_recall_report_row(
                    "truth_probe_1", AlignmentAssessment.PRIMARY_CORRECT, gt_conf=10
                ),
                create_recall_report_row(
                    "truth_probe_1", AlignmentAssessment.SUPPLEMENTARY_CORRECT, gt_conf=100
                ),
                create_recall_report_row(
                    "truth_probe_2", AlignmentAssessment.PRIMARY_CORRECT, gt_conf=20
                ),
                create_recall_report_row(
                    "truth_probe_2", AlignmentAssessment.SECONDARY_CORRECT, gt_conf=100
                ),
                create_recall_report_row(
                    "truth_probe_3", AlignmentAssessment.PRIMARY_CORRECT, gt_conf=100
                ),
                create_recall_report_row(
                    "truth_probe_4", AlignmentAssessment.UNMAPPED, gt_conf=100
                ),
                create_recall_report_row(
                    "truth_probe_5", AlignmentAssessment.SUPPLEMENTARY_INCORRECT, gt_conf=10
                ),
                create_recall_report_row(
                    "truth_probe_6", AlignmentAssessment.PRIMARY_INCORRECT, gt_conf=100
                ),
            ],
            columns=columns,
        )
        report = RecallReport([df])
        calculator = RecallCalculator(report)

        threshold = 60

        actual = calculator._calculate_recall_for_a_given_confidence(conf_threshold=threshold)

        assert actual.recall == 0.5
        assert actual.true_positives == 3.0
        assert actual.total == 6.0


    def test_calculateRecall_oneReportAllTruePositivesAllBelowThresholdReturnsZero(
        self
    ):
        columns = ["sample", "query_probe_header", "ref_probe_header", "classification"]
        df = pd.DataFrame(
            data=[
                create_recall_report_row(
                    "truth_probe_1", AlignmentAssessment.PRIMARY_CORRECT, gt_conf=50
                ),
                create_recall_report_row(
                    "truth_probe_2", AlignmentAssessment.PRIMARY_CORRECT, gt_conf=80
                ),
            ],
            columns=columns,
        )
        report = RecallReport([df])
        calculator = RecallCalculator(report)

        threshold = 100

        actual = calculator._calculate_recall_for_a_given_confidence(conf_threshold=threshold)

        assert actual.recall == 0.0
        assert actual.true_positives == 0.0
        assert actual.total == 2.0

    def test_calculateRecall_oneReportAllTruePositivesAllBelowThresholdOneTruthProbeReturnsZero(
        self
    ):
        columns = ["sample", "query_probe_header", "ref_probe_header", "classification"]
        df = pd.DataFrame(
            data=[
                create_recall_report_row(
                    "truth_probe_1", AlignmentAssessment.PRIMARY_CORRECT, gt_conf=50
                ),
                create_recall_report_row(
                    "truth_probe_1", AlignmentAssessment.PRIMARY_CORRECT, gt_conf=80
                ),
            ],
            columns=columns,
        )
        report = RecallReport([df])
        calculator = RecallCalculator(report)

        threshold = 100

        actual = calculator._calculate_recall_for_a_given_confidence(conf_threshold=threshold)

        assert actual.recall == 0.0
        assert actual.true_positives == 0.0
        assert actual.total == 1.0


    def test_calculateRecall_AllAlignmentsAreIncorrectAndOneIsCorrect(
        self
    ):
        columns = ["sample", "query_probe_header", "ref_probe_header", "classification"]
        df = pd.DataFrame(
            data=[
                create_recall_report_row(
                    "truth_probe_1", AlignmentAssessment.PRIMARY_INCORRECT, gt_conf=100
                ),
                create_recall_report_row(
                    "truth_probe_2", AlignmentAssessment.PRIMARY_INCORRECT, gt_conf=100
                ),
                create_recall_report_row(
                    "truth_probe_3", AlignmentAssessment.PRIMARY_INCORRECT, gt_conf=100
                ),
                create_recall_report_row(
                    "truth_probe_4", AlignmentAssessment.PRIMARY_CORRECT, gt_conf=100
                ),
            ],
            columns=columns,
        )
        report = RecallReport([df])
        calculator = RecallCalculator(report)

        threshold = 60

        actual = calculator._calculate_recall_for_a_given_confidence(conf_threshold=threshold)

        assert actual.recall == 0.25
        assert actual.true_positives == 1.0
        assert actual.total == 4.0


    def test_calculateRecall_AllAlignmentsAreIncorrectAndOneIsCorrectTwoTruthProbes(
        self
    ):
        columns = ["sample", "query_probe_header", "ref_probe_header", "classification"]
        df = pd.DataFrame(
            data=[
                create_recall_report_row(
                    "truth_probe_1", AlignmentAssessment.PRIMARY_INCORRECT, gt_conf=100
                ),
                create_recall_report_row(
                    "truth_probe_2", AlignmentAssessment.PRIMARY_INCORRECT, gt_conf=100
                ),
                create_recall_report_row(
                    "truth_probe_1", AlignmentAssessment.PRIMARY_INCORRECT, gt_conf=100
                ),
                create_recall_report_row(
                    "truth_probe_2", AlignmentAssessment.PRIMARY_CORRECT, gt_conf=100
                ),
            ],
            columns=columns,
        )
        report = RecallReport([df])
        calculator = RecallCalculator(report)

        threshold = 60

        actual = calculator._calculate_recall_for_a_given_confidence(conf_threshold=threshold)

        assert actual.recall == 0.5
        assert actual.true_positives == 1.0
        assert actual.total == 2.0


    def test_calculateRecall_twoReportsHalfTruePositiveHalfFalseNegativeReturnsFifty(
        self
    ):
        columns = ["sample", "query_probe_header", "ref_probe_header", "classification"]
        report1 = pd.DataFrame(
            data=[
                create_recall_report_row(
                    "truth_probe_1", AlignmentAssessment.PRIMARY_CORRECT, gt_conf=10
                ),
                create_recall_report_row(
                    "truth_probe_2", AlignmentAssessment.PRIMARY_CORRECT, gt_conf=100
                ),
                create_recall_report_row(
                    "truth_probe_3", AlignmentAssessment.PRIMARY_CORRECT, gt_conf=100
                ),
                create_recall_report_row(
                    "truth_probe_4", AlignmentAssessment.UNMAPPED, gt_conf=100
                ),
                create_recall_report_row(
                    "truth_probe_5", AlignmentAssessment.SUPPLEMENTARY_INCORRECT, gt_conf=10
                ),
                create_recall_report_row(
                    "truth_probe_6", AlignmentAssessment.SECONDARY_CORRECT, gt_conf=100
                ),
                create_recall_report_row(
                    "truth_probe_7", AlignmentAssessment.PRIMARY_INCORRECT, gt_conf=100
                ),
            ],
            columns=columns,
        )
        report2 = pd.DataFrame(
            data=[
                create_recall_report_row(
                    "truth_probe_8", AlignmentAssessment.PRIMARY_CORRECT, gt_conf=100
                ),
                create_recall_report_row(
                    "truth_probe_9", AlignmentAssessment.PRIMARY_CORRECT, gt_conf=100
                ),
                create_recall_report_row(
                    "truth_probe_10", AlignmentAssessment.PRIMARY_INCORRECT, gt_conf=100
                ),
            ],
            columns=columns,
        )
        report = RecallReport([report1, report2])
        calculator = RecallCalculator(report)

        threshold = 60

        actual = calculator._calculate_recall_for_a_given_confidence(conf_threshold=threshold)

        assert actual.recall == 0.5
        assert actual.true_positives == 5.0
        assert actual.total == 10.0


