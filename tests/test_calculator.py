import pandas as pd
from evaluate.calculator import (
    RecallCalculator,
    PrecisionCalculator,
    EmptyReportError,
)
import pytest
from unittest.mock import patch, Mock
from evaluate.report import (
    Report,
    PrecisionReport,
    RecallReport
)
from tests.common import create_precision_report_row
from io import StringIO

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
    @patch.object(Report, Report.get_classifications_as_list.__name__, return_value=[
        "unmapped", "partially_mapped", "primary_correct", "primary_incorrect",
        "secondary_correct", "secondary_incorrect", "supplementary_correct",
        "supplementary_incorrect"
    ])
    @patch.object(RecallReport, RecallReport._create_helper_columns.__name__)
    @patch.object(RecallReport, RecallReport.assure_there_are_no_duplicated_evaluation.__name__)
    @patch.object(RecallReport, RecallReport.get_number_of_truth_probes.__name__, return_value=8)
    def test____calculate_info_wrt_truth_probes___one_classification_of_each(self, *mocks):
        report = RecallReport([pd.DataFrame()], False)
        true_positives, number_of_truth_probes = RecallCalculator._calculate_info_wrt_truth_probes(report)
        assert true_positives==3 and number_of_truth_probes==8

    @patch.object(Report, Report.get_classifications_as_list.__name__, return_value=[
        "unmapped", "partially_mapped", "primary_correct", "primary_incorrect",
        "secondary_correct", "secondary_incorrect", "supplementary_correct",
        "supplementary_incorrect", "partially_mapped", "partially_mapped",
        "primary_correct", "primary_correct", "primary_correct",
        "supplementary_incorrect", "supplementary_incorrect", "supplementary_incorrect",
        "unmapped", "unmapped", "unmapped",
    ])
    @patch.object(RecallReport, RecallReport._create_helper_columns.__name__)
    @patch.object(RecallReport, RecallReport.assure_there_are_no_duplicated_evaluation.__name__)
    @patch.object(RecallReport, RecallReport.get_number_of_truth_probes.__name__, return_value=19)
    def test____calculate_info_wrt_truth_probes___some_duplicated_classifications(self, *mocks):
        report = RecallReport([pd.DataFrame()], False)
        true_positives, number_of_truth_probes = RecallCalculator._calculate_info_wrt_truth_probes(report)
        assert true_positives == 6 and number_of_truth_probes == 19


    @patch.object(RecallReport, RecallReport.get_proportion_of_allele_seqs_found_for_each_variant.__name__,
                  return_value=[1.0, 0.5, 0.8, 1.0, 0.9, 1.0, 0.0, 0.1, 1.0])
    @patch.object(RecallReport, RecallReport.get_proportion_of_alleles_found_for_each_variant.__name__,
                  return_value=[0.0, 0.1, 0.2, 0.3, 1.0, 0.9, 0.8, 0.7, 0.6])
    @patch.object(RecallReport, RecallReport.get_number_of_variants.__name__, return_value=20)
    @patch.object(RecallReport, RecallReport._create_helper_columns.__name__)
    @patch.object(RecallReport, RecallReport.assure_there_are_no_duplicated_evaluation.__name__)
    def test____calculate_info_wrt_variants(self, *mocks):
        report = RecallReport([pd.DataFrame()], False)
        nb_variants_where_all_allele_seqs_were_found, nb_variants_found_wrt_alleles, variants_total = \
            RecallCalculator._calculate_info_wrt_variants(report)
        assert nb_variants_where_all_allele_seqs_were_found == 4 and \
               nb_variants_found_wrt_alleles == 4.6 and \
               variants_total == 20

    @patch.object(RecallReport, RecallReport._create_helper_columns.__name__)
    @patch.object(RecallReport, RecallReport.assure_there_are_no_duplicated_evaluation.__name__)
    @patch.object(Report, Report.get_report_satisfying_confidence_threshold.__name__)
    @patch.object(RecallCalculator, RecallCalculator._calculate_info_wrt_truth_probes.__name__, return_value=(5, 10))
    @patch.object(RecallCalculator, RecallCalculator._calculate_info_wrt_variants.__name__, return_value=(4, 8, 10))
    def test____calculate_recall_for_a_given_confidence(self, calculate_info_wrt_variants_mock,
                                                        calculate_info_wrt_truth_probes_mock,
                                                        get_report_satisfying_confidence_threshold_mock,
                                                        *other_mocks):
        # setup
        report_satisfying_confidence_threshold_mock = Mock()
        get_report_satisfying_confidence_threshold_mock.return_value = report_satisfying_confidence_threshold_mock
        report = RecallReport([pd.DataFrame()], False)
        calculator = RecallCalculator(report)

        recall_info_actual = calculator._calculate_recall_for_a_given_confidence(100)

        get_report_satisfying_confidence_threshold_mock.assert_called_once_with(100)
        calculate_info_wrt_truth_probes_mock.assert_called_once_with(report_satisfying_confidence_threshold_mock)
        calculate_info_wrt_variants_mock.assert_called_once_with(report_satisfying_confidence_threshold_mock)

        assert recall_info_actual.truth_probes_true_positives == 5
        assert recall_info_actual.truth_probes_total == 10
        assert recall_info_actual.nb_variants_where_all_allele_seqs_were_found == 4
        assert recall_info_actual.nb_variants_found_wrt_alleles == 8
        assert recall_info_actual.variants_total == 10
        assert recall_info_actual.recall_wrt_truth_probes == 0.5
        assert recall_info_actual.recall_wrt_variants_where_all_allele_seqs_were_found == 0.4
        assert recall_info_actual.recall_wrt_variants_found_wrt_alleles == 0.8


    @patch.object(RecallReport, RecallReport.get_proportion_of_allele_seqs_found_for_each_variant_with_nb_of_samples.__name__,
                  return_value=
    pd.read_csv(StringIO(
        """PANGENOME_VARIATION_ID,proportion_of_allele_seqs_found,NB_OF_SAMPLES
        0,1.0,3
        1,0.5,5
        2,0.0,7
        3,1.0,5
        4,0.0,5
        5,1.0,3
        """
    ), index_col="PANGENOME_VARIATION_ID"))
    @patch.object(RecallReport, RecallReport._create_helper_columns.__name__)
    @patch.object(RecallReport, RecallReport.assure_there_are_no_duplicated_evaluation.__name__)
    def test___get_recall_vs_nb_of_samples_report(self, *mocks):
        report = RecallReport([pd.DataFrame()], False)
        calculator = RecallCalculator(report)
        actual = calculator.get_recall_vs_nb_of_samples_report(list(range(2, 8)))
        expected = pd.read_csv(StringIO(
            """NB_OF_SAMPLES,recall
            2,0.0
            3,1.0
            4,0.0
            5,0.5
            6,0.0
            7,0.0
            """
        ))

        assert actual.equals(expected)


    @patch.object(RecallReport, RecallReport.get_proportion_of_allele_seqs_found_for_each_variant_with_nb_of_samples.__name__,
                  return_value=
    pd.read_csv(StringIO(
        """PANGENOME_VARIATION_ID,proportion_of_allele_seqs_found,NB_OF_SAMPLES
        0,1.0,3
        1,0.5,5
        2,0.0,7
        3,1.0,5
        4,0.0,5
        5,1.0,3
        """
    ), index_col="PANGENOME_VARIATION_ID"))
    @patch.object(RecallReport, RecallReport._create_helper_columns.__name__)
    @patch.object(RecallReport, RecallReport.assure_there_are_no_duplicated_evaluation.__name__)
    def test___get_recall_vs_nb_of_samples_report___return_only_the_samples_given_in_parameter(self, *mocks):
        report = RecallReport([pd.DataFrame()], False)
        calculator = RecallCalculator(report)
        actual = calculator.get_recall_vs_nb_of_samples_report([2, 5])
        expected = pd.read_csv(StringIO(
            """NB_OF_SAMPLES,recall
            2,0.0
            5,0.5
            """
        ))

        assert actual.equals(expected)