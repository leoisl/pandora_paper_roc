import pandas as pd

from evaluate.probe import ProbeHeader
from evaluate.recall import RecallCalculator, StatisticalClassification
from evaluate.classification import AlignmentAssessment


def create_report_row(
    classification: AlignmentAssessment, gt_conf: float = 0, sample: str = "sample1"
) -> pd.Series:
    truth_probe_header = ProbeHeader()
    vcf_probe_header = ProbeHeader(gt_conf=gt_conf)
    data = {
        "sample": sample,
        "truth_probe_header": str(truth_probe_header),
        "vcf_probe_header": str(vcf_probe_header),
        "classification": classification,
    }
    return pd.Series(data=data)


class TestRecallCalculator:
    def test_statisticalClassification_unmappedReturnsFalseNegative(self):
        classification = "unmapped"
        row = create_report_row(classification)

        actual = RecallCalculator.statistical_classification(row)
        expected = StatisticalClassification.FALSE_NEGATIVE

        assert actual == expected

    def test_statisticalClassification_partiallyMappedReturnsFalseNegative(self):
        classification = "partially_mapped"
        row = create_report_row(classification)

        actual = RecallCalculator.statistical_classification(row)
        expected = StatisticalClassification.FALSE_NEGATIVE

        assert actual == expected

    def test_statisticalClassification_incorrectAboveConfReturnsFalsePositive(self):
        classification = AlignmentAssessment.PRIMARY_INCORRECT
        gt_conf = 10
        row = create_report_row(classification, gt_conf)
        threshold = 5

        actual = RecallCalculator.statistical_classification(
            row, conf_threshold=threshold
        )
        expected = StatisticalClassification.FALSE_POSITIVE

        assert actual == expected

    def test_statisticalClassification_incorrectBelowConfReturnsFalseNegative(self):
        classification = AlignmentAssessment.PRIMARY_INCORRECT
        gt_conf = 10
        row = create_report_row(classification, gt_conf)
        threshold = 50

        actual = RecallCalculator.statistical_classification(
            row, conf_threshold=threshold
        )
        expected = StatisticalClassification.FALSE_NEGATIVE

        assert actual == expected

    def test_statisticalClassification_correctAboveConfReturnsTruePositive(self):
        classification = AlignmentAssessment.PRIMARY_CORRECT
        gt_conf = 10
        row = create_report_row(classification, gt_conf)
        threshold = 5

        actual = RecallCalculator.statistical_classification(
            row, conf_threshold=threshold
        )
        expected = StatisticalClassification.TRUE_POSITIVE

        assert actual == expected

    def test_statisticalClassification_correctBelowConfReturnsFalseNegative(self):
        classification = AlignmentAssessment.PRIMARY_CORRECT
        gt_conf = 10
        row = create_report_row(classification, gt_conf)
        threshold = 50

        actual = RecallCalculator.statistical_classification(
            row, conf_threshold=threshold
        )
        expected = StatisticalClassification.FALSE_NEGATIVE

        assert actual == expected

    def test_statisticalClassification_secondaryIncorrectAboveConfReturnsFalsePositive(
        self
    ):
        classification = AlignmentAssessment.SECONDARY_INCORRECT
        gt_conf = 10
        row = create_report_row(classification, gt_conf)
        threshold = 5

        actual = RecallCalculator.statistical_classification(
            row, conf_threshold=threshold
        )
        expected = StatisticalClassification.FALSE_POSITIVE

        assert actual == expected

    def test_statisticalClassification_secondaryIncorrectBelowConfReturnsFalseNegative(
        self
    ):
        classification = AlignmentAssessment.SECONDARY_INCORRECT
        gt_conf = 10
        row = create_report_row(classification, gt_conf)
        threshold = 50

        actual = RecallCalculator.statistical_classification(
            row, conf_threshold=threshold
        )
        expected = StatisticalClassification.FALSE_NEGATIVE

        assert actual == expected

    def test_statisticalClassification_secondaryCorrectAboveConfReturnsTruePositive(
        self
    ):
        classification = AlignmentAssessment.SECONDARY_CORRECT
        gt_conf = 10
        row = create_report_row(classification, gt_conf)
        threshold = 5

        actual = RecallCalculator.statistical_classification(
            row, conf_threshold=threshold
        )
        expected = StatisticalClassification.TRUE_POSITIVE

        assert actual == expected

    def test_statisticalClassification_secondaryCorrectBelowConfReturnsFalseNegative(
        self
    ):
        classification = AlignmentAssessment.SECONDARY_CORRECT
        gt_conf = 10
        row = create_report_row(classification, gt_conf)
        threshold = 50

        actual = RecallCalculator.statistical_classification(
            row, conf_threshold=threshold
        )
        expected = StatisticalClassification.FALSE_NEGATIVE

        assert actual == expected

    def test_statisticalClassification_supplementaryIncorrectAboveConfReturnsFalsePositive(
        self
    ):
        classification = AlignmentAssessment.SUPPLEMENTARY_INCORRECT
        gt_conf = 10
        row = create_report_row(classification, gt_conf)
        threshold = 5

        actual = RecallCalculator.statistical_classification(
            row, conf_threshold=threshold
        )
        expected = StatisticalClassification.FALSE_POSITIVE

        assert actual == expected

    def test_statisticalClassification_supplementaryIncorrectBelowConfReturnsFalseNegative(
        self
    ):
        classification = AlignmentAssessment.SUPPLEMENTARY_INCORRECT
        gt_conf = 10
        row = create_report_row(classification, gt_conf)
        threshold = 50

        actual = RecallCalculator.statistical_classification(
            row, conf_threshold=threshold
        )
        expected = StatisticalClassification.FALSE_NEGATIVE

        assert actual == expected

    def test_statisticalClassification_supplementaryCorrectAboveConfReturnsTruePositive(
        self
    ):
        classification = AlignmentAssessment.SUPPLEMENTARY_CORRECT
        gt_conf = 10
        row = create_report_row(classification, gt_conf)
        threshold = 5

        actual = RecallCalculator.statistical_classification(
            row, conf_threshold=threshold
        )
        expected = StatisticalClassification.TRUE_POSITIVE

        assert actual == expected

    def test_statisticalClassification_supplementaryCorrectBelowConfReturnsFalseNegative(
        self
    ):
        classification = AlignmentAssessment.SUPPLEMENTARY_CORRECT
        gt_conf = 10
        row = create_report_row(classification, gt_conf)
        threshold = 50

        actual = RecallCalculator.statistical_classification(
            row, conf_threshold=threshold
        )
        expected = StatisticalClassification.FALSE_NEGATIVE

        assert actual == expected

    def test_calculateRecall_noReportsReturnsZero(self):
        columns = ["sample", "truth_probe_header", "vcf_probe_header", "classification"]
        report = pd.DataFrame(columns=columns)
        calculator = RecallCalculator([report])
        threshold = 0

        actual = calculator.calculate_recall(conf_threshold=threshold)
        expected = 0

        assert actual == expected

    def test_calculateRecall_oneReportNoTruePositivesReturnsZero(self):
        columns = ["sample", "truth_probe_header", "vcf_probe_header", "classification"]
        report = pd.DataFrame(
            data=[
                create_report_row(AlignmentAssessment.UNMAPPED, gt_conf=100),
                create_report_row(AlignmentAssessment.UNMAPPED, gt_conf=100),
                create_report_row(AlignmentAssessment.PRIMARY_CORRECT, gt_conf=10),
                create_report_row(AlignmentAssessment.PRIMARY_INCORRECT, gt_conf=100),
            ],
            columns=columns,
        )
        calculator = RecallCalculator([report])
        threshold = 60

        actual = calculator.calculate_recall(conf_threshold=threshold)
        expected = 0

        assert actual == expected

    def test_calculateRecall_oneReportNoFalseNegativesReturnsOne(self):
        columns = ["sample", "truth_probe_header", "vcf_probe_header", "classification"]
        report = pd.DataFrame(
            data=[
                create_report_row(AlignmentAssessment.PRIMARY_CORRECT, gt_conf=100),
                create_report_row(AlignmentAssessment.PRIMARY_CORRECT, gt_conf=100),
                create_report_row(AlignmentAssessment.PRIMARY_INCORRECT, gt_conf=100),
            ],
            columns=columns,
        )
        calculator = RecallCalculator([report])
        threshold = 60

        actual = calculator.calculate_recall(conf_threshold=threshold)
        expected = 1

        assert actual == expected

    def test_calculateRecall_oneReportHalfTruePositiveHalfFalseNegativeReturnsFifty(
        self
    ):
        columns = ["sample", "truth_probe_header", "vcf_probe_header", "classification"]
        report = pd.DataFrame(
            data=[
                create_report_row(AlignmentAssessment.PRIMARY_CORRECT, gt_conf=10),
                create_report_row(AlignmentAssessment.PRIMARY_CORRECT, gt_conf=100),
                create_report_row(AlignmentAssessment.PRIMARY_CORRECT, gt_conf=100),
                create_report_row(AlignmentAssessment.UNMAPPED, gt_conf=100),
                create_report_row("supplementary_incorrect", gt_conf=10),
                create_report_row("secondary_correct", gt_conf=100),
                create_report_row(AlignmentAssessment.PRIMARY_INCORRECT, gt_conf=100),
            ],
            columns=columns,
        )
        calculator = RecallCalculator([report])
        threshold = 60

        actual = calculator.calculate_recall(conf_threshold=threshold)
        expected = 0.5

        assert actual == expected

    def test_calculateRecall_oneReportNoTruePositivesOrFalseNegativesReturnsZero(self):
        columns = ["sample", "truth_probe_header", "vcf_probe_header", "classification"]
        report = pd.DataFrame(
            data=[
                create_report_row("supplementary_incorrect", gt_conf=100),
                create_report_row(AlignmentAssessment.PRIMARY_INCORRECT, gt_conf=100),
            ],
            columns=columns,
        )
        calculator = RecallCalculator([report])
        threshold = 60

        actual = calculator.calculate_recall(conf_threshold=threshold)
        expected = 0

        assert actual == expected

    def test_calculateRecall_twoReportsHalfTruePositiveHalfFalseNegativeReturnsFifty(
        self
    ):
        columns = ["sample", "truth_probe_header", "vcf_probe_header", "classification"]
        report1 = pd.DataFrame(
            data=[
                create_report_row(AlignmentAssessment.PRIMARY_CORRECT, gt_conf=10),
                create_report_row(AlignmentAssessment.PRIMARY_CORRECT, gt_conf=100),
                create_report_row(AlignmentAssessment.PRIMARY_CORRECT, gt_conf=100),
                create_report_row(AlignmentAssessment.UNMAPPED, gt_conf=100),
                create_report_row("supplementary_incorrect", gt_conf=10),
                create_report_row("secondary_correct", gt_conf=100),
                create_report_row(AlignmentAssessment.PRIMARY_INCORRECT, gt_conf=100),
            ],
            columns=columns,
        )
        report2 = pd.DataFrame(
            data=[
                create_report_row(AlignmentAssessment.PRIMARY_CORRECT, gt_conf=100),
                create_report_row(AlignmentAssessment.UNMAPPED, gt_conf=100),
                create_report_row(AlignmentAssessment.PRIMARY_INCORRECT, gt_conf=100),
            ],
            columns=columns,
        )
        calculator = RecallCalculator([report1, report2])
        threshold = 60

        actual = calculator.calculate_recall(conf_threshold=threshold)
        expected = 0.5

        assert actual == expected
