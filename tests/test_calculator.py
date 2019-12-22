import pandas as pd
import tempfile
from evaluate.probe import ProbeHeader
from evaluate.calculator import (
    RecallCalculator,
    StatisticalClassification,
    PrecisionCalculator,
    Calculator,
    EmptyReportError,
)
from evaluate.classification import AlignmentAssessment
from pathlib import Path

import pytest
from io import StringIO
import math


def create_tmp_file(contents: str) -> Path:
    with tempfile.NamedTemporaryFile(mode="r+", delete=False) as tmp:
        tmp.write(contents)
        tmp.truncate()

    return Path(tmp.name)


def create_recall_report_row(
    truth_probe_header:str, classification: AlignmentAssessment, gt_conf: float = 0, sample: str = "sample1", with_gt_conf=False
) -> pd.Series:
    vcf_probe_header = ProbeHeader(gt_conf=gt_conf)
    data = {
        "sample": sample,
        "query_probe_header": str(truth_probe_header),
        "ref_probe_header": str(vcf_probe_header),
        "classification": classification.value,
    }
    if with_gt_conf:
        data["gt_conf"] = gt_conf

    return pd.Series(data=data)


def create_precision_report_row(
    classification: float, gt_conf: float = 0, sample: str = "sample1"
) -> pd.Series:
    ref_probe_header = ProbeHeader()
    pandora_probe_header = ProbeHeader(gt_conf=gt_conf)
    data = {
        "sample": sample,
        "query_probe_header": str(pandora_probe_header),
        "ref_probe_header": str(ref_probe_header),
        "classification": classification,
    }
    return pd.Series(data=data)


class TestCalculator:
    def test_getMaximumGtConf_no_gt_conf_columnRaisesKeyError(self):
        calculator = Calculator([pd.DataFrame()])
        with pytest.raises(KeyError):
            calculator.get_maximum_gt_conf()

    def test_getMaximumGtConf_emptyReportReturnsNaN(self):
        calculator = Calculator([pd.DataFrame(data={"gt_conf": []})])
        actual = calculator.get_maximum_gt_conf()

        assert math.isnan(actual)

    def test_getMaximumGtConf_oneGTConfInReportReturnsGTConf(self):
        calculator = Calculator([pd.DataFrame(data={"gt_conf": [1.5]})])
        actual = calculator.get_maximum_gt_conf()
        expected = 1.5

        assert actual == expected

    def test_getMaximumGtConf_threeGTConfsInReportReturnsHighest(self):
        calculator = Calculator([pd.DataFrame(data={"gt_conf": [1.5, 10.5, 5.0]})])
        actual = calculator.get_maximum_gt_conf()
        expected = 10.5

        assert actual == expected


    def test_getMinimumGtConf_no_gt_conf_columnRaisesKeyError(self):
        calculator = Calculator([pd.DataFrame()])
        with pytest.raises(KeyError):
            calculator.get_minimum_gt_conf()

    def test_getMinimumGtConf_emptyReportReturnsNaN(self):
        calculator = Calculator([pd.DataFrame(data={"gt_conf": []})])
        actual = calculator.get_minimum_gt_conf()

        assert math.isnan(actual)

    def test_getMinimumGtConf_oneGTConfInReportReturnsGTConf(self):
        calculator = Calculator([pd.DataFrame(data={"gt_conf": [1.5]})])
        actual = calculator.get_minimum_gt_conf()
        expected = 1.5

        assert actual == expected

    def test_getMinimumGtConf_threeGTConfsInReportReturnsHighest(self):
        calculator = Calculator([pd.DataFrame(data={"gt_conf": [10.5, 5.0, 0.2]})])
        actual = calculator.get_minimum_gt_conf()
        expected = 0.2

        assert actual == expected


class TestRecallCalculator:
    def test_fromFiles_TwoFilesReturnsValidRecallCalculator(self):
        contents_1 = """sample	query_probe_header	ref_probe_header	classification
CFT073	>CHROM=1;POS=1246;INTERVAL=[20,30);		unmapped
CFT073	>CHROM=1;POS=1248;INTERVAL=[30,40);	>CHROM=GC00005358_3;SAMPLE=CFT073;POS=1;INTERVAL=[0,17);SVTYPE=PH_SNPs;MEAN_FWD_COVG=3;MEAN_REV_COVG=6;GT_CONF=60.1133;	primary_correct
CFT073	>CHROM=1;POS=1252;INTERVAL=[40,50);		unmapped
"""
        contents_2 = """sample	query_probe_header	ref_probe_header	classification
CFT073	>CHROM=1;POS=1260;INTERVAL=[50,60);	>CHROM=GC00000578_3;SAMPLE=CFT073;POS=165;INTERVAL=[25,29);SVTYPE=PH_SNPs;MEAN_FWD_COVG=3;MEAN_REV_COVG=3;GT_CONF=3.22199;	primary_incorrect
CFT073	>CHROM=1;POS=1262;INTERVAL=[60,70);		unmapped
CFT073	>CHROM=1;POS=1281;INTERVAL=[70,80);		unmapped
"""
        path_1 = create_tmp_file(contents_1)
        path_2 = create_tmp_file(contents_2)

        contents_1_input = StringIO(contents_1)
        contents_2_input = StringIO(contents_2)
        dataframes = [
            pd.read_csv(contents_1_input, sep="\t", keep_default_na=False),
            pd.read_csv(contents_2_input, sep="\t", keep_default_na=False),
        ]

        actual = RecallCalculator.from_files([path_1, path_2])
        expected = RecallCalculator(dataframes)

        path_1.unlink()
        path_2.unlink()

        assert actual == expected

    def test_init_gtconfIsExtractedCorrectly(self):
        columns = ["sample", "query_probe_header", "ref_probe_header", "classification"]
        report = pd.DataFrame(
            data=[
                create_recall_report_row("truth_probe_1", AlignmentAssessment.UNMAPPED, gt_conf=100),
                create_recall_report_row("truth_probe_2", AlignmentAssessment.UNMAPPED, gt_conf=100),
                create_recall_report_row("truth_probe_3",
                    AlignmentAssessment.PRIMARY_CORRECT, gt_conf=10
                ),
                create_recall_report_row("truth_probe_4",
                    AlignmentAssessment.PRIMARY_INCORRECT, gt_conf=100
                ),
            ],
            columns=columns,
        )
        calculator = RecallCalculator([report])
        actual = calculator.report.gt_conf

        expected = pd.Series([100.0, 100.0, 10.0, 100.0])

        assert actual.equals(expected)


    def test__getBestMappingForTruthProbe_hasPrimaryMapping(self):
        report = pd.DataFrame(
            data=[
                create_recall_report_row("truth_probe_1", AlignmentAssessment.UNMAPPED, gt_conf=100, with_gt_conf=True),
                create_recall_report_row("truth_probe_1", AlignmentAssessment.PARTIALLY_MAPPED, gt_conf=100, with_gt_conf=True),
                create_recall_report_row("truth_probe_1", AlignmentAssessment.PRIMARY_INCORRECT, gt_conf=100, with_gt_conf=True),
                create_recall_report_row("truth_probe_1", AlignmentAssessment.SECONDARY_INCORRECT, gt_conf=100, with_gt_conf=True),
                create_recall_report_row("truth_probe_1", AlignmentAssessment.SUPPLEMENTARY_INCORRECT, gt_conf=100, with_gt_conf=True),
                create_recall_report_row("truth_probe_1", AlignmentAssessment.PRIMARY_CORRECT, gt_conf=100, with_gt_conf=True),
            ],
        )
        report = RecallCalculator._get_truth_probe_to_all_mappings_dfs(report)
        actual = RecallCalculator._get_best_mapping_for_truth_probe(report, "truth_probe_1")
        expected = create_recall_report_row("truth_probe_1", AlignmentAssessment.PRIMARY_CORRECT, gt_conf=100, with_gt_conf=True)
        expected.classification = StatisticalClassification.TRUE_POSITIVE

        assert actual.equals(expected)

    def test__getBestMappingForTruthProbe_hasSecondaryMapping(self):
        report = pd.DataFrame(
            data=[
                create_recall_report_row("truth_probe_1", AlignmentAssessment.UNMAPPED, gt_conf=100, with_gt_conf=True),
                create_recall_report_row("truth_probe_1", AlignmentAssessment.PARTIALLY_MAPPED, gt_conf=100, with_gt_conf=True),
                create_recall_report_row("truth_probe_1", AlignmentAssessment.PRIMARY_INCORRECT, gt_conf=100, with_gt_conf=True),
                create_recall_report_row("truth_probe_1", AlignmentAssessment.SECONDARY_INCORRECT, gt_conf=100, with_gt_conf=True),
                create_recall_report_row("truth_probe_1", AlignmentAssessment.SUPPLEMENTARY_INCORRECT, gt_conf=100, with_gt_conf=True),
                create_recall_report_row("truth_probe_1", AlignmentAssessment.SECONDARY_CORRECT, gt_conf=100, with_gt_conf=True),
            ],
        )
        report = RecallCalculator._get_truth_probe_to_all_mappings_dfs(report)
        actual = RecallCalculator._get_best_mapping_for_truth_probe(report, "truth_probe_1")
        expected = create_recall_report_row("truth_probe_1", AlignmentAssessment.SECONDARY_CORRECT, gt_conf=100, with_gt_conf=True)
        expected.classification = StatisticalClassification.TRUE_POSITIVE

        assert actual.equals(expected)

    def test__getBestMappingForTruthProbe_hasSupplementaryMapping(self):
        report = pd.DataFrame(
            data=[
                create_recall_report_row("truth_probe_1", AlignmentAssessment.UNMAPPED, gt_conf=100, with_gt_conf=True),
                create_recall_report_row("truth_probe_1", AlignmentAssessment.PARTIALLY_MAPPED, gt_conf=100, with_gt_conf=True),
                create_recall_report_row("truth_probe_1", AlignmentAssessment.PRIMARY_INCORRECT, gt_conf=100, with_gt_conf=True),
                create_recall_report_row("truth_probe_1", AlignmentAssessment.SECONDARY_INCORRECT, gt_conf=100, with_gt_conf=True),
                create_recall_report_row("truth_probe_1", AlignmentAssessment.SUPPLEMENTARY_INCORRECT, gt_conf=100, with_gt_conf=True),
                create_recall_report_row("truth_probe_1", AlignmentAssessment.SUPPLEMENTARY_CORRECT, gt_conf=100, with_gt_conf=True),
            ],
        )
        report = RecallCalculator._get_truth_probe_to_all_mappings_dfs(report)
        actual = RecallCalculator._get_best_mapping_for_truth_probe(report, "truth_probe_1")
        expected = create_recall_report_row("truth_probe_1", AlignmentAssessment.SUPPLEMENTARY_CORRECT, gt_conf=100, with_gt_conf=True)
        expected.classification = StatisticalClassification.TRUE_POSITIVE

        assert actual.equals(expected)

    def test__getBestMappingForTruthProbe_hasPrimarySecondaryAndSupplementaryMapping_ChoosesTheOneWithHighestGTConf(self):
        report = pd.DataFrame(
            data=[
                create_recall_report_row("truth_probe_1", AlignmentAssessment.UNMAPPED, gt_conf=100, with_gt_conf=True),
                create_recall_report_row("truth_probe_1", AlignmentAssessment.PARTIALLY_MAPPED, gt_conf=100, with_gt_conf=True),
                create_recall_report_row("truth_probe_1", AlignmentAssessment.PRIMARY_INCORRECT, gt_conf=100, with_gt_conf=True),
                create_recall_report_row("truth_probe_1", AlignmentAssessment.SECONDARY_INCORRECT, gt_conf=100, with_gt_conf=True),
                create_recall_report_row("truth_probe_1", AlignmentAssessment.SUPPLEMENTARY_INCORRECT, gt_conf=100, with_gt_conf=True),
                create_recall_report_row("truth_probe_1", AlignmentAssessment.PRIMARY_CORRECT, gt_conf=100, with_gt_conf=True),
                create_recall_report_row("truth_probe_1", AlignmentAssessment.SECONDARY_CORRECT, gt_conf=200, with_gt_conf=True),
                create_recall_report_row("truth_probe_1", AlignmentAssessment.SUPPLEMENTARY_CORRECT, gt_conf=150, with_gt_conf=True),
            ],
        )
        report = RecallCalculator._get_truth_probe_to_all_mappings_dfs(report)
        actual = RecallCalculator._get_best_mapping_for_truth_probe(report, "truth_probe_1")
        expected = create_recall_report_row("truth_probe_1", AlignmentAssessment.SECONDARY_CORRECT, gt_conf=200, with_gt_conf=True)
        expected.classification = StatisticalClassification.TRUE_POSITIVE

        assert actual.equals(expected)


    def test__getBestMappingForTruthProbe_hasNoCorrectMapping_ChoosesTheOneWithHighestGTConf(self):
        report = pd.DataFrame(
            data=[
                create_recall_report_row("truth_probe_1", AlignmentAssessment.UNMAPPED, gt_conf=100, with_gt_conf=True),
                create_recall_report_row("truth_probe_1", AlignmentAssessment.PARTIALLY_MAPPED, gt_conf=140, with_gt_conf=True),
                create_recall_report_row("truth_probe_1", AlignmentAssessment.PRIMARY_INCORRECT, gt_conf=150, with_gt_conf=True),
                create_recall_report_row("truth_probe_1", AlignmentAssessment.SECONDARY_INCORRECT, gt_conf=110, with_gt_conf=True),
                create_recall_report_row("truth_probe_1", AlignmentAssessment.SUPPLEMENTARY_INCORRECT, gt_conf=120, with_gt_conf=True),
            ],
        )
        report = RecallCalculator._get_truth_probe_to_all_mappings_dfs(report)
        actual = RecallCalculator._get_best_mapping_for_truth_probe(report, "truth_probe_1")
        expected = create_recall_report_row("truth_probe_1", AlignmentAssessment.PRIMARY_INCORRECT, gt_conf=150, with_gt_conf=True)
        expected.classification = StatisticalClassification.FALSE_NEGATIVE

        assert actual.equals(expected)


    def test__getBestMappingForTruthProbe_hasPrimaryCorrectMappingWithLowGTConf_ChoosesPrimaryCorrectMapping(self):
        report = pd.DataFrame(
            data=[
                create_recall_report_row("truth_probe_1", AlignmentAssessment.UNMAPPED, gt_conf=100, with_gt_conf=True),
                create_recall_report_row("truth_probe_1", AlignmentAssessment.PARTIALLY_MAPPED, gt_conf=140, with_gt_conf=True),
                create_recall_report_row("truth_probe_1", AlignmentAssessment.PRIMARY_INCORRECT, gt_conf=150, with_gt_conf=True),
                create_recall_report_row("truth_probe_1", AlignmentAssessment.SECONDARY_INCORRECT, gt_conf=110, with_gt_conf=True),
                create_recall_report_row("truth_probe_1", AlignmentAssessment.SUPPLEMENTARY_INCORRECT, gt_conf=120, with_gt_conf=True),
                create_recall_report_row("truth_probe_1", AlignmentAssessment.PRIMARY_CORRECT, gt_conf=1, with_gt_conf=True),
            ],
        )
        report = RecallCalculator._get_truth_probe_to_all_mappings_dfs(report)
        actual = RecallCalculator._get_best_mapping_for_truth_probe(report, "truth_probe_1")
        expected = create_recall_report_row("truth_probe_1", AlignmentAssessment.PRIMARY_CORRECT, gt_conf=1, with_gt_conf=True)
        expected.classification = StatisticalClassification.TRUE_POSITIVE

        assert actual.equals(expected)


    def test_calculateRecall_noReportsRaisesEmptyReportError(self):
        columns = ["sample", "query_probe_header", "ref_probe_header", "classification"]
        report = pd.DataFrame(columns=columns)
        calculator = RecallCalculator([report])
        threshold = 0

        with pytest.raises(EmptyReportError):
            calculator.calculate_recall(conf_threshold=threshold)

    def test_calculateRecall_oneReportNoTruePositivesReturnsZero(self):
        columns = ["sample", "query_probe_header", "ref_probe_header", "classification"]
        report = pd.DataFrame(
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
        calculator = RecallCalculator([report])
        threshold = 60

        actual = calculator.calculate_recall(conf_threshold=threshold)

        assert actual.recall == 0.0
        assert actual.true_positives == 0.0
        assert actual.total == 4.0

    def test_calculateRecall_oneReportNoTruePositivesTwoTruthProbesReturnsZero(self):
        columns = ["sample", "query_probe_header", "ref_probe_header", "classification"]
        report = pd.DataFrame (
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
        calculator = RecallCalculator([report])
        threshold = 60

        actual = calculator.calculate_recall(conf_threshold=threshold)

        assert actual.recall == 0.0
        assert actual.true_positives == 0.0
        assert actual.total == 2.0

    def test_calculateRecall_oneReportNoFalseNegativesReturnsOne(self):
        columns = ["sample", "query_probe_header", "ref_probe_header", "classification"]
        report = pd.DataFrame(
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
        calculator = RecallCalculator([report])
        threshold = 60

        actual = calculator.calculate_recall(conf_threshold=threshold)

        assert actual.recall == 1.0
        assert actual.true_positives == 4.0
        assert actual.total == 4.0


    def test_calculateRecall_oneReportNoFalseNegativesTwoProbesReturnsOne(self):
        columns = ["sample", "query_probe_header", "ref_probe_header", "classification"]
        report = pd.DataFrame(
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
        calculator = RecallCalculator([report])
        threshold = 60

        actual = calculator.calculate_recall(conf_threshold=threshold)

        assert actual.recall == 1.0
        assert actual.true_positives == 2.0
        assert actual.total == 2.0


    def test_calculateRecall_oneReportHalfTruePositiveHalfFalseNegativeReturnsFifty(
        self
    ):
        columns = ["sample", "query_probe_header", "ref_probe_header", "classification"]
        report = pd.DataFrame(
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
        calculator = RecallCalculator([report])
        threshold = 60

        actual = calculator.calculate_recall(conf_threshold=threshold)

        assert actual.recall == 0.5
        assert actual.true_positives == 3.0
        assert actual.total == 6.0


    def test_calculateRecall_oneReportAllTruePositivesAllBelowThresholdReturnsZero(
        self
    ):
        columns = ["sample", "query_probe_header", "ref_probe_header", "classification"]
        report = pd.DataFrame(
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
        calculator = RecallCalculator([report])
        threshold = 100

        actual = calculator.calculate_recall(conf_threshold=threshold)

        assert actual.recall == 0.0
        assert actual.true_positives == 0.0
        assert actual.total == 2.0

    def test_calculateRecall_oneReportAllTruePositivesAllBelowThresholdOneTruthProbeReturnsZero(
        self
    ):
        columns = ["sample", "query_probe_header", "ref_probe_header", "classification"]
        report = pd.DataFrame(
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
        calculator = RecallCalculator([report])
        threshold = 100

        actual = calculator.calculate_recall(conf_threshold=threshold)

        assert actual.recall == 0.0
        assert actual.true_positives == 0.0
        assert actual.total == 1.0


    def test_calculateRecall_AllAlignmentsAreIncorrectAndOneIsCorrect(
        self
    ):
        columns = ["sample", "query_probe_header", "ref_probe_header", "classification"]
        report = pd.DataFrame(
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
        calculator = RecallCalculator([report])
        threshold = 60

        actual = calculator.calculate_recall(conf_threshold=threshold)

        assert actual.recall == 0.25
        assert actual.true_positives == 1.0
        assert actual.total == 4.0


    def test_calculateRecall_AllAlignmentsAreIncorrectAndOneIsCorrectTwoTruthProbes(
        self
    ):
        columns = ["sample", "query_probe_header", "ref_probe_header", "classification"]
        report = pd.DataFrame(
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
        calculator = RecallCalculator([report])
        threshold = 60

        actual = calculator.calculate_recall(conf_threshold=threshold)

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
        calculator = RecallCalculator([report1, report2])
        threshold = 60

        actual = calculator.calculate_recall(conf_threshold=threshold)

        assert actual.recall == 0.5
        assert actual.true_positives == 5.0
        assert actual.total == 10.0


class TestPrecisionCalculator:
    def test_init_gtconfIsExtractedCorrectly(self):
        columns = ["sample", "query_probe_header", "ref_probe_header", "classification"]
        report = pd.DataFrame(
            data=[
                create_precision_report_row(0.0, gt_conf=100),
                create_precision_report_row(0.0, gt_conf=101),
                create_precision_report_row(0.0, gt_conf=10),
                create_precision_report_row(0.0, gt_conf=102),
            ],
            columns=columns,
        )
        calculator = PrecisionCalculator([report])
        actual = calculator.report.gt_conf

        expected = pd.Series([100.0, 101.0, 10.0, 102.0])

        assert actual.equals(expected)

    def test_calculatePrecision_NoReportsRaisesEmptyReportError(self):
        columns = ["sample", "query_probe_header", "ref_probe_header", "classification"]
        report = pd.DataFrame(columns=columns)
        calculator = PrecisionCalculator([report])

        with pytest.raises(EmptyReportError):
            calculator.calculate_precision()

    def test_calculatePrecision_OneReportWithOneRowCompletelyCorrectReturnsOne(self):
        columns = ["sample", "query_probe_header", "ref_probe_header", "classification"]
        report = pd.DataFrame(
            data=[create_precision_report_row(1.0, gt_conf=100)], columns=columns
        )
        calculator = PrecisionCalculator([report])

        actual = calculator.calculate_precision()

        assert actual.precision == 1.0
        assert actual.true_positives == 1.0
        assert actual.total == 1.0

    def test_calculatePrecision_OneReportWithOneRowCompletelyIncorrectReturnsZero(self):
        columns = ["sample", "query_probe_header", "ref_probe_header", "classification"]
        report = pd.DataFrame(
            data=[create_precision_report_row(0.0, gt_conf=100)], columns=columns
        )
        calculator = PrecisionCalculator([report])

        actual = calculator.calculate_precision()

        assert actual.precision == 0.0
        assert actual.true_positives == 0.0
        assert actual.total == 1.0

    def test_calculatePrecision_OneReportWithOneRowCompletelyCorrectBelowConfThreasholdRaisesEmptyReportError(
        self
    ):
        columns = ["sample", "query_probe_header", "ref_probe_header", "classification"]
        report = pd.DataFrame(
            data=[create_precision_report_row(1.0, gt_conf=10)], columns=columns
        )
        calculator = PrecisionCalculator([report])
        confidence_threshold = 60

        with pytest.raises(EmptyReportError):
            calculator.calculate_precision(confidence_threshold)

    def test_calculatePrecision_OneReportWithOneRowCompletelyCorrectEqualConfThreasholdReturnsOne(
        self
    ):
        columns = ["sample", "query_probe_header", "ref_probe_header", "classification"]
        report = pd.DataFrame(
            data=[create_precision_report_row(1.0, gt_conf=60)], columns=columns
        )
        calculator = PrecisionCalculator([report])
        confidence_threshold = 60

        actual = calculator.calculate_precision(confidence_threshold)

        assert actual.precision == 1.0
        assert actual.true_positives == 1.0
        assert actual.total == 1.0

    def test_calculatePrecision_OneReportWithTwoRowsPartiallyCorrect(self):
        columns = ["sample", "query_probe_header", "ref_probe_header", "classification"]
        report = pd.DataFrame(
            data=[
                create_precision_report_row(0.5, gt_conf=100),
                create_precision_report_row(0.7, gt_conf=100),
            ],
            columns=columns,
        )
        calculator = PrecisionCalculator([report])

        actual = calculator.calculate_precision()

        assert actual.precision == 1.2/2
        assert actual.true_positives == 1.2
        assert actual.total == 2.0


    def test_calculatePrecision_OneReportWithThreeRowsTwoPartiallyCorrectOneBelowThreshold(
        self
    ):
        columns = ["sample", "query_probe_header", "ref_probe_header", "classification"]
        report = pd.DataFrame(
            data=[
                create_precision_report_row(0.4, gt_conf=100),
                create_precision_report_row(0.8, gt_conf=20),
                create_precision_report_row(0.3, gt_conf=100),
            ],
            columns=columns,
        )
        calculator = PrecisionCalculator([report])
        confidence_threshold = 80

        actual = calculator.calculate_precision(confidence_threshold)

        assert actual.precision == 0.7/2.0
        assert actual.true_positives == 0.7
        assert actual.total == 2.0

