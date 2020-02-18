import pandas as pd
from pandas.testing import assert_frame_equal
from evaluate.report import (
    PrecisionReport,
    RecallReport,
    Report
)
from evaluate.classification import AlignmentAssessment
import pytest
from io import StringIO
import math
from tests.common import create_tmp_file, create_recall_report_row, create_precision_report_row

class TestReport:
    def test_getMaximumGtConf_no_gt_conf_columnRaisesKeyError(self):
        calculator = Report([pd.DataFrame()])
        with pytest.raises(KeyError):
            calculator.get_maximum_gt_conf()

    def test_getMaximumGtConf_emptyReportReturnsNaN(self):
        calculator = Report([pd.DataFrame(data={"gt_conf": []})])
        actual = calculator.get_maximum_gt_conf()

        assert math.isnan(actual)

    def test_getMaximumGtConf_oneGTConfInReportReturnsGTConf(self):
        calculator = Report([pd.DataFrame(data={"gt_conf": [1.5]})])
        actual = calculator.get_maximum_gt_conf()
        expected = 1.5

        assert actual == expected

    def test_getMaximumGtConf_threeGTConfsInReportReturnsHighest(self):
        calculator = Report([pd.DataFrame(data={"gt_conf": [1.5, 10.5, 5.0]})])
        actual = calculator.get_maximum_gt_conf()
        expected = 10.5

        assert actual == expected


    def test_getMinimumGtConf_no_gt_conf_columnRaisesKeyError(self):
        calculator = Report([pd.DataFrame()])
        with pytest.raises(KeyError):
            calculator.get_minimum_gt_conf()

    def test_getMinimumGtConf_emptyReportReturnsNaN(self):
        calculator = Report([pd.DataFrame(data={"gt_conf": []})])
        actual = calculator.get_minimum_gt_conf()

        assert math.isnan(actual)

    def test_getMinimumGtConf_oneGTConfInReportReturnsGTConf(self):
        calculator = Report([pd.DataFrame(data={"gt_conf": [1.5]})])
        actual = calculator.get_minimum_gt_conf()
        expected = 1.5

        assert actual == expected

    def test_getMinimumGtConf_threeGTConfsInReportReturnsHighest(self):
        calculator = Report([pd.DataFrame(data={"gt_conf": [10.5, 5.0, 0.2]})])
        actual = calculator.get_minimum_gt_conf()
        expected = 0.2

        assert actual == expected


class TestRecallReport:
    def test_fromFiles_TwoFilesReturnsValidRecallReport(self):
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

        actual = RecallReport.from_files([path_1, path_2])
        expected = RecallReport(dataframes)

        path_1.unlink()
        path_2.unlink()

        assert actual == expected

    def test_init_gtconfIsExtractedCorrectly(self):
        columns = ["sample", "query_probe_header", "ref_probe_header", "classification"]
        dfs = pd.DataFrame(
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
        calculator = RecallReport([dfs])
        actual = calculator.report.gt_conf

        expected = pd.Series([100.0, 100.0, 10.0, 100.0])

        assert actual.equals(expected)


    def test_checkIfOnlyBestMappingIsKept_hasPrimaryMapping(self):
        dfs = pd.DataFrame(
            data=[
                create_recall_report_row("truth_probe_1", AlignmentAssessment.UNMAPPED, gt_conf=100, with_gt_conf=True),
                create_recall_report_row("truth_probe_1", AlignmentAssessment.PARTIALLY_MAPPED, gt_conf=100, with_gt_conf=True),
                create_recall_report_row("truth_probe_1", AlignmentAssessment.PRIMARY_INCORRECT, gt_conf=100, with_gt_conf=True),
                create_recall_report_row("truth_probe_1", AlignmentAssessment.SECONDARY_INCORRECT, gt_conf=100, with_gt_conf=True),
                create_recall_report_row("truth_probe_1", AlignmentAssessment.SUPPLEMENTARY_INCORRECT, gt_conf=100, with_gt_conf=True),
                create_recall_report_row("truth_probe_1", AlignmentAssessment.PRIMARY_CORRECT, gt_conf=100, with_gt_conf=True),
            ],
        )

        report = RecallReport([dfs])
        actual = report.report
        expected = pd.DataFrame(data=[create_recall_report_row("truth_probe_1", AlignmentAssessment.PRIMARY_CORRECT, gt_conf=100, with_gt_conf=True)])
        assert_frame_equal(actual, expected, check_dtype=False)

    def test_checkIfOnlyBestMappingIsKept_hasSecondaryMapping(self):
        dfs = pd.DataFrame(
            data=[
                create_recall_report_row("truth_probe_1", AlignmentAssessment.UNMAPPED, gt_conf=100, with_gt_conf=True),
                create_recall_report_row("truth_probe_1", AlignmentAssessment.PARTIALLY_MAPPED, gt_conf=100, with_gt_conf=True),
                create_recall_report_row("truth_probe_1", AlignmentAssessment.PRIMARY_INCORRECT, gt_conf=100, with_gt_conf=True),
                create_recall_report_row("truth_probe_1", AlignmentAssessment.SECONDARY_INCORRECT, gt_conf=100, with_gt_conf=True),
                create_recall_report_row("truth_probe_1", AlignmentAssessment.SUPPLEMENTARY_INCORRECT, gt_conf=100, with_gt_conf=True),
                create_recall_report_row("truth_probe_1", AlignmentAssessment.SECONDARY_CORRECT, gt_conf=100, with_gt_conf=True),
            ],
        )

        report = RecallReport([dfs])
        actual = report.report
        expected = pd.DataFrame(data=[create_recall_report_row("truth_probe_1", AlignmentAssessment.SECONDARY_CORRECT, gt_conf=100, with_gt_conf=True)])
        assert_frame_equal(actual, expected, check_dtype=False)

    def test_checkIfOnlyBestMappingIsKept_hasSupplementaryMapping(self):
        dfs = pd.DataFrame(
            data=[
                create_recall_report_row("truth_probe_1", AlignmentAssessment.UNMAPPED, gt_conf=100, with_gt_conf=True),
                create_recall_report_row("truth_probe_1", AlignmentAssessment.PARTIALLY_MAPPED, gt_conf=100, with_gt_conf=True),
                create_recall_report_row("truth_probe_1", AlignmentAssessment.PRIMARY_INCORRECT, gt_conf=100, with_gt_conf=True),
                create_recall_report_row("truth_probe_1", AlignmentAssessment.SECONDARY_INCORRECT, gt_conf=100, with_gt_conf=True),
                create_recall_report_row("truth_probe_1", AlignmentAssessment.SUPPLEMENTARY_INCORRECT, gt_conf=100, with_gt_conf=True),
                create_recall_report_row("truth_probe_1", AlignmentAssessment.SUPPLEMENTARY_CORRECT, gt_conf=100, with_gt_conf=True),
            ],
        )

        report = RecallReport([dfs])
        actual = report.report
        expected = pd.DataFrame(data=[create_recall_report_row("truth_probe_1", AlignmentAssessment.SUPPLEMENTARY_CORRECT, gt_conf=100, with_gt_conf=True)])
        assert_frame_equal(actual, expected, check_dtype=False)

    def test_checkIfOnlyBestMappingIsKept_ChoosesTheOneWithHighestGTConf(self):
        dfs = pd.DataFrame(
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


        report = RecallReport([dfs])
        actual = report.report
        expected = pd.DataFrame(data=[create_recall_report_row("truth_probe_1", AlignmentAssessment.SECONDARY_CORRECT, gt_conf=200, with_gt_conf=True)])
        assert_frame_equal(actual, expected, check_dtype=False)


    def test_checkIfOnlyBestMappingIsKept_hasNoCorrectMapping_ChoosesTheOneWithHighestGTConf(self):
        dfs = pd.DataFrame(
            data=[
                create_recall_report_row("truth_probe_1", AlignmentAssessment.UNMAPPED, gt_conf=100, with_gt_conf=True),
                create_recall_report_row("truth_probe_1", AlignmentAssessment.PARTIALLY_MAPPED, gt_conf=140, with_gt_conf=True),
                create_recall_report_row("truth_probe_1", AlignmentAssessment.PRIMARY_INCORRECT, gt_conf=150, with_gt_conf=True),
                create_recall_report_row("truth_probe_1", AlignmentAssessment.SECONDARY_INCORRECT, gt_conf=110, with_gt_conf=True),
                create_recall_report_row("truth_probe_1", AlignmentAssessment.SUPPLEMENTARY_INCORRECT, gt_conf=120, with_gt_conf=True),
            ],
        )

        report = RecallReport([dfs])
        actual = report.report
        expected = pd.DataFrame(data=[create_recall_report_row("truth_probe_1", AlignmentAssessment.PRIMARY_INCORRECT, gt_conf=150, with_gt_conf=True)])
        assert_frame_equal(actual, expected, check_dtype=False)


    def test_checkIfOnlyBestMappingIsKept_hasPrimaryCorrectMappingWithLowGTConf_ChoosesPrimaryCorrectMapping(self):
        dfs = pd.DataFrame(
            data=[
                create_recall_report_row("truth_probe_1", AlignmentAssessment.UNMAPPED, gt_conf=100, with_gt_conf=True),
                create_recall_report_row("truth_probe_1", AlignmentAssessment.PARTIALLY_MAPPED, gt_conf=140, with_gt_conf=True),
                create_recall_report_row("truth_probe_1", AlignmentAssessment.PRIMARY_INCORRECT, gt_conf=150, with_gt_conf=True),
                create_recall_report_row("truth_probe_1", AlignmentAssessment.SECONDARY_INCORRECT, gt_conf=110, with_gt_conf=True),
                create_recall_report_row("truth_probe_1", AlignmentAssessment.SUPPLEMENTARY_INCORRECT, gt_conf=120, with_gt_conf=True),
                create_recall_report_row("truth_probe_1", AlignmentAssessment.PRIMARY_CORRECT, gt_conf=1, with_gt_conf=True),
            ],
        )

        report = RecallReport([dfs])
        actual = report.report
        expected = pd.DataFrame(data=[create_recall_report_row("truth_probe_1", AlignmentAssessment.PRIMARY_CORRECT, gt_conf=1, with_gt_conf=True)])
        assert_frame_equal(actual, expected, check_dtype=False)


    def test_checkIfOnlyBestMappingIsKept_hasPrimaryMapping_and_several_dfs(self):
        df_1 = pd.DataFrame(
            data=[
                create_recall_report_row("truth_probe_1", AlignmentAssessment.UNMAPPED, gt_conf=100, with_gt_conf=True),
                create_recall_report_row("truth_probe_1", AlignmentAssessment.PARTIALLY_MAPPED, gt_conf=100, with_gt_conf=True),
            ],
        )
        df_2 = pd.DataFrame(
            data=[
                create_recall_report_row("truth_probe_1", AlignmentAssessment.PRIMARY_INCORRECT, gt_conf=100,
                                         with_gt_conf=True),
                create_recall_report_row("truth_probe_1", AlignmentAssessment.SECONDARY_INCORRECT, gt_conf=100,
                                         with_gt_conf=True),
            ],
        )
        df_3 = pd.DataFrame(
            data=[
                create_recall_report_row("truth_probe_1", AlignmentAssessment.SUPPLEMENTARY_INCORRECT, gt_conf=100,
                                         with_gt_conf=True),
                create_recall_report_row("truth_probe_1", AlignmentAssessment.PRIMARY_CORRECT, gt_conf=100,
                                         with_gt_conf=True),
            ],
        )

        report = RecallReport([df_1, df_2, df_3])
        actual = report.report
        expected = pd.DataFrame(data=[create_recall_report_row("truth_probe_1", AlignmentAssessment.PRIMARY_CORRECT, gt_conf=100, with_gt_conf=True)])
        assert_frame_equal(actual, expected, check_dtype=False)

    def test_checkIfOnlyBestMappingIsKept_hasNoCorrectMapping_ChoosesTheOneWithHighestGTConf_with_several_dfs(self):
        dfs = [pd.DataFrame(data=[create_recall_report_row("truth_probe_1", AlignmentAssessment.UNMAPPED, gt_conf=100, with_gt_conf=True)]),
               pd.DataFrame(data=[create_recall_report_row("truth_probe_1", AlignmentAssessment.PARTIALLY_MAPPED, gt_conf=140, with_gt_conf=True)]),
               pd.DataFrame(data=[create_recall_report_row("truth_probe_1", AlignmentAssessment.PRIMARY_INCORRECT, gt_conf=150, with_gt_conf=True)]),
               pd.DataFrame(data=[create_recall_report_row("truth_probe_1", AlignmentAssessment.SECONDARY_INCORRECT, gt_conf=110, with_gt_conf=True)]),
               pd.DataFrame(data=[create_recall_report_row("truth_probe_1", AlignmentAssessment.SUPPLEMENTARY_INCORRECT, gt_conf=120, with_gt_conf=True)])]
        report = RecallReport(dfs)
        actual = report.report
        expected = pd.DataFrame(data=[create_recall_report_row("truth_probe_1", AlignmentAssessment.PRIMARY_INCORRECT, gt_conf=150, with_gt_conf=True)])
        assert_frame_equal(actual, expected, check_dtype=False)


class TestPrecisionReporter:
    def test_init_gtconfIsExtractedCorrectly(self):
        columns = ["sample", "query_probe_header", "ref_probe_header", "classification"]
        dfs = pd.DataFrame(
            data=[
                create_precision_report_row(0.0, gt_conf=100),
                create_precision_report_row(0.0, gt_conf=100),
                create_precision_report_row(0.0, gt_conf=10),
                create_precision_report_row(0.0, gt_conf=100),
            ],
            columns=columns,
        )
        calculator = PrecisionReport([dfs])
        actual = calculator.report.gt_conf

        expected = pd.Series([100.0, 100.0, 10.0, 100.0])

        assert actual.equals(expected)

