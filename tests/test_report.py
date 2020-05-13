import pandas as pd
from pandas.testing import assert_frame_equal
from evaluate.report import (
    PrecisionReport,
    RecallReport,
    Report,
    DelimNotFoundError,
    ReturnTypeDoesNotMatchError
)
from evaluate.classification import AlignmentAssessment
import pytest
from io import StringIO
import math
from tests.common import create_tmp_file, create_recall_report_row, create_precision_report_row
from unittest.mock import patch

class TestReport:
    def test___get_report_satisfying_confidence_threshold(self):
        report = Report([
            pd.read_csv(StringIO(
"""id,GT_CONF
0,2
1,1
2,3
""")),
            pd.read_csv(StringIO(
"""id,GT_CONF
4,3
5,1
6,2
"""))
        ])

        actual_report = report.get_report_satisfying_confidence_threshold(2)
        expected_report = Report([
            pd.read_csv(StringIO(
"""id,GT_CONF
0,2
2,3
4,3
6,2
"""))])
        assert actual_report==expected_report

    def test___get_value_from_header_fast___field_is_in_header(self):
        actual_value = Report.get_value_from_header_fast("FIELD_1=10;", "FIELD_1", int, -1, delim=";")
        expected_value = 10
        assert actual_value==expected_value

    def test___get_value_from_header_fast___field_is_in_header_between_two_other_fields(self):
        actual_value = Report.get_value_from_header_fast("DUMMY_1=asd;FIELD_1=10;DUMMY_2=99;", "FIELD_1", int, -1, delim=";")
        expected_value = 10
        assert actual_value==expected_value

    def test___get_value_from_header_fast___field_is_first_before_two_other_fields(self):
        actual_value = Report.get_value_from_header_fast("FIELD_1=10;DUMMY_1=asd;DUMMY_2=99;", "FIELD_1", int, -1, delim=";")
        expected_value = 10
        assert actual_value==expected_value

    def test___get_value_from_header_fast___field_is_last_after_two_other_fields(self):
        actual_value = Report.get_value_from_header_fast("DUMMY_1=asd;DUMMY_2=99;FIELD_1=10;", "FIELD_1", int, -1, delim=";")
        expected_value = 10
        assert actual_value==expected_value

    def test___get_value_from_header_fast___field_is_not_in_header(self):
        actual_value = Report.get_value_from_header_fast("DUMMY_1=asd;FIELD_1=10;DUMMY_2=99;", "FIELD_2", int, -1, delim=";")
        expected_value = -1
        assert actual_value==expected_value

    def test___get_value_from_header_fast___field_is_in_header___return_type_does_not_match(self):
        with pytest.raises(ReturnTypeDoesNotMatchError):
            Report.get_value_from_header_fast("DUMMY_1=asd;FIELD_1=asd;DUMMY_2=99;", "FIELD_1", int, -1, delim=";")

    def test___get_value_from_header_fast___field_is_in_header___delim_is_not(self):
        with pytest.raises(DelimNotFoundError):
            Report.get_value_from_header_fast("DUMMY_1=asd;FIELD_1=asd;DUMMY_2=99;", "FIELD_1", int, -1, delim="~")

    def test____create_field_from_header(self):
        report = Report([
            pd.read_csv(StringIO(
"""id,header
1,SEQ=ACGT;LEN=4;
2,SEQ=TG;LEN=2;
3,dummy
"""))])
        report._create_field_from_header("SEQ", "header", str, "A")
        report._create_field_from_header("LEN", "header", int, 1)

        expected_report = Report([
            pd.read_csv(StringIO(
"""id,header,SEQ,LEN
1,SEQ=ACGT;LEN=4;,ACGT,4
2,SEQ=TG;LEN=2;,TG,2
3,dummy,A,1
"""))])
        assert report==expected_report


    def test____create_good_eval_column(self):
        report = Report([
            pd.read_csv(StringIO(
"""classification
primary_correct
whatever
secondary_correct
dummy
supplementary_correct
woot
"""))])
        report._create_good_eval_column()

        expected_report = Report([
            pd.read_csv(StringIO(
"""classification,good_eval
primary_correct,True
whatever,False
secondary_correct,True
dummy,False
supplementary_correct,True
woot,False
"""))])
        assert report==expected_report

    def test_getMaximumGtConf_no_gt_conf_columnRaisesKeyError(self):
        report = Report([pd.DataFrame()])
        with pytest.raises(KeyError):
            report.get_maximum_gt_conf()

    def test_getMaximumGtConf_emptyReportReturnsNaN(self):
        report = Report([pd.DataFrame(data={"GT_CONF": []})])
        actual = report.get_maximum_gt_conf()

        assert math.isnan(actual)

    def test_getMaximumGtConf_oneGTConfInReportReturnsGTConf(self):
        report = Report([pd.DataFrame(data={"GT_CONF": [1.5]})])
        actual = report.get_maximum_gt_conf()
        expected = 1.5

        assert actual == expected

    def test_getMaximumGtConf_threeGTConfsInReportReturnsHighest(self):
        report = Report([pd.DataFrame(data={"GT_CONF": [1.5, 10.5, 5.0]})])
        actual = report.get_maximum_gt_conf()
        expected = 10.5

        assert actual == expected


    def test_getMinimumGtConf_no_gt_conf_columnRaisesKeyError(self):
        report = Report([pd.DataFrame()])
        with pytest.raises(KeyError):
            report.get_minimum_gt_conf()

    def test_getMinimumGtConf_emptyReportReturnsNaN(self):
        report = Report([pd.DataFrame(data={"GT_CONF": []})])
        actual = report.get_minimum_gt_conf()

        assert math.isnan(actual)

    def test_getMinimumGtConf_oneGTConfInReportReturnsGTConf(self):
        report = Report([pd.DataFrame(data={"GT_CONF": [1.5]})])
        actual = report.get_minimum_gt_conf()
        expected = 1.5

        assert actual == expected

    def test_getMinimumGtConf_threeGTConfsInReportReturnsHighest(self):
        report = Report([pd.DataFrame(data={"GT_CONF": [10.5, 5.0, 0.2]})])
        actual = report.get_minimum_gt_conf()
        expected = 0.2

        assert actual == expected




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
        report = PrecisionReport([dfs])
        actual = report.report.GT_CONF

        expected = pd.Series([100.0, 100.0, 10.0, 100.0])

        assert actual.equals(expected)

    def test_fromFiles_TwoFilesReturnsValidRecallReport(self):
        contents_1 = """sample	query_probe_header	ref_probe_header	classification
CFT073	>CHROM=1;POS=1246;INTERVAL=[20,30);PANGENOME_VARIATION_ID=1;NUMBER_OF_ALLELES=1;ALLELE_ID=1;NUMBER_OF_DIFFERENT_ALLELE_SEQUENCES=1;ALLELE_SEQUENCE_ID=1;	>GT_CONF=1;	unmapped
CFT073	>CHROM=1;POS=1248;INTERVAL=[30,40);PANGENOME_VARIATION_ID=2;NUMBER_OF_ALLELES=2;ALLELE_ID=2;NUMBER_OF_DIFFERENT_ALLELE_SEQUENCES=2;ALLELE_SEQUENCE_ID=2;	>CHROM=GC00005358_3;SAMPLE=CFT073;POS=1;INTERVAL=[0,17);SVTYPE=PH_SNPs;MEAN_FWD_COVG=3;MEAN_REV_COVG=6;GT_CONF=60.1133;	primary_correct
CFT073	>CHROM=1;POS=1252;INTERVAL=[40,50);PANGENOME_VARIATION_ID=3;NUMBER_OF_ALLELES=3;ALLELE_ID=3;NUMBER_OF_DIFFERENT_ALLELE_SEQUENCES=3;ALLELE_SEQUENCE_ID=3;	>GT_CONF=3;	unmapped
"""
        contents_2 = """sample	query_probe_header	ref_probe_header	classification
CFT073	>CHROM=1;POS=1260;INTERVAL=[50,60);PANGENOME_VARIATION_ID=4;NUMBER_OF_ALLELES=4;ALLELE_ID=4;NUMBER_OF_DIFFERENT_ALLELE_SEQUENCES=4;ALLELE_SEQUENCE_ID=4;	>CHROM=GC00000578_3;SAMPLE=CFT073;POS=165;INTERVAL=[25,29);SVTYPE=PH_SNPs;MEAN_FWD_COVG=3;MEAN_REV_COVG=3;GT_CONF=3.22199;	primary_incorrect
CFT073	>CHROM=1;POS=1262;INTERVAL=[60,70);PANGENOME_VARIATION_ID=5;NUMBER_OF_ALLELES=5;ALLELE_ID=5;NUMBER_OF_DIFFERENT_ALLELE_SEQUENCES=5;ALLELE_SEQUENCE_ID=5;	>GT_CONF=5;	unmapped
CFT073	>CHROM=1;POS=1281;INTERVAL=[70,80);PANGENOME_VARIATION_ID=6;NUMBER_OF_ALLELES=6;ALLELE_ID=6;NUMBER_OF_DIFFERENT_ALLELE_SEQUENCES=6;ALLELE_SEQUENCE_ID=6;	>GT_CONF=6;	unmapped
"""
        path_1 = create_tmp_file(contents_1)
        path_2 = create_tmp_file(contents_2)

        contents_1_input = StringIO(contents_1)
        contents_2_input = StringIO(contents_2)
        dataframes = [
            pd.read_csv(contents_1_input, sep="\t", keep_default_na=False),
            pd.read_csv(contents_2_input, sep="\t", keep_default_na=False),
        ]

        actual = PrecisionReport.from_files([path_1, path_2])
        expected = PrecisionReport(dataframes)

        path_1.unlink()
        path_2.unlink()

        assert actual == expected



class TestRecallReport:
    def test_fromFiles_TwoFilesReturnsValidRecallReport(self):
        contents_1 = """sample	query_probe_header	ref_probe_header	classification
CFT073	>CHROM=1;POS=1246;INTERVAL=[20,30);PANGENOME_VARIATION_ID=1;NUMBER_OF_ALLELES=1;ALLELE_ID=1;NUMBER_OF_DIFFERENT_ALLELE_SEQUENCES=1;ALLELE_SEQUENCE_ID=1;	>GT_CONF=1;	unmapped
CFT073	>CHROM=1;POS=1248;INTERVAL=[30,40);PANGENOME_VARIATION_ID=2;NUMBER_OF_ALLELES=2;ALLELE_ID=2;NUMBER_OF_DIFFERENT_ALLELE_SEQUENCES=2;ALLELE_SEQUENCE_ID=2;	>CHROM=GC00005358_3;SAMPLE=CFT073;POS=1;INTERVAL=[0,17);SVTYPE=PH_SNPs;MEAN_FWD_COVG=3;MEAN_REV_COVG=6;GT_CONF=60.1133;	primary_correct
CFT073	>CHROM=1;POS=1252;INTERVAL=[40,50);PANGENOME_VARIATION_ID=3;NUMBER_OF_ALLELES=3;ALLELE_ID=3;NUMBER_OF_DIFFERENT_ALLELE_SEQUENCES=3;ALLELE_SEQUENCE_ID=3;	>GT_CONF=3;	unmapped
"""
        contents_2 = """sample	query_probe_header	ref_probe_header	classification
CFT073	>CHROM=1;POS=1260;INTERVAL=[50,60);PANGENOME_VARIATION_ID=4;NUMBER_OF_ALLELES=4;ALLELE_ID=4;NUMBER_OF_DIFFERENT_ALLELE_SEQUENCES=4;ALLELE_SEQUENCE_ID=4;	>CHROM=GC00000578_3;SAMPLE=CFT073;POS=165;INTERVAL=[25,29);SVTYPE=PH_SNPs;MEAN_FWD_COVG=3;MEAN_REV_COVG=3;GT_CONF=3.22199;	primary_incorrect
CFT073	>CHROM=1;POS=1262;INTERVAL=[60,70);PANGENOME_VARIATION_ID=5;NUMBER_OF_ALLELES=5;ALLELE_ID=5;NUMBER_OF_DIFFERENT_ALLELE_SEQUENCES=5;ALLELE_SEQUENCE_ID=5;	>GT_CONF=5;	unmapped
CFT073	>CHROM=1;POS=1281;INTERVAL=[70,80);PANGENOME_VARIATION_ID=6;NUMBER_OF_ALLELES=6;ALLELE_ID=6;NUMBER_OF_DIFFERENT_ALLELE_SEQUENCES=6;ALLELE_SEQUENCE_ID=6;	>GT_CONF=6;	unmapped
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

    def test_init(self):
        contents_1 = """sample	query_probe_header	ref_probe_header	classification
CFT073	>CHROM=1;POS=1246;INTERVAL=[20,30);PANGENOME_VARIATION_ID=1;NUMBER_OF_ALLELES=1;ALLELE_ID=1;NUMBER_OF_DIFFERENT_ALLELE_SEQUENCES=1;ALLELE_SEQUENCE_ID=1;	>GT_CONF=1;	unmapped
CFT073	>CHROM=1;POS=1248;INTERVAL=[30,40);PANGENOME_VARIATION_ID=2;NUMBER_OF_ALLELES=2;ALLELE_ID=2;NUMBER_OF_DIFFERENT_ALLELE_SEQUENCES=2;ALLELE_SEQUENCE_ID=2;	>CHROM=GC00005358_3;SAMPLE=CFT073;POS=1;INTERVAL=[0,17);SVTYPE=PH_SNPs;MEAN_FWD_COVG=3;MEAN_REV_COVG=6;GT_CONF=60.1133;	primary_correct
CFT073	>CHROM=1;POS=1252;INTERVAL=[40,50);PANGENOME_VARIATION_ID=3;NUMBER_OF_ALLELES=3;ALLELE_ID=3;NUMBER_OF_DIFFERENT_ALLELE_SEQUENCES=3;ALLELE_SEQUENCE_ID=3;	>GT_CONF=3;	unmapped
"""
        contents_1_input = StringIO(contents_1)
        dataframes = [pd.read_csv(contents_1_input, sep="\t", keep_default_na=False)]
        report = RecallReport(dataframes)
        actual_df = report.report
        expected_df = pd.read_csv(StringIO(
"""sample	query_probe_header	ref_probe_header	classification	GT_CONF	PANGENOME_VARIATION_ID	NUMBER_OF_ALLELES	ALLELE_ID	NUMBER_OF_DIFFERENT_ALLELE_SEQUENCES	ALLELE_SEQUENCE_ID	good_eval
CFT073	>CHROM=1;POS=1246;INTERVAL=[20,30);PANGENOME_VARIATION_ID=1;NUMBER_OF_ALLELES=1;ALLELE_ID=1;NUMBER_OF_DIFFERENT_ALLELE_SEQUENCES=1;ALLELE_SEQUENCE_ID=1;	>GT_CONF=1;	unmapped	1.0	1	1	1	1	1	False
CFT073	>CHROM=1;POS=1248;INTERVAL=[30,40);PANGENOME_VARIATION_ID=2;NUMBER_OF_ALLELES=2;ALLELE_ID=2;NUMBER_OF_DIFFERENT_ALLELE_SEQUENCES=2;ALLELE_SEQUENCE_ID=2;	>CHROM=GC00005358_3;SAMPLE=CFT073;POS=1;INTERVAL=[0,17);SVTYPE=PH_SNPs;MEAN_FWD_COVG=3;MEAN_REV_COVG=6;GT_CONF=60.1133;	primary_correct	60.1133	2	2	2	2	2	True
CFT073	>CHROM=1;POS=1252;INTERVAL=[40,50);PANGENOME_VARIATION_ID=3;NUMBER_OF_ALLELES=3;ALLELE_ID=3;NUMBER_OF_DIFFERENT_ALLELE_SEQUENCES=3;ALLELE_SEQUENCE_ID=3;	>GT_CONF=3;	unmapped	3.0	3	3	3	3	3	False
"""), sep="\t")

        assert actual_df.equals(expected_df)



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



    def test_simple_concatenation_with_several_dfs(self):
        df_1 = pd.DataFrame(
            data=[
                create_recall_report_row("truth_probe_1", AlignmentAssessment.UNMAPPED, gt_conf=100, with_gt_conf=True),
                create_recall_report_row("truth_probe_2", AlignmentAssessment.PARTIALLY_MAPPED, gt_conf=100, with_gt_conf=True),
            ],
        )
        df_2 = pd.DataFrame(
            data=[
                create_recall_report_row("truth_probe_3", AlignmentAssessment.PRIMARY_INCORRECT, gt_conf=100,
                                         with_gt_conf=True),
                create_recall_report_row("truth_probe_4", AlignmentAssessment.SECONDARY_INCORRECT, gt_conf=100,
                                         with_gt_conf=True),
            ],
        )
        df_3 = pd.DataFrame(
            data=[
                create_recall_report_row("truth_probe_5", AlignmentAssessment.SUPPLEMENTARY_INCORRECT, gt_conf=100,
                                         with_gt_conf=True),
                create_recall_report_row("truth_probe_6", AlignmentAssessment.PRIMARY_CORRECT, gt_conf=100,
                                         with_gt_conf=True),
            ],
        )

        report = RecallReport([df_1, df_2, df_3], concatenate_dfs_one_by_one_keeping_only_best_mappings=False)
        actual = report.report
        expected = pd.DataFrame(data=[
            create_recall_report_row("truth_probe_1", AlignmentAssessment.UNMAPPED, gt_conf=100, with_gt_conf=True),
            create_recall_report_row("truth_probe_2", AlignmentAssessment.PARTIALLY_MAPPED, gt_conf=100,
                                     with_gt_conf=True),
            create_recall_report_row("truth_probe_3", AlignmentAssessment.PRIMARY_INCORRECT, gt_conf=100,
                                     with_gt_conf=True),
            create_recall_report_row("truth_probe_4", AlignmentAssessment.SECONDARY_INCORRECT, gt_conf=100,
                                     with_gt_conf=True),
            create_recall_report_row("truth_probe_5", AlignmentAssessment.SUPPLEMENTARY_INCORRECT, gt_conf=100,
                                     with_gt_conf=True),
            create_recall_report_row("truth_probe_6", AlignmentAssessment.PRIMARY_CORRECT, gt_conf=100,
                                     with_gt_conf=True),
        ])
        assert_frame_equal(actual, expected, check_dtype=False)


    @patch.object(RecallReport, RecallReport._create_helper_columns.__name__)
    def test____get_id_to_nb_of_allele_sequences_found(self, *mocks):
        contents = StringIO(
"""sample,query_probe_header,PANGENOME_VARIATION_ID,ALLELE_SEQUENCE_ID,good_eval
S1,0,2,0,True
S2,1,0,2,False
S3,2,1,1,True
S4,3,0,2,True
S5,4,1,1,False
S6,5,1,2,False
S7,6,2,1,True
S8,7,1,2,True
S1,8,2,2,True
S1,9,0,2,False
S1,10,2,3,True
S1,11,1,3,False
S1,12,2,4,False
S1,13,2,5,False
S1,14,2,6,False
S1,15,3,0,False
S1,16,3,1,False
""")
        report = RecallReport([pd.read_csv(contents)], False)
        actual=report._get_id_to_nb_of_allele_sequences_found()
        expected=pd.read_csv(StringIO(
"""PANGENOME_VARIATION_ID,NB_OF_ALLELE_SEQUENCE_ID_FOUND
0,1
1,2
2,4
3,0
"""), index_col="PANGENOME_VARIATION_ID")

        assert actual.equals(expected)

    @patch.object(RecallReport, RecallReport._create_helper_columns.__name__)
    def test____get_id_to_nb_of_different_allele_sequences(self, *mocks):
        contents = StringIO(
"""sample,query_probe_header,PANGENOME_VARIATION_ID,NUMBER_OF_DIFFERENT_ALLELE_SEQUENCES,good_eval
S1,0,2,10,True
S2,1,0,1,False
S3,2,1,3,True
S4,3,0,1,True
S5,4,1,3,False
S6,5,1,3,False
S7,6,2,10,True
S8,7,1,3,True
S1,8,2,10,True
S1,9,0,1,False
S1,10,2,10,True
S1,11,1,3,False
S1,12,2,10,False
S1,13,2,10,False
S1,14,2,10,False
S1,15,3,2,False
S1,16,3,2,False
""")
        report = RecallReport([pd.read_csv(contents)], False)
        actual = report._get_id_to_nb_of_different_allele_sequences()
        expected = pd.read_csv(StringIO(
"""PANGENOME_VARIATION_ID,NUMBER_OF_DIFFERENT_ALLELE_SEQUENCES
0,1
1,3
2,10
3,2
"""), index_col="PANGENOME_VARIATION_ID")

        assert actual.equals(expected)

    @patch.object(RecallReport, RecallReport._create_helper_columns.__name__)
    def test____get_id_to_nb_of_samples(self, *mocks):
        contents = StringIO(
            """sample,query_probe_header,PANGENOME_VARIATION_ID,NB_OF_SAMPLES
            S1,0,2,3
            S2,1,0,4
            S3,2,1,10
            """)
        report = RecallReport([pd.read_csv(contents)], False)
        actual = report._get_id_to_nb_of_samples()
        expected = pd.read_csv(StringIO(
            """PANGENOME_VARIATION_ID,NB_OF_SAMPLES
            0,4
            1,10
            2,3
            """), index_col="PANGENOME_VARIATION_ID")

        assert actual.equals(expected)



    @patch.object(RecallReport, RecallReport._create_helper_columns.__name__)
    def test___get_proportion_of_allele_seqs_found_for_each_variant(self, *mocks):
        contents = StringIO(
"""sample,query_probe_header,PANGENOME_VARIATION_ID,ALLELE_SEQUENCE_ID,NUMBER_OF_DIFFERENT_ALLELE_SEQUENCES,good_eval,ALLELE_ID,NUMBER_OF_ALLELES
S1,1,2,0,10,True,0,0
S2,2,0,2,1,False,0,0
S3,3,1,1,3,True,0,0
S4,4,0,2,1,True,0,0
S5,5,1,1,3,False,0,0
S6,6,1,2,3,False,0,0
S7,7,2,1,10,True,0,0
S8,8,1,2,3,True,0,0
S1,9,2,2,10,True,0,0
S1,10,0,2,1,False,0,0
S1,11,2,3,10,True,0,0
S1,12,1,3,3,False,0,0
S1,13,2,4,10,False,0,0
S1,14,2,5,10,False,0,0
S1,15,2,6,10,False,0,0
S1,16,3,0,2,False,0,0
S1,17,3,1,2,False,0,0
""")
        report = RecallReport([pd.read_csv(contents)], False)
        actual = report.get_proportion_of_allele_seqs_found_for_each_variant()
        expected = [1/1, 2/3, 4/10, 0/2]

        assert actual == expected

    @patch.object(RecallReport, RecallReport._create_helper_columns.__name__)
    def test___get_proportion_of_allele_seqs_found_for_each_variant_with_nb_of_samples(self, *mocks):
        contents = StringIO(
            """sample,query_probe_header,PANGENOME_VARIATION_ID,ALLELE_SEQUENCE_ID,NUMBER_OF_DIFFERENT_ALLELE_SEQUENCES,good_eval,ALLELE_ID,NUMBER_OF_ALLELES,NB_OF_SAMPLES
            S1,1,2,0,10,True,0,0,20
            S2,2,0,2,1,False,0,0,1
            S3,3,1,1,3,True,0,0,10
            S4,4,0,2,1,True,0,0,1
            S5,5,1,1,3,False,0,0,10
            S6,6,1,2,3,False,0,0,10
            S7,7,2,1,10,True,0,0,20
            S8,8,1,2,3,True,0,0,10
            S1,9,2,2,10,True,0,0,20
            S1,10,0,2,1,False,0,0,1
            S1,11,2,3,10,True,0,0,20
            S1,12,1,3,3,False,0,0,10
            S1,13,2,4,10,False,0,0,20
            S1,14,2,5,10,False,0,0,20
            S1,15,2,6,10,False,0,0,20
            S1,16,3,0,2,False,0,0,30
            S1,17,3,1,2,False,0,0,30
            """)
        report = RecallReport([pd.read_csv(contents)], False)
        actual = report.get_proportion_of_allele_seqs_found_for_each_variant_with_nb_of_samples()
        expected = pd.read_csv(StringIO(
            """PANGENOME_VARIATION_ID,proportion_of_allele_seqs_found,NB_OF_SAMPLES
            0,1.0,1
            1,0.6666666666666666,10
            2,0.4,20
            3,0.0,30
            """
        ), index_col="PANGENOME_VARIATION_ID")

        assert actual.equals(expected)

    @patch.object(RecallReport, RecallReport._create_helper_columns.__name__)
    def test___get_proportion_of_allele_seqs_found_for_each_variant___duplicated_evaluation_is_disregarded(self, *mocks):
        contents = StringIO(
            """sample,query_probe_header,PANGENOME_VARIATION_ID,ALLELE_SEQUENCE_ID,NUMBER_OF_DIFFERENT_ALLELE_SEQUENCES,good_eval,ALLELE_ID,NUMBER_OF_ALLELES
            S1,1,0,0,5,True,0,0
            S1,2,0,1,5,True,0,0
            S1,3,0,0,5,True,0,0
            S1,4,0,0,5,True,0,0
            S1,5,0,1,5,True,0,0
            """)
        report = RecallReport([pd.read_csv(contents)], False)
        actual = report.get_proportion_of_allele_seqs_found_for_each_variant()
        expected = [2 / 5]

        assert actual == expected

    @patch.object(RecallReport, RecallReport._create_helper_columns.__name__)
    def test___get_proportion_of_alleles_found_for_each_variant(self, *mocks):
        contents = StringIO(
"""sample,query_probe_header,PANGENOME_VARIATION_ID,ALLELE_ID,NUMBER_OF_ALLELES,good_eval,ALLELE_SEQUENCE_ID,NUMBER_OF_DIFFERENT_ALLELE_SEQUENCES
S1,0,2,0,10,True,0,0
S2,1,0,2,1,False,0,0
S3,2,1,1,3,True,0,0
S4,3,0,2,1,True,0,0
S5,4,1,1,3,False,0,0
S6,5,1,2,3,False,0,0
S7,6,2,1,10,True,0,0
S8,7,1,2,3,True,0,0
S1,8,2,2,10,True,0,0
S1,9,0,2,1,False,0,0
S1,10,2,3,10,True,0,0
S1,11,1,3,3,False,0,0
S1,12,2,4,10,False,0,0
S1,13,2,5,10,False,0,0
S1,14,2,6,10,False,0,0
S1,15,3,0,2,False,0,0
S1,16,3,1,2,False,0,0
""")
        report = RecallReport([pd.read_csv(contents)], False)
        actual = report.get_proportion_of_alleles_found_for_each_variant()
        expected = [1/1, 2/3, 4/10, 0/2]

        assert actual == expected


    @patch.object(RecallReport, RecallReport._create_helper_columns.__name__)
    def test___get_proportion_of_alleles_found_for_each_variant___duplicated_evaluation_is_disregarded(self, *mocks):
        contents = StringIO(
            """sample,query_probe_header,PANGENOME_VARIATION_ID,ALLELE_ID,NUMBER_OF_ALLELES,good_eval,ALLELE_SEQUENCE_ID,NUMBER_OF_DIFFERENT_ALLELE_SEQUENCES
            S1,1,0,0,5,True,0,0
            S1,2,0,1,5,True,0,0
            S1,3,0,0,5,True,0,0
            S1,4,0,0,5,True,0,0
            S1,5,0,1,5,True,0,0
            """)
        report = RecallReport([pd.read_csv(contents)], False)
        actual = report.get_proportion_of_alleles_found_for_each_variant()
        expected = [2 / 5]

        assert actual == expected