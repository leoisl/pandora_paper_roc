import tempfile
from io import StringIO

import pandas as pd
import pysam

from evaluate.probe import Interval, ProbeHeader
from evaluate.recall import (
    RecallReporter,
    RecallClassifier,
    RecallClassification,
    RecallCalculator,
    StatisticalClassification,
)
from .test_evaluation import create_sam_header


def create_tmp_sam(contents: str) -> pysam.AlignmentFile:
    with tempfile.NamedTemporaryFile(mode="r+") as tmp:
        tmp.write(contents)
        tmp.truncate()
        return pysam.AlignmentFile(tmp.name)


def create_unmapped_sam_record() -> pysam.AlignedSegment:
    header = create_sam_header(
        "GC00000422_2_SAMPLE=CFT073_POS=603_CALL_INTERVAL=[25,32)_SVTYPE=PH_SNPs_MEAN_FWD_COVG=23_MEAN_REV_COVG=13_GT_CONF=89.5987",
        57,
    )
    record = pysam.AlignedSegment.fromstring(
        ""
        "3_POS=14788_CALL_INTERVAL=[21,22)\t4\t*\t0\t0\t*\t*\t0\t0\tCGCGAAAGCCCTGACCATCTGCACCGTGTCTGACCACATCCGC\t*\tAS:i:0\tXS:i:0",
        header,
    )
    return record


def create_partially_mapped_sam_record() -> pysam.AlignedSegment:
    ref_name = "reference"
    ref_length = 59
    header = create_sam_header(ref_name, ref_length)
    flag = 0
    cigar = "30M38S"
    nm = "NM:i:0"
    md = "MD:Z:30"
    mapq = 60
    pos = 6
    query_name = "INTERVAL=[23,33);"
    sequence = "AAAAAAAAAAAAAAAAAAAAAAACGGCTCGCATAGACACGACGACGACACGTACGATCGATCAGTCAT"
    sam_string = f"{query_name}\t{flag}\t{ref_name}\t{pos}\t{mapq}\t{cigar}\t*\t0\t0\t{sequence}\t*\t{nm}\t{md}\tAS:i:0\tXS:i:0"
    record = pysam.AlignedSegment.fromstring(sam_string, header)
    return record


def create_incorrect_secondary_sam_record() -> pysam.AlignedSegment:
    flag = 256
    cigar = "43M"
    nm = "NM:i:1"
    md = "MD:Z:21T21"
    mapq = 0
    pos = 5
    query_header = ProbeHeader(chrom="3", pos=14788, interval=Interval(21, 22))
    ref_header = ProbeHeader(
        chrom="GC00000422_2",
        sample="CFT073",
        pos=603,
        interval=Interval(25, 32),
        svtype="PH_SNPs",
        mean_fwd_covg=23,
        mean_rev_covg=13,
        gt_conf=89.5987,
    )
    sequence = "CGCGAAAGCCCTGACCATCTGCACCGTGTCTGACCACATCCGC"
    header = create_sam_header(str(ref_header), 57)
    record = pysam.AlignedSegment.fromstring(
        f"{query_header}\t{flag}\t{ref_header}\t{pos}\t{mapq}\t{cigar}\t*\t0\t0\t{sequence}\t*\t{nm}\t{md}\tAS:i:43\tXS:i:32",
        header,
    )
    return record


def create_correct_secondary_sam_record() -> pysam.AlignedSegment:
    flag = 256
    cigar = "43M"
    nm = "NM:i:1"
    md = "MD:Z:19T23"
    mapq = 0
    pos = 5
    query_header = ProbeHeader(chrom="3", pos=14788, interval=Interval(21, 22))
    ref_header = ProbeHeader(
        chrom="GC00000422_2",
        sample="CFT073",
        pos=603,
        interval=Interval(25, 32),
        svtype="PH_SNPs",
        mean_fwd_covg=23,
        mean_rev_covg=13,
        gt_conf=89.5987,
    )
    sequence = "CGCGAAAGCCCTGACCATCTGCACCGTGTCTGACCACATCCGC"
    header = create_sam_header(str(ref_header), 57)
    record = pysam.AlignedSegment.fromstring(
        f"{query_header}\t{flag}\t{ref_header}\t{pos}\t{mapq}\t{cigar}\t*\t0\t0\t{sequence}\t*\t{nm}\t{md}\tAS:i:43\tXS:i:32",
        header,
    )
    return record


def create_incorrect_supplementary_sam_record() -> pysam.AlignedSegment:
    flag = 2048
    cigar = "43M"
    nm = "NM:i:1"
    md = "MD:Z:21T21"
    mapq = 0
    pos = 5
    query_header = ProbeHeader(chrom="3", pos=14788, interval=Interval(21, 22))
    ref_header = ProbeHeader(
        chrom="GC00000422_2",
        sample="CFT073",
        pos=603,
        interval=Interval(25, 32),
        svtype="PH_SNPs",
        mean_fwd_covg=23,
        mean_rev_covg=13,
        gt_conf=89.5987,
    )
    sequence = "CGCGAAAGCCCTGACCATCTGCACCGTGTCTGACCACATCCGC"
    header = create_sam_header(str(ref_header), 57)
    record = pysam.AlignedSegment.fromstring(
        f"{query_header}\t{flag}\t{ref_header}\t{pos}\t{mapq}\t{cigar}\t*\t0\t0\t{sequence}\t*\t{nm}\t{md}\tAS:i:43\tXS:i:32",
        header,
    )
    return record


def create_correct_supplementary_sam_record() -> pysam.AlignedSegment:
    flag = 2048
    cigar = "43M"
    nm = "NM:i:1"
    md = "MD:Z:19T23"
    mapq = 0
    pos = 5
    query_header = ProbeHeader(chrom="3", pos=14788, interval=Interval(21, 22))
    ref_header = ProbeHeader(
        chrom="GC00000422_2",
        sample="CFT073",
        pos=603,
        interval=Interval(25, 32),
        svtype="PH_SNPs",
        mean_fwd_covg=23,
        mean_rev_covg=13,
        gt_conf=89.5987,
    )
    sequence = "CGCGAAAGCCCTGACCATCTGCACCGTGTCTGACCACATCCGC"
    header = create_sam_header(str(ref_header), 57)
    record = pysam.AlignedSegment.fromstring(
        f"{query_header}\t{flag}\t{ref_header}\t{pos}\t{mapq}\t{cigar}\t*\t0\t0\t{sequence}\t*\t{nm}\t{md}\tAS:i:43\tXS:i:32",
        header,
    )
    return record


def create_correct_primary_sam_record() -> pysam.AlignedSegment:
    flag = 0
    cigar = "56M"
    nm = "NM:i:0"
    md = "MD:Z:56"
    mapq = 60
    pos = 1
    query_header = ProbeHeader(interval=Interval(12, 17))
    ref_header = ProbeHeader(
        chrom="GC00000422_2",
        sample="CFT073",
        pos=603,
        interval=Interval(25, 32),
        svtype="PH_SNPs",
        mean_fwd_covg=23,
        mean_rev_covg=13,
        gt_conf=89.5987,
    )
    sequence = "AAAAAAAAAAACGGCTCGCATAGACACGACGACGACACGTACGATCGATCAGTCAT"
    header = create_sam_header(str(ref_header), 64)
    record = pysam.AlignedSegment.fromstring(
        f"{query_header}\t{flag}\t{ref_header}\t{pos}\t{mapq}\t{cigar}\t*\t0\t0\t{sequence}\t*\t{nm}\t{md}\tAS:i:43\tXS:i:32",
        header,
    )
    return record


def create_incorrect_primary_sam_record() -> pysam.AlignedSegment:
    flag = 0
    cigar = "56M"
    nm = "NM:i:1"
    md = "MD:Z:12T43"
    mapq = 60
    pos = 1
    query_header = ProbeHeader(interval=Interval(12, 13))
    ref_header = ProbeHeader(
        chrom="GC00000422_2",
        sample="CFT073",
        pos=603,
        interval=Interval(25, 32),
        svtype="PH_SNPs",
        mean_fwd_covg=23,
        mean_rev_covg=13,
        gt_conf=89.5987,
    )
    sequence = "AAAAAAAAAAACGGCTCGCATAGACACGACGACGACACGTACGATCGATCAGTCAT"
    header = create_sam_header(str(ref_header), 64)
    record = pysam.AlignedSegment.fromstring(
        f"{query_header}\t{flag}\t{ref_header}\t{pos}\t{mapq}\t{cigar}\t*\t0\t0\t{sequence}\t*\t{nm}\t{md}\tAS:i:43\tXS:i:32",
        header,
    )
    return record


def create_recall_classifier_with_two_entries() -> RecallClassifier:
    flag = 0
    cigar = "56M"
    nm = "NM:i:0"
    md = "MD:Z:56"
    mapq = 60
    pos = 1
    query_header = ProbeHeader(interval=Interval(12, 17))
    sequence = "AAAAAAAAAAACGGCTCGCATAGACACGACGACGACACGTACGATCGATCAGTCAT"
    ref_header = ProbeHeader(
        chrom="GC00000422_2",
        sample="CFT073",
        pos=603,
        interval=Interval(25, 32),
        svtype="PH_SNPs",
        mean_fwd_covg=23,
        mean_rev_covg=13,
        gt_conf=89.5987,
    )
    header = create_sam_header(str(ref_header), 64)
    contents = str(header) + "\n"
    record1 = pysam.AlignedSegment.fromstring(
        f"{query_header}\t{flag}\t{ref_header}\t{pos}\t{mapq}\t{cigar}\t*\t0\t0\t{sequence}\t*\t{nm}\t{md}\tAS:i:43\tXS:i:32",
        header,
    )
    contents += record1.to_string() + "\n"

    flag = 2048
    cigar = "43M"
    nm = "NM:i:1"
    md = "MD:Z:21T21"
    mapq = 0
    pos = 5
    query_header = ProbeHeader(chrom="3", pos=14788, interval=Interval(21, 22))
    sequence = "CGCGAAAGCCCTGACCATCTGCACCGTGTCTGACCACATCCGC"
    header = create_sam_header(str(ref_header), 57)
    record2 = pysam.AlignedSegment.fromstring(
        f"{query_header}\t{flag}\t{ref_header}\t{pos}\t{mapq}\t{cigar}\t*\t0\t0\t{sequence}\t*\t{nm}\t{md}\tAS:i:43\tXS:i:32",
        header,
    )
    contents += record2.to_string() + "\n"
    sam = create_tmp_sam(contents)
    return RecallClassifier(sam)


class TestRecallClassifier:
    def test_classify_noAlignmentFileReturnsEmpty(self):
        classifier = RecallClassifier()

        actual = classifier.classify()
        expected = []

        assert actual == expected

    def test_classify(self):
        classifier = create_recall_classifier_with_two_entries()

        actual = classifier.classify()
        expected = [
            RecallClassification(record=create_correct_primary_sam_record()),
            RecallClassification(record=create_incorrect_supplementary_sam_record()),
        ]

        assert actual == expected

        expected_classifications = ["correct", "supplementary_incorrect"]
        assert [x.assessment() for x in actual] == expected_classifications


class TestRecallReporter:
    def test_generateReport_noClassifierReturnsEmpty(self):
        sample = "sample"
        classifier = RecallClassifier(name=sample)
        reporter = RecallReporter(classifiers=[classifier])

        actual = reporter.generate_report()
        expected = pd.DataFrame(
            [],
            columns=[
                "sample",
                "truth_probe_header",
                "vcf_probe_header",
                "classification",
            ],
        )

        assert actual.equals(expected)

    def test_generateReport_twoClassificationsReturnsDataframeWithTwoEntries(self):
        classifier = create_recall_classifier_with_two_entries()
        sample = "sample"
        classifier.name = sample
        reporter = RecallReporter(classifiers=[classifier])

        actual = reporter.generate_report()
        expected_data = []
        for assessment, record in [
            ("correct", create_correct_primary_sam_record()),
            ("supplementary_incorrect", create_incorrect_supplementary_sam_record()),
        ]:
            expected_data.append(
                [sample, record.query_name, record.reference_name, assessment]
            )
        expected = pd.DataFrame(
            expected_data,
            columns=[
                "sample",
                "truth_probe_header",
                "vcf_probe_header",
                "classification",
            ],
        )

        assert actual.equals(expected)

    def test_save_emptyReporterReturnsHeadersOnly(self):
        delim = "\t"
        reporter = RecallReporter(classifiers=[RecallClassifier()], delim=delim)
        fh = StringIO(newline="")
        reporter.save(fh)

        fh.seek(0)
        actual = fh.read()
        expected = delim.join(reporter.columns) + "\n"

        assert actual == expected

    def test_save_reporterWithTwoClassificationsWritesHeadersAndTwoRows(self):
        primary_correct_record = create_correct_primary_sam_record()
        suppl_incorrect_record = create_incorrect_supplementary_sam_record()
        delim = "\t"
        classifier = create_recall_classifier_with_two_entries()
        sample = "sample"
        classifier.name = sample
        reporter = RecallReporter(classifiers=[classifier], delim=delim)

        fh = StringIO(newline="")
        reporter.save(fh)

        fh.seek(0)
        actual = fh.read()
        expected_data = []
        for assessment, record in [
            ("correct", primary_correct_record),
            ("supplementary_incorrect", suppl_incorrect_record),
        ]:
            expected_data.append(
                [sample, record.query_name, record.reference_name, assessment]
            )
        expected = StringIO(newline="")
        pd.DataFrame(
            expected_data,
            columns=[
                "sample",
                "truth_probe_header",
                "vcf_probe_header",
                "classification",
            ],
        ).to_csv(expected, sep=delim, header=True, index=False)
        expected.seek(0)
        expected = expected.read()

        assert actual == expected

    def test_save_reporterWithTwoClassificationsWritesHeadersAndTwoRowsWithCommaDelim(
        self
    ):
        primary_correct_record = create_correct_primary_sam_record()
        suppl_incorrect_record = create_incorrect_supplementary_sam_record()
        delim = ","
        classifier = create_recall_classifier_with_two_entries()
        sample = "sample"
        classifier.name = sample
        reporter = RecallReporter(classifiers=[classifier], delim=delim)

        fh = StringIO(newline="")
        reporter.save(fh)

        fh.seek(0)
        actual = fh.read()
        expected_data = []
        for assessment, record in [
            ("correct", primary_correct_record),
            ("supplementary_incorrect", suppl_incorrect_record),
        ]:
            expected_data.append(
                [sample, record.query_name, record.reference_name, assessment]
            )
        expected = StringIO(newline="")
        pd.DataFrame(
            expected_data,
            columns=[
                "sample",
                "truth_probe_header",
                "vcf_probe_header",
                "classification",
            ],
        ).to_csv(expected, sep=delim, header=True, index=False)
        expected.seek(0)
        expected = expected.read()

        assert actual == expected

    def test_save_reporterWithTwoClassifiersWritesTwoSamplesWithTwoRows(self):
        primary_correct_record = create_correct_primary_sam_record()
        suppl_incorrect_record = create_incorrect_supplementary_sam_record()
        delim = ","
        classifier1 = create_recall_classifier_with_two_entries()
        sample = "sample"
        classifier1.name = sample
        classifier2 = create_recall_classifier_with_two_entries()
        sample2 = "sample2"
        classifier2.name = sample2

        reporter = RecallReporter(classifiers=[classifier1, classifier2], delim=delim)

        fh = StringIO(newline="")
        reporter.save(fh)

        fh.seek(0)
        actual = fh.read()
        expected_data = []
        for s in [sample, sample2]:
            for assessment, record in [
                ("correct", primary_correct_record),
                ("supplementary_incorrect", suppl_incorrect_record),
            ]:
                expected_data.append(
                    [s, record.query_name, record.reference_name, assessment]
                )
        expected = StringIO(newline="")
        pd.DataFrame(
            expected_data,
            columns=[
                "sample",
                "truth_probe_header",
                "vcf_probe_header",
                "classification",
            ],
        ).to_csv(expected, sep=delim, header=True, index=False)
        expected.seek(0)
        expected = expected.read()

        assert actual == expected


def create_report_row(
    classification: str, gt_conf: float = 0, sample: str = "sample1"
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
        classification = "incorrect"
        gt_conf = 10
        row = create_report_row(classification, gt_conf)
        threshold = 5

        actual = RecallCalculator.statistical_classification(
            row, conf_threshold=threshold
        )
        expected = StatisticalClassification.FALSE_POSITIVE

        assert actual == expected

    def test_statisticalClassification_incorrectBelowConfReturnsFalseNegative(self):
        classification = "incorrect"
        gt_conf = 10
        row = create_report_row(classification, gt_conf)
        threshold = 50

        actual = RecallCalculator.statistical_classification(
            row, conf_threshold=threshold
        )
        expected = StatisticalClassification.FALSE_NEGATIVE

        assert actual == expected

    def test_statisticalClassification_correctAboveConfReturnsTruePositive(self):
        classification = "correct"
        gt_conf = 10
        row = create_report_row(classification, gt_conf)
        threshold = 5

        actual = RecallCalculator.statistical_classification(
            row, conf_threshold=threshold
        )
        expected = StatisticalClassification.TRUE_POSITIVE

        assert actual == expected

    def test_statisticalClassification_correctBelowConfReturnsFalseNegative(self):
        classification = "correct"
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
        classification = "secondary_incorrect"
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
        classification = "secondary_incorrect"
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
        classification = "secondary_correct"
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
        classification = "secondary_correct"
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
        classification = "supplementary_incorrect"
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
        classification = "supplementary_incorrect"
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
        classification = "supplementary_correct"
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
        classification = "supplementary_correct"
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
                create_report_row("unmapped", gt_conf=100),
                create_report_row("unmapped", gt_conf=100),
                create_report_row("correct", gt_conf=10),
                create_report_row("incorrect", gt_conf=100),
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
                create_report_row("correct", gt_conf=100),
                create_report_row("correct", gt_conf=100),
                create_report_row("incorrect", gt_conf=100),
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
                create_report_row("correct", gt_conf=10),
                create_report_row("correct", gt_conf=100),
                create_report_row("correct", gt_conf=100),
                create_report_row("unmapped", gt_conf=100),
                create_report_row("supplementary_incorrect", gt_conf=10),
                create_report_row("secondary_correct", gt_conf=100),
                create_report_row("incorrect", gt_conf=100),
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
                create_report_row("incorrect", gt_conf=100),
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
                create_report_row("correct", gt_conf=10),
                create_report_row("correct", gt_conf=100),
                create_report_row("correct", gt_conf=100),
                create_report_row("unmapped", gt_conf=100),
                create_report_row("supplementary_incorrect", gt_conf=10),
                create_report_row("secondary_correct", gt_conf=100),
                create_report_row("incorrect", gt_conf=100),
            ],
            columns=columns,
        )
        report2 = pd.DataFrame(
            data=[
                create_report_row("correct", gt_conf=100),
                create_report_row("unmapped", gt_conf=100),
                create_report_row("incorrect", gt_conf=100),
            ],
            columns=columns,
        )
        calculator = RecallCalculator([report1, report2])
        threshold = 60

        actual = calculator.calculate_recall(conf_threshold=threshold)
        expected = 0.5

        assert actual == expected
