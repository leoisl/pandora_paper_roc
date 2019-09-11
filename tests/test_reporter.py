from io import StringIO

import pandas as pd

from evaluate.classification import AlignmentAssessment
from evaluate.classifier import RecallClassifier
from evaluate.reporter import RecallReporter
from tests.common import (
    create_classifier_with_two_entries,
    create_correct_primary_sam_record,
    create_incorrect_supplementary_sam_record,
)


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
                "query_probe_header",
                "ref_probe_header",
                "classification",
            ],
        )

        assert actual.equals(expected)

    def test_generateReport_twoClassificationsReturnsDataframeWithTwoEntries(self):
        classifier = create_classifier_with_two_entries(RecallClassifier)
        sample = "sample"
        classifier.name = sample
        reporter = RecallReporter(classifiers=[classifier])

        actual = reporter.generate_report()
        expected_data = []
        for assessment, record in [
            (AlignmentAssessment.PRIMARY_CORRECT, create_correct_primary_sam_record()),
            (
                AlignmentAssessment.SUPPLEMENTARY_INCORRECT,
                create_incorrect_supplementary_sam_record(),
            ),
        ]:
            expected_data.append(
                [sample, record.query_name, record.reference_name, assessment]
            )
        expected = pd.DataFrame(
            expected_data,
            columns=[
                "sample",
                "query_probe_header",
                "ref_probe_header",
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
        classifier = create_classifier_with_two_entries(RecallClassifier)
        sample = "sample"
        classifier.name = sample
        reporter = RecallReporter(classifiers=[classifier], delim=delim)

        fh = StringIO(newline="")
        reporter.save(fh)

        fh.seek(0)
        actual = fh.read()
        expected_data = []
        for assessment, record in [
            (AlignmentAssessment.PRIMARY_CORRECT, primary_correct_record),
            (AlignmentAssessment.SUPPLEMENTARY_INCORRECT, suppl_incorrect_record),
        ]:
            expected_data.append(
                [sample, record.query_name, record.reference_name, assessment]
            )
        expected = StringIO(newline="")
        pd.DataFrame(
            expected_data,
            columns=[
                "sample",
                "query_probe_header",
                "ref_probe_header",
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
        classifier = create_classifier_with_two_entries(RecallClassifier)
        sample = "sample"
        classifier.name = sample
        reporter = RecallReporter(classifiers=[classifier], delim=delim)

        fh = StringIO(newline="")
        reporter.save(fh)

        fh.seek(0)
        actual = fh.read()
        expected_data = []
        for assessment, record in [
            (AlignmentAssessment.PRIMARY_CORRECT, primary_correct_record),
            (AlignmentAssessment.SUPPLEMENTARY_INCORRECT, suppl_incorrect_record),
        ]:
            expected_data.append(
                [sample, record.query_name, record.reference_name, assessment]
            )
        expected = StringIO(newline="")
        pd.DataFrame(
            expected_data,
            columns=[
                "sample",
                "query_probe_header",
                "ref_probe_header",
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
        classifier1 = create_classifier_with_two_entries(RecallClassifier)
        sample = "sample"
        classifier1.name = sample
        classifier2 = create_classifier_with_two_entries(RecallClassifier)
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
                (AlignmentAssessment.PRIMARY_CORRECT, primary_correct_record),
                (AlignmentAssessment.SUPPLEMENTARY_INCORRECT, suppl_incorrect_record),
            ]:
                expected_data.append(
                    [s, record.query_name, record.reference_name, assessment]
                )
        expected = StringIO(newline="")
        pd.DataFrame(
            expected_data,
            columns=[
                "sample",
                "query_probe_header",
                "ref_probe_header",
                "classification",
            ],
        ).to_csv(expected, sep=delim, header=True, index=False)
        expected.seek(0)
        expected = expected.read()

        assert actual == expected
