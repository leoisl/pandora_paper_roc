from evaluate.classification import (
    RecallClassification,
    PrecisionClassification,
    AlignmentAssessment,
)
from evaluate.classifier import Classifier, PrecisionClassifier, RecallClassifier
from tests.common import (
    create_classifier_with_two_entries,
    create_incorrect_supplementary_sam_record,
    create_correct_primary_sam_record,
)


class TestClassifier:
    def test_classify_noAlignmentFileReturnsEmpty(self):
        classifier = Classifier()

        actual = classifier.classify()
        expected = []

        assert actual == expected


class TestRecallClassifier:
    def test_classify(self):
        classifier = create_classifier_with_two_entries(RecallClassifier)

        actual = classifier.classify()
        expected = [
            RecallClassification(record=create_correct_primary_sam_record()),
            RecallClassification(record=create_incorrect_supplementary_sam_record()),
        ]

        assert actual == expected

        expected_classifications = [
            AlignmentAssessment.PRIMARY_CORRECT,
            AlignmentAssessment.SUPPLEMENTARY_INCORRECT,
        ]
        assert [x.assessment() for x in actual] == expected_classifications


class TestPrecisionClassifier:
    def test_classify(self):
        classifier = create_classifier_with_two_entries(PrecisionClassifier)

        actual = classifier.classify()
        expected = [
            PrecisionClassification(record=create_correct_primary_sam_record()),
            PrecisionClassification(record=create_incorrect_supplementary_sam_record()),
        ]

        assert actual == expected

        expected_classifications = [1.0, 0.0]
        assert [x.assessment() for x in actual] == expected_classifications
