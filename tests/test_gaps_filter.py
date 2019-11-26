from unittest.mock import MagicMock
from evaluate.gaps_filter import GapsFilter


class TestGapsFilter:
    def test_recordShouldBeFilteredOut_justBeforeThresholdReturnsFalse(self):
        gaps_filter = GapsFilter(0.75)
        record_mock = MagicMock(gaps=0.74)

        actual = gaps_filter.record_should_be_filtered_out(record_mock)
        expected = False

        assert actual == expected

    def test_recordShouldBeFilteredOut_equalsThresholdReturnsFalse(self):
        gaps_filter = GapsFilter(0.75)
        record_mock = MagicMock(gaps=0.75)

        actual = gaps_filter.record_should_be_filtered_out(record_mock)
        expected = False

        assert actual == expected

    def test_recordShouldBeFilteredOut_justAfterThresholdReturnsTrue(self):
        gaps_filter = GapsFilter(0.75)
        record_mock = MagicMock(gaps=0.76)

        actual = gaps_filter.record_should_be_filtered_out(record_mock)
        expected = True

        assert actual == expected
