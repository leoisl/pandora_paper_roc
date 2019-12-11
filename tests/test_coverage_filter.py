from unittest.mock import MagicMock
from evaluate.coverage_filter import CoverageFilter


class TestCoverageFilter:
    def test_recordShouldBeFilteredOut_coverageBelowThresholdReturnsTrue(self):
        record = MagicMock(coverage=5)
        coverage_filter = CoverageFilter(10.0)

        actual = coverage_filter.record_should_be_filtered_out(record)
        expected = True

        assert actual == expected

    def test_recordShouldBeFilteredOut_coverageAboveThresholdReturnsFalse(self):
        record = MagicMock(coverage=15)
        coverage_filter = CoverageFilter(10.0)

        actual = coverage_filter.record_should_be_filtered_out(record)
        expected = False

        assert actual == expected

    def test_recordShouldBeFilteredOut_coverageEqualsThresholdReturnsFalse(self):
        record = MagicMock(coverage=10)
        coverage_filter = CoverageFilter(10.0)

        actual = coverage_filter.record_should_be_filtered_out(record)
        expected = False

        assert actual == expected
