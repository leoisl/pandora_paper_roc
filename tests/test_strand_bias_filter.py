from unittest.mock import MagicMock
from evaluate.strand_bias_filter import StrandBiasFilter


class TestStrandBiasFilter:
    def test_recordShouldBeFilteredOut_meanCoverageZeroReturnsTrue(self):
        mocked_vcf = MagicMock(mean_coverage=0)
        strand_bias_filter = StrandBiasFilter(0.1)

        actual = strand_bias_filter.record_should_be_filtered_out(mocked_vcf)
        expected = True

        assert actual == expected

    def test_recordShouldBeFilteredOut_strandRatioLessThanThresholdReturnsTrue(self):
        mocked_vcf = MagicMock(mean_coverage=100, mean_coverage_forward=9)
        strand_bias_filter = StrandBiasFilter(0.1)

        actual = strand_bias_filter.record_should_be_filtered_out(mocked_vcf)
        expected = True

        assert actual == expected

    def test_recordShouldBeFilteredOut_strandRatioEqualsThresholdReturnsTrue(self):
        mocked_vcf = MagicMock(mean_coverage=100, mean_coverage_forward=10)
        strand_bias_filter = StrandBiasFilter(0.1)

        actual = strand_bias_filter.record_should_be_filtered_out(mocked_vcf)
        expected = True

        assert actual == expected

    def test_recordShouldBeFilteredOut_strandRatioMoreThanUpperThresholdReturnsTrue(
        self
    ):
        mocked_vcf = MagicMock(mean_coverage=100, mean_coverage_forward=91)
        strand_bias_filter = StrandBiasFilter(0.1)

        actual = strand_bias_filter.record_should_be_filtered_out(mocked_vcf)
        expected = True

        assert actual == expected

    def test_recordShouldBeFilteredOut_strandRatioEqualsUpperThresholdReturnsTrue(self):
        mocked_vcf = MagicMock(mean_coverage=100, mean_coverage_forward=90)
        strand_bias_filter = StrandBiasFilter(0.1)

        actual = strand_bias_filter.record_should_be_filtered_out(mocked_vcf)
        expected = True

        assert actual == expected

    def test_recordShouldBeFilteredOut_strandRatioJustAfterThresholdReturnsFalse(self):
        mocked_vcf = MagicMock(mean_coverage=100, mean_coverage_forward=11)
        strand_bias_filter = StrandBiasFilter(0.1)

        actual = strand_bias_filter.record_should_be_filtered_out(mocked_vcf)
        expected = False

        assert actual == expected

    def test_recordShouldBeFilteredOut_strandRatioJustBeforeUpperThresholdReturnsFalse(self):
        mocked_vcf = MagicMock(mean_coverage=100, mean_coverage_forward=89)
        strand_bias_filter = StrandBiasFilter(0.1)

        actual = strand_bias_filter.record_should_be_filtered_out(mocked_vcf)
        expected = False

        assert actual == expected