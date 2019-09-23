from evaluate.vcf_filters import VCF_Filters
import pytest
from unittest.mock import MagicMock


@pytest.fixture
def vcf_filter_that_always_filter_out():
    vcf_filter = MagicMock(record_should_be_filtered_out=MagicMock(return_value=True))
    return vcf_filter


@pytest.fixture
def vcf_filter_that_never_filter_out():
    vcf_filter = MagicMock(record_should_be_filtered_out=MagicMock(return_value=False))
    return vcf_filter

class Test_VCF_Filters:
    def test_recordShouldBeFilteredOut_recordDoesNotPassAnyFilterReturnsTrue(self, vcf_filter_that_always_filter_out):
        vcf_filters = VCF_Filters([vcf_filter_that_always_filter_out]*3)
        record = MagicMock()

        actual = vcf_filters.record_should_be_filtered_out(record)
        expected = True

        assert actual == expected

    def test_recordShouldNotBeFilteredOut_recordPassAllFiltersReturnsFalse(self, vcf_filter_that_never_filter_out):
        vcf_filters = VCF_Filters([vcf_filter_that_never_filter_out]*3)
        record = MagicMock()

        actual = vcf_filters.record_should_be_filtered_out(record)
        expected = False

        assert actual == expected

    def test_recordShouldBeFilteredOut_recordPassAllFiltersButOneReturnsTrue(self, vcf_filter_that_never_filter_out, vcf_filter_that_always_filter_out):
        vcf_filters = VCF_Filters([vcf_filter_that_never_filter_out]*3 + [vcf_filter_that_always_filter_out])
        record = MagicMock()

        actual = vcf_filters.record_should_be_filtered_out(record)
        expected = True

        assert actual == expected