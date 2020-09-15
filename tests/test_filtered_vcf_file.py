from unittest.mock import Mock, patch, PropertyMock
from evaluate.filtered_vcf_file import FilteredVCFFile
from evaluate.vcf_filters import VCF_Filters
from evaluate.vcf_file import VCFFile
import pytest

@pytest.fixture
def remove_record_filter_mock():
    remove_record_filter_mock = Mock()
    remove_record_filter_mock.record_should_be_filtered_out.return_value = True
    return remove_record_filter_mock

@pytest.fixture
def keep_record_filter_mock():
    keep_record_filter_mock = Mock()
    keep_record_filter_mock.record_should_be_filtered_out.return_value = False
    return keep_record_filter_mock

@pytest.fixture
def sample_to_gene_to_VCFs_all_records():
    sample_to_gene_to_VCFs_all_records = {
        "sample_1": {
            "gene_0": [Mock(), Mock()],
            "gene_1": [Mock()],
            "gene_2": [Mock(), Mock(), Mock()],
        },
        "sample_2": {
            "gene_0": [Mock()],
            "gene_1": [Mock(), Mock(), Mock()],
            "gene_2": [Mock(), Mock(), Mock(), Mock(), Mock(), Mock()],
        }
    }
    return sample_to_gene_to_VCFs_all_records

class TestFilteredVCFFile:
    filter_records_mock_return_value = Mock()
    sample_to_gene_to_VCFs_mock = Mock()

    @patch.object(FilteredVCFFile, "sample_to_gene_to_VCFs", new_callable=PropertyMock, return_value=sample_to_gene_to_VCFs_mock)
    @patch.object(VCFFile, VCFFile.__init__.__name__)
    @patch.object(FilteredVCFFile, FilteredVCFFile._filter_records.__name__, return_value=filter_records_mock_return_value)
    def test___constructor(self, filter_records_mock, VCFFile_init_mock, sample_to_gene_to_VCFs_property_mock):
        pysam_variant_file_mock = Mock()
        filters_mock = Mock()
        VCF_creator_method_mock = Mock()

        filtered_vcf_file = FilteredVCFFile(pysam_variant_file_mock, filters_mock, VCF_creator_method_mock)

        VCFFile_init_mock.assert_called_once_with(filtered_vcf_file, pysam_variant_file_mock, VCF_creator_method_mock)
        filter_records_mock.assert_called_once_with(TestFilteredVCFFile.sample_to_gene_to_VCFs_mock, filters_mock)


    def test___constructor___empty_medaka_file(self):
        import pysam
        from evaluate.vcf import VCFFactory
        with pysam.VariantFile("tests/test_cases/sample_medaka_empty_vcf.expected.vcf") as pysam_variant_file:
            filtered_vcf_file = FilteredVCFFile(pysam_variant_file=pysam_variant_file, filters=[],
                                                VCF_creator_method=VCFFactory.create_Medaka_VCF_from_VariantRecord_and_Sample)
        assert filtered_vcf_file._sample_to_gene_to_VCFs == {}


    def test_filter_records_noFiltersReturnsAllRecords(self, sample_to_gene_to_VCFs_all_records):
        filters = VCF_Filters()

        actual = FilteredVCFFile._filter_records(sample_to_gene_to_VCFs_all_records, filters)
        expected = sample_to_gene_to_VCFs_all_records

        assert actual == expected

    def test_filter_records_severalFiltersNothingIsFilteredReturnsAllRecords(self,
                                                                             sample_to_gene_to_VCFs_all_records,
                                                                             keep_record_filter_mock):
        filters = VCF_Filters([keep_record_filter_mock]*3)

        actual = FilteredVCFFile._filter_records(sample_to_gene_to_VCFs_all_records, filters)
        expected = sample_to_gene_to_VCFs_all_records

        assert actual == expected

    def test_filter_records_severalFiltersEverythingIsFilteredReturnsNothing(self,
                                                                             sample_to_gene_to_VCFs_all_records,
                                                                             remove_record_filter_mock):
        filters = VCF_Filters([remove_record_filter_mock]*3)

        actual = FilteredVCFFile._filter_records(sample_to_gene_to_VCFs_all_records, filters)
        expected = {}

        assert actual == expected


    def test_filter_records_severalFiltersOneFiltersEverythingReturnsNothing(self,
                                                                             sample_to_gene_to_VCFs_all_records,
                                                                             keep_record_filter_mock,
                                                                             remove_record_filter_mock):
        filters = VCF_Filters([keep_record_filter_mock]*4 + [remove_record_filter_mock])

        actual = FilteredVCFFile._filter_records(sample_to_gene_to_VCFs_all_records, filters)
        expected = {}

        assert actual == expected

    def test_filter_records_twoFiltersOneFiltersGene0OtherFiltersGene2ReturnsRecordsInGene1(self,
                                                                                            sample_to_gene_to_VCFs_all_records):
        filter_gene_0_mock = Mock()
        filter_gene_0_mock.record_should_be_filtered_out.side_effect = \
            lambda vcf_record : \
            vcf_record in sample_to_gene_to_VCFs_all_records["sample_1"]["gene_0"] or \
            vcf_record in sample_to_gene_to_VCFs_all_records["sample_2"]["gene_0"]
        filter_gene_2_mock = Mock()
        filter_gene_2_mock.record_should_be_filtered_out.side_effect = \
            lambda vcf_record: \
                vcf_record in sample_to_gene_to_VCFs_all_records["sample_1"]["gene_2"] or \
                vcf_record in sample_to_gene_to_VCFs_all_records["sample_2"]["gene_2"]
        filters = VCF_Filters([filter_gene_0_mock, filter_gene_2_mock])

        actual = FilteredVCFFile._filter_records(sample_to_gene_to_VCFs_all_records, filters)
        expected = {
            "sample_1": {
                "gene_1": sample_to_gene_to_VCFs_all_records["sample_1"]["gene_1"],
            },
            "sample_2": {
                "gene_1": sample_to_gene_to_VCFs_all_records["sample_2"]["gene_1"],
            }
        }

        assert actual == expected