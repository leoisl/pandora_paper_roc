from unittest.mock import MagicMock, PropertyMock, patch
from evaluate.filtered_vcf_file import FilteredVCFFile
from .test_vcf_file import build_test_input_and_output
from evaluate.vcf_filters import VCF_Filters
from evaluate.vcf_file import VCFFile

_, sample_to_gene_to_VCFs_all_records = build_test_input_and_output(
    nb_of_samples=3, nb_of_records_in_each_gene=[2, 3, 1]
)

remove_record_filter_mock = MagicMock()
remove_record_filter_mock.record_should_be_filtered_out.return_value = True
keep_record_filter_mock = MagicMock()
keep_record_filter_mock.record_should_be_filtered_out.return_value = False

class TestFilteredVCFFile:
    def test_filter_records_noFiltersReturnsAllRecords(self):
        filters = VCF_Filters()

        actual = FilteredVCFFile._filter_records(sample_to_gene_to_VCFs_all_records, filters)
        expected = sample_to_gene_to_VCFs_all_records

        assert actual == expected

    def test_filter_records_severalFiltersNothingIsFilteredReturnsAllRecords(self):
        filters = VCF_Filters([keep_record_filter_mock]*3)

        actual = FilteredVCFFile._filter_records(sample_to_gene_to_VCFs_all_records, filters)
        expected = sample_to_gene_to_VCFs_all_records

        assert actual == expected

    def test_filter_records_severalFiltersEverythingIsFilteredReturnsNothing(self):
        filters = VCF_Filters([remove_record_filter_mock]*3)

        actual = FilteredVCFFile._filter_records(sample_to_gene_to_VCFs_all_records, filters)
        expected = {}

        assert actual == expected


    def test_filter_records_severalFiltersOneFiltersEverythingReturnsNothing(self):
        filters = VCF_Filters([keep_record_filter_mock]*4 + [remove_record_filter_mock])

        actual = FilteredVCFFile._filter_records(sample_to_gene_to_VCFs_all_records, filters)
        expected = {}

        assert actual == expected

    def test_filter_records_twoFiltersOneFiltersGene0OtherFiltersGene2ReturnsRecordsInGene1(self):
        filter_gene_0_mock = MagicMock()
        filter_gene_0_mock.record_should_be_filtered_out.side_effect = lambda vcf_record : vcf_record.chrom == "gene_0"
        filter_gene_2_mock = MagicMock()
        filter_gene_2_mock.record_should_be_filtered_out.side_effect = lambda vcf_record : vcf_record.chrom == "gene_2"
        filters = VCF_Filters([filter_gene_0_mock, filter_gene_2_mock])

        actual = FilteredVCFFile._filter_records(sample_to_gene_to_VCFs_all_records, filters)
        expected = {"samples_0": {"gene_1": sample_to_gene_to_VCFs_all_records["samples_0"]["gene_1"]},
                    "samples_1": {"gene_1": sample_to_gene_to_VCFs_all_records["samples_1"]["gene_1"]},
                    "samples_2": {"gene_1": sample_to_gene_to_VCFs_all_records["samples_2"]["gene_1"]}}

        assert actual == expected


    @staticmethod
    def fun_test(self, vcf_filepath):
        self._sample_to_gene_to_VCFs = None
        return None
    @patch.object(VCFFile, VCFFile.__init__.__name__, fun_test)
    @patch.object(FilteredVCFFile, FilteredVCFFile._filter_records.__name__, return_value=["test_1", "test_2"])
    def test_constructor_checksIfFiltersAreBeingApplied(self, *mocks):
        filtered_vcf_file = FilteredVCFFile(None, None)

        actual = filtered_vcf_file.sample_to_gene_to_VCFs
        expected = ["test_1", "test_2"]

        assert actual == expected