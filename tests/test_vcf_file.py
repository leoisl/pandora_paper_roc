from unittest.mock import patch, Mock, PropertyMock
from evaluate.vcf_file import VCFFile
from evaluate.vcf import NullVCFError, BuggedVCFError, VCFFactory
import pytest


@pytest.fixture
def pysam_variant_record_mock_that_maps_to_chrom_1_and_one_sample():
    pysam_variant_record_mock_that_maps_to_chrom_1_and_one_sample = Mock()
    pysam_variant_record_mock_that_maps_to_chrom_1_and_one_sample.samples = ["sample_1"]
    pysam_variant_record_mock_that_maps_to_chrom_1_and_one_sample.chrom = "chrom_1"
    return pysam_variant_record_mock_that_maps_to_chrom_1_and_one_sample

@pytest.fixture
def pysam_variant_record_mock_that_maps_to_chrom_2_and_one_sample():
    pysam_variant_record_mock_that_maps_to_chrom_2_and_one_sample = Mock()
    pysam_variant_record_mock_that_maps_to_chrom_2_and_one_sample.samples = ["sample_1"]
    pysam_variant_record_mock_that_maps_to_chrom_2_and_one_sample.chrom = "chrom_2"
    return pysam_variant_record_mock_that_maps_to_chrom_2_and_one_sample

@pytest.fixture
def pysam_variant_record_mock_that_maps_to_chrom_1_and_two_samples():
    pysam_variant_record_mock_that_maps_to_chrom_1_and_two_samples = Mock()
    pysam_variant_record_mock_that_maps_to_chrom_1_and_two_samples.samples = ["sample_1", "sample_2"]
    pysam_variant_record_mock_that_maps_to_chrom_1_and_two_samples.chrom = "chrom_1"
    return pysam_variant_record_mock_that_maps_to_chrom_1_and_two_samples

@pytest.fixture
def vcf_record_1_mock():
    vcf_record_1_mock = Mock()
    return vcf_record_1_mock

@pytest.fixture
def vcf_record_2_mock():
    vcf_record_2_mock = Mock()
    return vcf_record_2_mock

@pytest.fixture
def vcf_record_3_mock():
    vcf_record_3_mock = Mock()
    return vcf_record_3_mock


sample_with_some_genes = {
    "sample_1": {
        "gene_1": [1, 2, 3, 4],
        "gene_2": [5, 6],
    },
    "sample_2": {
        "gene_1": [7, 8, 9],
        "gene_2": [10],
    }
}

def chrom_1_raises_NullVCFError_others_are_fine(pysam_variant_record, sample):
    if pysam_variant_record.chrom == "chrom_1":
        raise NullVCFError()
    else:
        return vcf_record_2_mock


class Test_VCFFile:
    def test___constructor___no_records_in_VCF_returns_nothing(self):
        vcf_file = VCFFile([], VCFFactory.create_Pandora_VCF_from_VariantRecord_and_Sample)
        actual = vcf_file.sample_to_gene_to_VCFs

        expected = {}
        assert actual == expected

    @patch.object(VCFFactory, VCFFactory.create_Pandora_VCF_from_VariantRecord_and_Sample.__name__, return_value=vcf_record_1_mock)
    def test___constructor___one_record_in_one_sample_and_one_gene(self, from_VariantRecord_and_Sample_Mock,
                                                                   pysam_variant_record_mock_that_maps_to_chrom_1_and_one_sample):
        vcf_file = VCFFile([pysam_variant_record_mock_that_maps_to_chrom_1_and_one_sample], VCFFactory.create_Pandora_VCF_from_VariantRecord_and_Sample)
        actual = vcf_file.sample_to_gene_to_VCFs

        expected = {"sample_1": {"chrom_1": [vcf_record_1_mock]}}
        assert actual == expected

    @patch.object(VCFFactory, VCFFactory.create_Pandora_VCF_from_VariantRecord_and_Sample.__name__,
                  return_value=vcf_record_1_mock)
    def test___constructor___one_record_in_two_samples_and_one_gene(self, from_VariantRecord_and_Sample_Mock,
                                                                    pysam_variant_record_mock_that_maps_to_chrom_1_and_two_samples):
        vcf_file = VCFFile([pysam_variant_record_mock_that_maps_to_chrom_1_and_two_samples], VCFFactory.create_Pandora_VCF_from_VariantRecord_and_Sample)
        actual = vcf_file.sample_to_gene_to_VCFs

        expected = {"sample_1": {"chrom_1": [vcf_record_1_mock]},
                    "sample_2": {"chrom_1": [vcf_record_1_mock]}}
        assert actual == expected

    @patch.object(VCFFactory, VCFFactory.create_Pandora_VCF_from_VariantRecord_and_Sample.__name__,
                  side_effect=[vcf_record_1_mock, vcf_record_2_mock])
    def test___constructor___two_records_in_one_sample_and_two_genes(self, from_VariantRecord_and_Sample_Mock,
                                                                     pysam_variant_record_mock_that_maps_to_chrom_1_and_one_sample,
                                                                     pysam_variant_record_mock_that_maps_to_chrom_2_and_one_sample):
        vcf_file = VCFFile([pysam_variant_record_mock_that_maps_to_chrom_1_and_one_sample,
                            pysam_variant_record_mock_that_maps_to_chrom_2_and_one_sample], VCFFactory.create_Pandora_VCF_from_VariantRecord_and_Sample)
        actual = vcf_file.sample_to_gene_to_VCFs

        expected = {"sample_1": {"chrom_1": [vcf_record_1_mock], "chrom_2": [vcf_record_2_mock]}}
        assert actual == expected

    @patch.object(VCFFactory, VCFFactory.create_Pandora_VCF_from_VariantRecord_and_Sample.__name__,
                  side_effect=[vcf_record_1_mock, vcf_record_2_mock])
    def test___constructor___two_records_in_one_sample_and_one_gene(self, from_VariantRecord_and_Sample_Mock,
                                                                    pysam_variant_record_mock_that_maps_to_chrom_1_and_one_sample):
        vcf_file = VCFFile([pysam_variant_record_mock_that_maps_to_chrom_1_and_one_sample,
                            pysam_variant_record_mock_that_maps_to_chrom_1_and_one_sample], VCFFactory.create_Pandora_VCF_from_VariantRecord_and_Sample)
        actual = vcf_file.sample_to_gene_to_VCFs

        expected = {"sample_1": {"chrom_1": [vcf_record_1_mock, vcf_record_2_mock]}}
        assert actual == expected

    @patch.object(VCFFactory, VCFFactory.create_Pandora_VCF_from_VariantRecord_and_Sample.__name__,
                  side_effect=chrom_1_raises_NullVCFError_others_are_fine)
    def test___constructor___two_records_in_one_sample_and_two_genes___first_is_null_and_is_not_added(self, from_VariantRecord_and_Sample_Mock,
                                                                     pysam_variant_record_mock_that_maps_to_chrom_1_and_one_sample,
                                                                     pysam_variant_record_mock_that_maps_to_chrom_2_and_one_sample):
        vcf_file = VCFFile([pysam_variant_record_mock_that_maps_to_chrom_1_and_one_sample,
                            pysam_variant_record_mock_that_maps_to_chrom_2_and_one_sample], VCFFactory.create_Pandora_VCF_from_VariantRecord_and_Sample)
        actual = vcf_file.sample_to_gene_to_VCFs

        expected = {"sample_1": {"chrom_2": [vcf_record_2_mock]}}
        assert actual == expected


    @patch.object(VCFFactory, VCFFactory.create_Pandora_VCF_from_VariantRecord_and_Sample.__name__,
                  side_effect=chrom_1_raises_NullVCFError_others_are_fine)
    def test___constructor___two_records_in_one_sample_and_two_genes___second_is_null_and_is_not_added(self, from_VariantRecord_and_Sample_Mock,
                                                                     pysam_variant_record_mock_that_maps_to_chrom_1_and_one_sample,
                                                                     pysam_variant_record_mock_that_maps_to_chrom_2_and_one_sample):
        vcf_file = VCFFile([pysam_variant_record_mock_that_maps_to_chrom_2_and_one_sample,
                            pysam_variant_record_mock_that_maps_to_chrom_1_and_one_sample,], VCFFactory.create_Pandora_VCF_from_VariantRecord_and_Sample)
        actual = vcf_file.sample_to_gene_to_VCFs

        expected = {"sample_1": {"chrom_2": [vcf_record_2_mock]}}
        assert actual == expected

    @patch.object(VCFFactory, VCFFactory.create_Pandora_VCF_from_VariantRecord_and_Sample.__name__,
                  side_effect=BuggedVCFError())
    @pytest.mark.xfail(strict=True)
    def test___constructor___two_records_in_one_sample_and_two_genes___first_is_bugged___expects_death(self, from_VariantRecord_and_Sample_Mock,
                                                                     pysam_variant_record_mock_that_maps_to_chrom_1_and_one_sample,
                                                                     pysam_variant_record_mock_that_maps_to_chrom_2_and_one_sample):
        vcf_file = VCFFile([pysam_variant_record_mock_that_maps_to_chrom_1_and_one_sample,
                            pysam_variant_record_mock_that_maps_to_chrom_2_and_one_sample], VCFFactory.create_Pandora_VCF_from_VariantRecord_and_Sample)


    def test___constructor___several_records_in_several_samples_and_several_genes(self):
        pysam_variant_record_mock_that_maps_to_chrom_1_and_sample_1 = Mock()
        pysam_variant_record_mock_that_maps_to_chrom_1_and_sample_1.samples = ["sample_1"]
        pysam_variant_record_mock_that_maps_to_chrom_1_and_sample_1.chrom = "chrom_1"

        pysam_variant_record_mock_that_maps_to_chrom_1_and_sample_1_2_3 = Mock()
        pysam_variant_record_mock_that_maps_to_chrom_1_and_sample_1_2_3.samples = ["sample_1", "sample_2", "sample_3"]
        pysam_variant_record_mock_that_maps_to_chrom_1_and_sample_1_2_3.chrom = "chrom_1"

        pysam_variant_record_mock_that_maps_to_chrom_2_and_sample_1_2 = Mock()
        pysam_variant_record_mock_that_maps_to_chrom_2_and_sample_1_2.samples = ["sample_1", "sample_2"]
        pysam_variant_record_mock_that_maps_to_chrom_2_and_sample_1_2.chrom = "chrom_2"

        pysam_variant_record_mock_that_maps_to_chrom_2_and_sample_2 = Mock()
        pysam_variant_record_mock_that_maps_to_chrom_2_and_sample_2.samples = ["sample_2"]
        pysam_variant_record_mock_that_maps_to_chrom_2_and_sample_2.chrom = "chrom_2"

        another_pysam_variant_record_mock_that_maps_to_chrom_2_and_sample_2 = Mock()
        another_pysam_variant_record_mock_that_maps_to_chrom_2_and_sample_2.samples = ["sample_2"]
        another_pysam_variant_record_mock_that_maps_to_chrom_2_and_sample_2.chrom = "chrom_2"

        vcf_record_1_mock = Mock(name="vcf_record_1_mock")
        vcf_record_2_mock = Mock(name="vcf_record_2_mock")
        vcf_record_3_mock = Mock(name="vcf_record_3_mock")
        vcf_record_4_mock = Mock(name="vcf_record_4_mock")
        vcf_record_5_mock = Mock(name="vcf_record_5_mock")
        vcf_record_6_mock = Mock(name="vcf_record_6_mock")
        vcf_record_7_mock = Mock(name="vcf_record_7_mock")
        vcf_record_8_mock = Mock(name="vcf_record_8_mock")



        with patch.object(VCFFactory, VCFFactory.create_Pandora_VCF_from_VariantRecord_and_Sample.__name__,
                          side_effect=[vcf_record_1_mock, vcf_record_2_mock, vcf_record_3_mock,
                               vcf_record_4_mock, vcf_record_5_mock, vcf_record_6_mock,
                               vcf_record_7_mock, vcf_record_8_mock]):
            vcf_file = VCFFile([pysam_variant_record_mock_that_maps_to_chrom_1_and_sample_1,
                                pysam_variant_record_mock_that_maps_to_chrom_1_and_sample_1_2_3,
                                pysam_variant_record_mock_that_maps_to_chrom_2_and_sample_1_2,
                                pysam_variant_record_mock_that_maps_to_chrom_2_and_sample_2,
                                another_pysam_variant_record_mock_that_maps_to_chrom_2_and_sample_2], VCFFactory.create_Pandora_VCF_from_VariantRecord_and_Sample)
            actual = vcf_file.sample_to_gene_to_VCFs

            expected = {"sample_1": {"chrom_1": [vcf_record_1_mock, vcf_record_2_mock],
                                     "chrom_2": [vcf_record_5_mock]},
                        "sample_2": {"chrom_1": [vcf_record_3_mock],
                                     "chrom_2": [vcf_record_6_mock, vcf_record_7_mock, vcf_record_8_mock]},
                        "sample_3": {"chrom_1": [vcf_record_4_mock]}}
            assert actual == expected

    @patch.object(VCFFile, "sample_to_gene_to_VCFs", new_callable=PropertyMock, return_value = sample_with_some_genes)
    def test___get_VCF_records_given_sample_and_gene___sample_1_gene_1(self, *mocks):
        vcf_file = VCFFile([], VCFFactory.create_Pandora_VCF_from_VariantRecord_and_Sample)
        actual = vcf_file.get_VCF_records_given_sample_and_gene("sample_1", "gene_1")
        expected = [1,2,3,4]
        assert actual == expected

    @patch.object(VCFFile, "sample_to_gene_to_VCFs", new_callable=PropertyMock, return_value=sample_with_some_genes)
    def test___get_VCF_records_given_sample_and_gene___sample_1_gene_2(self, *mocks):
        vcf_file = VCFFile([], VCFFactory.create_Pandora_VCF_from_VariantRecord_and_Sample)
        actual = vcf_file.get_VCF_records_given_sample_and_gene("sample_1", "gene_2")
        expected = [5,6]
        assert actual == expected

    @patch.object(VCFFile, "sample_to_gene_to_VCFs", new_callable=PropertyMock, return_value = sample_with_some_genes)
    def test___get_VCF_records_given_sample_and_gene___sample_2_gene_1(self, *mocks):
        vcf_file = VCFFile([], VCFFactory.create_Pandora_VCF_from_VariantRecord_and_Sample)
        actual = vcf_file.get_VCF_records_given_sample_and_gene("sample_2", "gene_1")
        expected = [7, 8, 9]
        assert actual == expected

    @patch.object(VCFFile, "sample_to_gene_to_VCFs", new_callable=PropertyMock, return_value=sample_with_some_genes)
    def test___get_VCF_records_given_sample_and_gene___sample_2_gene_2(self, *mocks):
        vcf_file = VCFFile([], VCFFactory.create_Pandora_VCF_from_VariantRecord_and_Sample)
        actual = vcf_file.get_VCF_records_given_sample_and_gene("sample_2", "gene_2")
        expected = [10]
        assert actual == expected
