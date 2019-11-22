from unittest.mock import patch, Mock
from evaluate.vcf_file import VCFFile
from evaluate.vcf import VCF, NullVCFError, BuggedVCFError
import pytest

# TODO : put this in a Fixture
pysam_variant_record_mock_that_maps_to_chrom_1_and_one_sample = Mock()
pysam_variant_record_mock_that_maps_to_chrom_1_and_one_sample.samples = ["sample_1"]
pysam_variant_record_mock_that_maps_to_chrom_1_and_one_sample.chrom = "chrom_1"

pysam_variant_record_mock_that_maps_to_chrom_2_and_one_sample = Mock()
pysam_variant_record_mock_that_maps_to_chrom_2_and_one_sample.samples = ["sample_1"]
pysam_variant_record_mock_that_maps_to_chrom_2_and_one_sample.chrom = "chrom_2"

pysam_variant_record_mock_that_maps_to_chrom_1_and_two_samples = Mock()
pysam_variant_record_mock_that_maps_to_chrom_1_and_two_samples.samples = ["sample_1", "sample_2"]
pysam_variant_record_mock_that_maps_to_chrom_1_and_two_samples.chrom = "chrom_1"


vcf_record_1_mock = Mock()
vcf_record_2_mock = Mock()
vcf_record_3_mock = Mock()


def chrom_1_raises_NullVCFError_others_are_fine(pysam_variant_record, sample):
    if pysam_variant_record.chrom == "chrom_1":
        raise NullVCFError()
    else:
        return vcf_record_2_mock


class Test_VCFFile:
    def test___constructor___no_records_in_VCF_returns_nothing(self):
        vcf_file = VCFFile([])
        actual = vcf_file.sample_to_gene_to_VCFs

        expected = {}
        assert actual == expected

    @patch.object(VCF, VCF.from_VariantRecord_and_Sample.__name__, return_value=vcf_record_1_mock)
    def test___constructor___one_record_in_one_sample_and_one_gene(self, *mocks):
        vcf_file = VCFFile([pysam_variant_record_mock_that_maps_to_chrom_1_and_one_sample])
        actual = vcf_file.sample_to_gene_to_VCFs

        expected = {"sample_1": {"chrom_1": [vcf_record_1_mock]}}
        assert actual == expected

    @patch.object(VCF, VCF.from_VariantRecord_and_Sample.__name__,
                  return_value=vcf_record_1_mock)
    def test___constructor___one_record_in_two_samples_and_one_gene(self, *mocks):
        vcf_file = VCFFile([pysam_variant_record_mock_that_maps_to_chrom_1_and_two_samples])
        actual = vcf_file.sample_to_gene_to_VCFs

        expected = {"sample_1": {"chrom_1": [vcf_record_1_mock]},
                    "sample_2": {"chrom_1": [vcf_record_1_mock]}}
        assert actual == expected

    @patch.object(VCF, VCF.from_VariantRecord_and_Sample.__name__,
                  side_effect=[vcf_record_1_mock, vcf_record_2_mock])
    def test___constructor___two_records_in_one_sample_and_two_genes(self, *mocks):
        vcf_file = VCFFile([pysam_variant_record_mock_that_maps_to_chrom_1_and_one_sample,
                            pysam_variant_record_mock_that_maps_to_chrom_2_and_one_sample])
        actual = vcf_file.sample_to_gene_to_VCFs

        expected = {"sample_1": {"chrom_1": [vcf_record_1_mock], "chrom_2": [vcf_record_2_mock]}}
        assert actual == expected

    @patch.object(VCF, VCF.from_VariantRecord_and_Sample.__name__,
                  side_effect=[vcf_record_1_mock, vcf_record_2_mock])
    def test___constructor___two_records_in_one_sample_and_one_gene(self, *mocks):
        vcf_file = VCFFile([pysam_variant_record_mock_that_maps_to_chrom_1_and_one_sample,
                            pysam_variant_record_mock_that_maps_to_chrom_1_and_one_sample])
        actual = vcf_file.sample_to_gene_to_VCFs

        expected = {"sample_1": {"chrom_1": [vcf_record_1_mock, vcf_record_2_mock]}}
        assert actual == expected

    @patch.object(VCF, VCF.from_VariantRecord_and_Sample.__name__,
                  side_effect=chrom_1_raises_NullVCFError_others_are_fine)
    def test___constructor___two_records_in_one_sample_and_two_genes___first_is_null_and_is_not_added(self, *mocks):
        vcf_file = VCFFile([pysam_variant_record_mock_that_maps_to_chrom_1_and_one_sample,
                            pysam_variant_record_mock_that_maps_to_chrom_2_and_one_sample])
        actual = vcf_file.sample_to_gene_to_VCFs

        expected = {"sample_1": {"chrom_2": [vcf_record_2_mock]}}
        assert actual == expected


    @patch.object(VCF, VCF.from_VariantRecord_and_Sample.__name__,
                  side_effect=chrom_1_raises_NullVCFError_others_are_fine)
    def test___constructor___two_records_in_one_sample_and_two_genes___second_is_null_and_is_not_added(self, *mocks):
        vcf_file = VCFFile([pysam_variant_record_mock_that_maps_to_chrom_2_and_one_sample,
                            pysam_variant_record_mock_that_maps_to_chrom_1_and_one_sample,])
        actual = vcf_file.sample_to_gene_to_VCFs

        expected = {"sample_1": {"chrom_2": [vcf_record_2_mock]}}
        assert actual == expected

    @patch.object(VCF, VCF.from_VariantRecord_and_Sample.__name__,
                  side_effect=BuggedVCFError())
    @pytest.mark.xfail(strict=True)
    def test___constructor___two_records_in_one_sample_and_two_genes___first_is_bugged___expects_death(self, *mocks):
        vcf_file = VCFFile([pysam_variant_record_mock_that_maps_to_chrom_1_and_one_sample,
                            pysam_variant_record_mock_that_maps_to_chrom_2_and_one_sample])




    def test___constructor___several_records_in_several_samples_and_several_genes(self, *mocks):
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



        with patch.object(VCF, VCF.from_VariantRecord_and_Sample.__name__,
                          side_effect=[vcf_record_1_mock, vcf_record_2_mock, vcf_record_3_mock,
                               vcf_record_4_mock, vcf_record_5_mock, vcf_record_6_mock,
                               vcf_record_7_mock, vcf_record_8_mock]):
            vcf_file = VCFFile([pysam_variant_record_mock_that_maps_to_chrom_1_and_sample_1,
                                pysam_variant_record_mock_that_maps_to_chrom_1_and_sample_1_2_3,
                                pysam_variant_record_mock_that_maps_to_chrom_2_and_sample_1_2,
                                pysam_variant_record_mock_that_maps_to_chrom_2_and_sample_2,
                                another_pysam_variant_record_mock_that_maps_to_chrom_2_and_sample_2])
            actual = vcf_file.sample_to_gene_to_VCFs

            expected = {"sample_1": {"chrom_1": [vcf_record_1_mock, vcf_record_2_mock],
                                     "chrom_2": [vcf_record_5_mock]},
                        "sample_2": {"chrom_1": [vcf_record_3_mock],
                                     "chrom_2": [vcf_record_6_mock, vcf_record_7_mock, vcf_record_8_mock]},
                        "sample_3": {"chrom_1": [vcf_record_4_mock]}}
            assert actual == expected
