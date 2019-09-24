from unittest.mock import patch, PropertyMock, MagicMock
import pysam
from evaluate.vcf_file import VCFFile
from pytest import fixture
from evaluate.vcf import VCF


variantMock1 = MagicMock()
variantMock2 = MagicMock()

class Test_VCFFile:
    @patch.object(VCFFile, VCFFile.__init__.__name__, return_value=None)
    @patch.object(VCFFile, VCFFile.subset_samples.__name__)
    @patch.object(VCFFile, VCFFile.fetch.__name__)
    def test_getVCFRecordsGivenSampleAndGene_noRecordsInVariantFileMatchingSampleAndGeneReturnsEmpty(self, mock_fetch, mock_subset_samples, mock_init):
        vcf_file = VCFFile()

        actual = vcf_file.get_VCF_records_given_sample_and_gene("sample", "gene")
        expected = []

        assert actual == expected
        mock_subset_samples.assert_called_once_with(["sample"])
        mock_fetch.assert_called_once_with("gene")


    @patch.object(VCFFile, VCFFile.__init__.__name__, return_value=None)
    @patch.object(VCFFile, VCFFile.subset_samples.__name__)
    @patch.object(VCFFile, VCFFile.fetch.__name__, return_value=[variantMock1])
    def test_getVCFRecordsGivenSampleAndGene_oneRecordInVariantFileMatchingSampleAndGeneReturnsRecord(self, mock_fetch, mock_subset_samples, mock_init):
        vcf_file = VCFFile()
        vcf_records = vcf_file.get_VCF_records_given_sample_and_gene("sample", "gene")

        actual = [vcf_record.__dict__ for vcf_record in vcf_records]
        expected = [VCF(variantMock1, "sample").__dict__]

        assert actual == expected
        mock_subset_samples.assert_called_once_with(["sample"])
        mock_fetch.assert_called_once_with("gene")


    @patch.object(VCFFile, VCFFile.__init__.__name__, return_value=None)
    @patch.object(VCFFile, VCFFile.subset_samples.__name__)
    @patch.object(VCFFile, VCFFile.fetch.__name__, return_value=[variantMock1, variantMock2])
    def test_getVCFRecordsGivenSampleAndGene_twoRecordsInVariantFileMatchingSampleAndGeneReturnsTwoRecords(self, mock_fetch, mock_subset_samples, mock_init):
        vcf_file = VCFFile()
        vcf_records = vcf_file.get_VCF_records_given_sample_and_gene("sample", "gene")

        actual = [vcf_record.__dict__ for vcf_record in vcf_records]
        expected = [VCF(variantMock1, "sample").__dict__, VCF(variantMock2, "sample").__dict__]

        assert actual == expected
        mock_subset_samples.assert_called_once_with(["sample"])
        mock_fetch.assert_called_once_with("gene")


    @patch.object(VCFFile, VCFFile.__init__.__name__, return_value=None)
    @patch.object(VCFFile, VCFFile.subset_samples.__name__)
    @patch.object(VCFFile, VCFFile.fetch.__name__, return_value=[variantMock2])
    def test_getVCFRecordsGivenSampleAndGene_twoRecordsInVariantFileOneMatchingSampleAndGeneReturnsOneRecord(self, mock_fetch, mock_subset_samples, mock_init):
        vcf_file = VCFFile()
        vcf_records = vcf_file.get_VCF_records_given_sample_and_gene("sample", "gene")

        actual = [vcf_record.__dict__ for vcf_record in vcf_records]
        expected = [VCF(variantMock2, "sample").__dict__]

        assert actual == expected
        mock_subset_samples.assert_called_once_with(["sample"])
        mock_fetch.assert_called_once_with("gene")
