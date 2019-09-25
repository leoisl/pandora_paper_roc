from unittest.mock import patch, Mock
from evaluate.vcf_file import VCFFile
from evaluate.vcf import VCF
from typing import List
from collections import defaultdict


class VCFMock(Mock):
    def __eq__(self, other):
        return self.variant == other.variant and self.sample == other.sample


def build_test_input_and_output(
    nb_of_samples: int, nb_of_records_in_each_gene: List[int]
):
    sample_to_gene_to_VCFs = defaultdict(lambda: defaultdict(list))
    list_of_variant_records = []
    samples = [f"samples_{sample_number}" for sample_number in range(nb_of_samples)]
    for gene_number, nb_of_records in enumerate(nb_of_records_in_each_gene):
        gene = f"gene_{gene_number}"
        for _ in range(nb_of_records):
            variant_record_mock = Mock(samples=samples, chrom=gene)
            list_of_variant_records.append(variant_record_mock)
            for sample in samples:
                sample_to_gene_to_VCFs[sample][gene].append(
                    VCFMock(variant=variant_record_mock, sample=sample)
                )

    return list_of_variant_records, sample_to_gene_to_VCFs


class Test_VCFFile:
    @patch("pysam.VariantFile")
    def test_constructor_noRecordsInVCFReturnsNothing(self, variant_file_mock):
        list_of_variant_records, expected = build_test_input_and_output(
            nb_of_samples=0, nb_of_records_in_each_gene=[]
        )
        variant_file_mock.return_value.__enter__.return_value.__iter__.return_value = (
            list_of_variant_records
        )

        vcf_file = VCFFile(None)
        actual = vcf_file.sample_to_gene_to_VCFs

        assert actual == expected

    @patch("pysam.VariantFile")
    @patch.object(
        VCF,
        VCF.from_VariantRecord_and_Sample.__name__,
        lambda variant, sample: VCFMock(variant=variant, sample=sample),
    )
    def test_constructor_oneRecordInVCFReturnsIndexedRecord(self, variant_file_mock):
        list_of_variant_records, expected = build_test_input_and_output(
            nb_of_samples=1, nb_of_records_in_each_gene=[1]
        )
        variant_file_mock.return_value.__enter__.return_value.__iter__.return_value = (
            list_of_variant_records
        )

        vcf_file = VCFFile(None)
        actual = vcf_file.sample_to_gene_to_VCFs

        assert actual == expected

    @patch("pysam.VariantFile")
    @patch.object(
        VCF,
        VCF.from_VariantRecord_and_Sample.__name__,
        lambda variant, sample: VCFMock(variant=variant, sample=sample),
    )
    def test_constructor_oneRecordInTwoSamplesAndOneGeneVCFReturnsIndexedRecord(
        self, variant_file_mock
    ):
        list_of_variant_records, expected = build_test_input_and_output(
            nb_of_samples=2, nb_of_records_in_each_gene=[1]
        )
        variant_file_mock.return_value.__enter__.return_value.__iter__.return_value = (
            list_of_variant_records
        )

        vcf_file = VCFFile(None)
        actual = vcf_file.sample_to_gene_to_VCFs

        assert actual == expected

    @patch("pysam.VariantFile")
    @patch.object(
        VCF,
        VCF.from_VariantRecord_and_Sample.__name__,
        lambda variant, sample: VCFMock(variant=variant, sample=sample),
    )
    def test_constructor_twoRecordsInOneSampleAndTwoGenesVCFReturnsIndexedRecords(
        self, variant_file_mock
    ):
        list_of_variant_records, expected = build_test_input_and_output(
            nb_of_samples=1, nb_of_records_in_each_gene=[1, 1]
        )
        variant_file_mock.return_value.__enter__.return_value.__iter__.return_value = (
            list_of_variant_records
        )

        vcf_file = VCFFile(None)
        actual = vcf_file.sample_to_gene_to_VCFs

        assert actual == expected

    @patch("pysam.VariantFile")
    @patch.object(
        VCF,
        VCF.from_VariantRecord_and_Sample.__name__,
        lambda variant, sample: VCFMock(variant=variant, sample=sample),
    )
    def test_constructor_twoRecordsInOneSampleAndOneGeneVCFReturnsIndexedRecords(
        self, variant_file_mock
    ):
        list_of_variant_records, expected = build_test_input_and_output(
            nb_of_samples=1, nb_of_records_in_each_gene=[2]
        )
        variant_file_mock.return_value.__enter__.return_value.__iter__.return_value = (
            list_of_variant_records
        )

        vcf_file = VCFFile(None)
        actual = vcf_file.sample_to_gene_to_VCFs

        assert actual == expected

    @patch("pysam.VariantFile")
    @patch.object(
        VCF,
        VCF.from_VariantRecord_and_Sample.__name__,
        lambda variant, sample: VCFMock(variant=variant, sample=sample),
    )
    def test_constructor_severalRecordsInSeveralSamplesAndSeveralGenesVCFReturnsIndexedRecords(
        self, variant_file_mock
    ):
        list_of_variant_records, expected = build_test_input_and_output(
            nb_of_samples=10, nb_of_records_in_each_gene=[5, 3, 4, 10, 15, 23, 41, 2]
        )
        variant_file_mock.return_value.__enter__.return_value.__iter__.return_value = (
            list_of_variant_records
        )

        vcf_file = VCFFile(None)
        actual = vcf_file.sample_to_gene_to_VCFs

        assert actual == expected
