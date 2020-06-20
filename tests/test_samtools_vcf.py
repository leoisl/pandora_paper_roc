from evaluate.vcf import VCFFactory
import pysam
import math
import pytest

def retrieve_samtools_entry_from_test_vcf(idx: int) -> pysam.VariantRecord:
    with pysam.VariantFile("tests/test_cases/sample_samtools_to_be_fixed.expected.vcf") as vcf:
        for i, record in enumerate(vcf):
            if i == idx:
                return record
    raise IndexError("You asked for an index that is beyond the number in the test VCF")


class Test_SamtoolsVCF:
    def test___samtools_first_real_record(self):
        samtools_variant_record = retrieve_samtools_entry_from_test_vcf(0)
        samtools_vcf = VCFFactory.create_Samtools_VCF_from_VariantRecord_and_Sample(samtools_variant_record, "sample_name")

        assert samtools_vcf.genotype == 1
        assert math.isclose(samtools_vcf.genotype_confidence, 100, abs_tol=0.0001)
        assert samtools_vcf.svtype == "SNP"
        assert samtools_vcf.coverage == 66

    def test___samtools_second_real_record(self):
        samtools_variant_record = retrieve_samtools_entry_from_test_vcf(1)
        samtools_vcf = VCFFactory.create_Samtools_VCF_from_VariantRecord_and_Sample(samtools_variant_record, "sample_name")

        assert samtools_vcf.genotype == 1
        assert math.isclose(samtools_vcf.genotype_confidence, 85.4, abs_tol=0.0001)
        assert samtools_vcf.svtype == "INDEL"
        assert samtools_vcf.coverage == 83

    def test___samtools_third_real_record(self):
        samtools_variant_record = retrieve_samtools_entry_from_test_vcf(2)
        samtools_vcf = VCFFactory.create_Samtools_VCF_from_VariantRecord_and_Sample(samtools_variant_record, "sample_name")

        assert samtools_vcf.genotype == 1
        assert math.isclose(samtools_vcf.genotype_confidence, 0.0, abs_tol=0.0001)
        assert samtools_vcf.svtype == "SNP"
        assert samtools_vcf.coverage == 82

    @pytest.mark.xfail(strict=True)
    def test___samtools_fails_with_genotype_0(self):
        samtools_variant_record = retrieve_samtools_entry_from_test_vcf(3)
        samtools_vcf = VCFFactory.create_Samtools_VCF_from_VariantRecord_and_Sample(samtools_variant_record, "sample_name")

    @pytest.mark.xfail(strict=True)
    def test___samtools_fails_with_genotype_2(self):
        samtools_variant_record = retrieve_samtools_entry_from_test_vcf(4)
        samtools_vcf = VCFFactory.create_Samtools_VCF_from_VariantRecord_and_Sample(samtools_variant_record, "sample_name")
