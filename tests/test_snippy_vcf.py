from evaluate.vcf import VCFFactory
import pysam
import math
import pytest

def retrieve_snippy_entry_from_test_vcf(idx: int) -> pysam.VariantRecord:
    with pysam.VariantFile("tests/test_cases/sample_snippy_H131800734_NZ_CP007265.1.vcf") as vcf:
        for i, record in enumerate(vcf):
            if i == idx:
                return record
    raise IndexError("You asked for an index that is beyond the number in the test VCF")


class Test_SnippyVCF:
    def test___snippy_first_real_record(self):
        snippy_variant_record = retrieve_snippy_entry_from_test_vcf(0)
        snippy_vcf = VCFFactory.create_Snippy_VCF_from_VariantRecord_and_Sample(snippy_variant_record, "H131800734")

        assert snippy_vcf.genotype == 1
        assert math.isclose(snippy_vcf.genotype_confidence, 2546.69, abs_tol=0.0001)
        assert snippy_vcf.svtype == "snp"
        assert snippy_vcf.coverage == 1000

    def test___snippy_second_real_record(self):
        snippy_variant_record = retrieve_snippy_entry_from_test_vcf(1)
        snippy_vcf = VCFFactory.create_Snippy_VCF_from_VariantRecord_and_Sample(snippy_variant_record, "H131800734")

        assert snippy_vcf.genotype == 1
        assert math.isclose(snippy_vcf.genotype_confidence, 1843.76, abs_tol=0.0001)
        assert snippy_vcf.svtype == "snp"
        assert snippy_vcf.coverage == 1000

    def test___snippy_third_real_record(self):
        snippy_variant_record = retrieve_snippy_entry_from_test_vcf(2)
        snippy_vcf = VCFFactory.create_Snippy_VCF_from_VariantRecord_and_Sample(snippy_variant_record, "H131800734")

        assert snippy_vcf.genotype == 1
        assert math.isclose(snippy_vcf.genotype_confidence, 1957.94, abs_tol=0.0001)
        assert snippy_vcf.svtype == "complex"
        assert snippy_vcf.coverage == 1000

    @pytest.mark.xfail(strict=True)
    def test___snippy_fails_with_genotype_1_0(self):
        snippy_variant_record = retrieve_snippy_entry_from_test_vcf(3)
        snippy_vcf = VCFFactory.create_Snippy_VCF_from_VariantRecord_and_Sample(snippy_variant_record, "H131800734")

    @pytest.mark.xfail(strict=True)
    def test___snippy_fails_with_genotype_0_1(self):
        snippy_variant_record = retrieve_snippy_entry_from_test_vcf(4)
        snippy_vcf = VCFFactory.create_Snippy_VCF_from_VariantRecord_and_Sample(snippy_variant_record, "H131800734")

    @pytest.mark.xfail(strict=True)
    def test___snippy_fails_with_genotype_0_0(self):
        snippy_variant_record = retrieve_snippy_entry_from_test_vcf(5)
        snippy_vcf = VCFFactory.create_Snippy_VCF_from_VariantRecord_and_Sample(snippy_variant_record, "H131800734")