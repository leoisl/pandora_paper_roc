from evaluate.vcf import VCFFactory, NullVCFError
import pysam
import math
import pytest

def retrieve_medaka_entry_from_test_vcf(idx: int) -> pysam.VariantRecord:
    with pysam.VariantFile("tests/test_cases/sample_medaka.vcf") as vcf:
        for i, record in enumerate(vcf):
            if i == idx:
                return record
    raise IndexError("You asked for an index that is beyond the number in the test VCF")


class Test_MedakaVCF:
    def test___medaka_first_record(self):
        medaka_variant_record = retrieve_medaka_entry_from_test_vcf(0)
        medaka_vcf = VCFFactory.create_Medaka_VCF_from_VariantRecord_and_Sample(medaka_variant_record, "SAMPLE")

        assert medaka_vcf.genotype == 1
        assert math.isclose(medaka_vcf.genotype_confidence, 41.857, abs_tol=0.0001)
        assert medaka_vcf.svtype == "NA"
        assert medaka_vcf.coverage == 1000

    def test___medaka_second_record(self):
        medaka_variant_record = retrieve_medaka_entry_from_test_vcf(1)
        medaka_vcf = VCFFactory.create_Medaka_VCF_from_VariantRecord_and_Sample(medaka_variant_record, "SAMPLE")

        assert medaka_vcf.genotype == 1
        assert math.isclose(medaka_vcf.genotype_confidence, 40.78, abs_tol=0.0001)
        assert medaka_vcf.svtype == "NA"
        assert medaka_vcf.coverage == 1000

    def test___medaka_third_record(self):
        medaka_variant_record = retrieve_medaka_entry_from_test_vcf(2)
        medaka_vcf = VCFFactory.create_Medaka_VCF_from_VariantRecord_and_Sample(medaka_variant_record, "SAMPLE")

        assert medaka_vcf.genotype == 1
        assert math.isclose(medaka_vcf.genotype_confidence, 53.109, abs_tol=0.0001)
        assert medaka_vcf.svtype == "NA"
        assert medaka_vcf.coverage == 1000

    @pytest.mark.xfail(strict=True)
    def test___medaka_fourth_record___gt_0___fails(self):
        medaka_variant_record = retrieve_medaka_entry_from_test_vcf(3)
        medaka_vcf = VCFFactory.create_Medaka_VCF_from_VariantRecord_and_Sample(medaka_variant_record, "SAMPLE")

    @pytest.mark.xfail(strict=True)
    def test___medaka_fifth_record___gt_2___fails(self):
        medaka_variant_record = retrieve_medaka_entry_from_test_vcf(4)
        medaka_vcf = VCFFactory.create_Medaka_VCF_from_VariantRecord_and_Sample(medaka_variant_record, "SAMPLE")
