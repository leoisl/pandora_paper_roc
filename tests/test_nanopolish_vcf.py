from evaluate.vcf import VCFFactory, NullVCFError
import pysam
import math
import pytest

def retrieve_nanopolish_entry_from_test_vcf(idx: int) -> pysam.VariantRecord:
    with pysam.VariantFile("tests/test_cases/sample_nanopolish.vcf") as vcf:
        for i, record in enumerate(vcf):
            if i == idx:
                return record
    raise IndexError("You asked for an index that is beyond the number in the test VCF")


class Test_NanopolishVCF:
    def test___nanopolish_first_record(self):
        nanopolish_variant_record = retrieve_nanopolish_entry_from_test_vcf(0)
        nanopolish_vcf = VCFFactory.create_Nanopolish_VCF_from_VariantRecord_and_Sample(nanopolish_variant_record, "SAMPLE")

        assert nanopolish_vcf.genotype == 1
        assert math.isclose(nanopolish_vcf.genotype_confidence, 121.8, abs_tol=0.0001)
        assert nanopolish_vcf.svtype == "NA"
        assert nanopolish_vcf.coverage == 26

    def test___nanopolish_second_record(self):
        nanopolish_variant_record = retrieve_nanopolish_entry_from_test_vcf(1)
        nanopolish_vcf = VCFFactory.create_Nanopolish_VCF_from_VariantRecord_and_Sample(nanopolish_variant_record, "SAMPLE")

        assert nanopolish_vcf.genotype == 1
        assert math.isclose(nanopolish_vcf.genotype_confidence, 652.1, abs_tol=0.0001)
        assert nanopolish_vcf.svtype == "NA"
        assert nanopolish_vcf.coverage == 63

    def test___nanopolish_third_record(self):
        nanopolish_variant_record = retrieve_nanopolish_entry_from_test_vcf(2)
        nanopolish_vcf = VCFFactory.create_Nanopolish_VCF_from_VariantRecord_and_Sample(nanopolish_variant_record, "SAMPLE")

        assert nanopolish_vcf.genotype == 1
        assert math.isclose(nanopolish_vcf.genotype_confidence, 560.2, abs_tol=0.0001)
        assert nanopolish_vcf.svtype == "NA"
        assert nanopolish_vcf.coverage == 81

    @pytest.mark.xfail(strict=True)
    def test___nanopolish_fourth_record___gt_0___fails(self):
        nanopolish_variant_record = retrieve_nanopolish_entry_from_test_vcf(3)
        nanopolish_vcf = VCFFactory.create_Nanopolish_VCF_from_VariantRecord_and_Sample(nanopolish_variant_record, "SAMPLE")

    @pytest.mark.xfail(strict=True)
    def test___nanopolish_fifth_record___gt_2___fails(self):
        nanopolish_variant_record = retrieve_nanopolish_entry_from_test_vcf(4)
        nanopolish_vcf = VCFFactory.create_Nanopolish_VCF_from_VariantRecord_and_Sample(nanopolish_variant_record, "SAMPLE")
