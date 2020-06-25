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
        assert math.isclose(medaka_vcf.genotype_confidence, 9.491, abs_tol=0.0001)
        assert medaka_vcf.svtype == "NA"
        assert medaka_vcf.coverage == 100

    def test___medaka_second_record(self):
        medaka_variant_record = retrieve_medaka_entry_from_test_vcf(1)
        medaka_vcf = VCFFactory.create_Medaka_VCF_from_VariantRecord_and_Sample(medaka_variant_record, "SAMPLE")

        assert medaka_vcf.genotype == 1
        assert math.isclose(medaka_vcf.genotype_confidence, 57.336, abs_tol=0.0001)
        assert medaka_vcf.svtype == "NA"
        assert medaka_vcf.coverage == 100

    def test___medaka_third_record(self):
        medaka_variant_record = retrieve_medaka_entry_from_test_vcf(2)
        medaka_vcf = VCFFactory.create_Medaka_VCF_from_VariantRecord_and_Sample(medaka_variant_record, "SAMPLE")

        assert medaka_vcf.genotype == 0
        assert math.isclose(medaka_vcf.genotype_confidence, 2.246, abs_tol=0.0001)
        assert medaka_vcf.svtype == "NA"
        assert medaka_vcf.coverage == 100

    def test___medaka_fourth_to_ninth_record___all_inconsistent_gts___thus_null_calls(self):
        for vcf_record_index in range(3, 9):
            with pytest.raises(NullVCFError):
                medaka_variant_record = retrieve_medaka_entry_from_test_vcf(vcf_record_index)
                VCFFactory.create_Medaka_VCF_from_VariantRecord_and_Sample(medaka_variant_record, "SAMPLE")

    def test___medaka_tenth_record(self):
        medaka_variant_record = retrieve_medaka_entry_from_test_vcf(9)
        medaka_vcf = VCFFactory.create_Medaka_VCF_from_VariantRecord_and_Sample(medaka_variant_record, "SAMPLE")

        assert medaka_vcf.genotype == 2
        assert math.isclose(medaka_vcf.genotype_confidence, 2428.6575000000003, abs_tol=0.0001)
        assert medaka_vcf.svtype == "NA"
        assert medaka_vcf.coverage == 100