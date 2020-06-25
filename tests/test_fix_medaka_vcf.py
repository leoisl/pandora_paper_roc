from pipeline.scripts.fix_medaka_vcf import FixMedakaVCF
import filecmp

class TestFixMedakaVCF:
    def test_big_bang_fix_medaka_vcf(self):
        fixer = FixMedakaVCF()
        fixer.process_vcf("tests/test_cases/sample_medaka_to_be_fixed.vcf",
                           "tests/test_cases/sample_medaka_to_be_fixed.corrected.vcf",
                            "sample_name")
        files_are_equal = \
            filecmp.cmp("tests/test_cases/sample_medaka_to_be_fixed.corrected.vcf",
                        "tests/test_cases/sample_medaka_to_be_fixed.expected.vcf")

        assert files_are_equal