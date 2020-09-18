from pipeline.scripts.fix_nanopolish_vcf import FixNanopolishVCF
import filecmp

class TestFixNanopolishVCF:
    def test_big_bang_fix_nanopolish_vcf(self):
        fixer = FixNanopolishVCF()
        fixer.process_vcf("tests/test_cases/sample_nanopolish_to_be_fixed.vcf",
                           "tests/test_cases/sample_nanopolish_to_be_fixed.corrected.vcf",
                            "sample_name")
        files_are_equal = \
            filecmp.cmp("tests/test_cases/sample_nanopolish_to_be_fixed.corrected.vcf",
                        "tests/test_cases/sample_nanopolish_to_be_fixed.expected.vcf")

        assert files_are_equal
