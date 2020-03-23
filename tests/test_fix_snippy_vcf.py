from pipeline.scripts.fix_snippy_vcf import FixSnippyVCF
import filecmp

class TestFixSnippyVCF:
    def test_big_bang_fix_snippy_vcf(self):
        fixer = FixSnippyVCF()
        fixer.process_snippy_vcf("tests/test_cases/sample_snippy_to_be_fixed.vcf",
                           "tests/test_cases/sample_snippy_to_be_fixed.corrected.vcf",
                            "sample_name")
        files_are_equal = \
            filecmp.cmp("tests/test_cases/sample_snippy_to_be_fixed.corrected.vcf",
                        "tests/test_cases/sample_snippy_to_be_fixed.expected.vcf")

        assert files_are_equal