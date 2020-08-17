from pipeline.scripts.fix_samtools_vcf import FixSamtoolsVCF
import filecmp

class TestFixSamtoolsVCF:
    def test_big_bang_fix_samtools_vcf(self):
        fixer = FixSamtoolsVCF()
        fixer.process_vcf("tests/test_cases/sample_samtools_to_be_fixed.vcf",
                           "tests/test_cases/sample_samtools_to_be_fixed.corrected.vcf",
                            "sample_name")
        files_are_equal = \
            filecmp.cmp("tests/test_cases/sample_samtools_to_be_fixed.corrected.vcf",
                        "tests/test_cases/sample_samtools_to_be_fixed.expected.vcf")

        assert files_are_equal