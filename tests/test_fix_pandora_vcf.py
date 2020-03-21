from pipeline.scripts.fix_pandora_vcf import process_pandora_vcf
import filecmp

class TestFixPandoraVCF:
    def test_big_bang_fix_pandora_vcf(self):
        process_pandora_vcf("tests/test_cases/pandora_multisample_genotyped_global.test.vcf",
                            "tests/test_cases/pandora_multisample_genotyped_global.test.vcf.corrected.vcf",
                            "illumina", "100x", "random")
        files_are_equal = \
            filecmp.cmp("tests/test_cases/pandora_multisample_genotyped_global.test.vcf.corrected.vcf",
                        "tests/test_cases/pandora_multisample_genotyped_global.test.vcf.expected.vcf")

        assert files_are_equal