from pathlib import Path
import sys
sys.path.append(str(Path().absolute()))
if __name__=="__main__":
    from fix_snippy_vcf import FixSnippyVCF
else:
    from pipeline.scripts.fix_snippy_vcf import FixSnippyVCF


# FixSamtoolsVCF is just like FixSnippyVCF
class FixSamtoolsVCF(FixSnippyVCF):
    pass

if __name__=="__main__":
    # setup
    samtools_original_vcf = snakemake.input.samtools_original_vcf
    samtools_vcf_corrected = snakemake.output.samtools_vcf_corrected
    sample = snakemake.wildcards.sample
    fixer = FixSamtoolsVCF()
    fixer.process_vcf(samtools_original_vcf, samtools_vcf_corrected, sample)