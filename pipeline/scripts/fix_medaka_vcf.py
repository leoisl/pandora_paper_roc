from pathlib import Path
import sys
sys.path.append(str(Path().absolute()))
if __name__=="__main__":
    from fix_snippy_vcf import FixSnippyVCF
else:
    from pipeline.scripts.fix_snippy_vcf import FixSnippyVCF


# FixMedakaVCF is just like FixSnippyVCF
class FixMedakaVCF(FixSnippyVCF):
    pass

if __name__=="__main__":
    # setup
    medaka_original_vcf = snakemake.input.medaka_original_vcf
    medaka_vcf_corrected = snakemake.output.medaka_vcf_corrected
    sample = snakemake.wildcards.sample
    fixer = FixMedakaVCF()
    fixer.process_vcf(medaka_original_vcf, medaka_vcf_corrected, sample)