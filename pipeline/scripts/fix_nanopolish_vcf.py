from pathlib import Path
import sys
sys.path.append(str(Path().absolute()))
if __name__=="__main__":
    from fix_snippy_vcf import FixSnippyVCF
else:
    from pipeline.scripts.fix_snippy_vcf import FixSnippyVCF


# FixNanopolishVCF is just like FixSnippyVCF
class FixNanopolishVCF(FixSnippyVCF):
    pass

if __name__=="__main__":
    # setup
    nanopolish_original_vcf = snakemake.input.nanopolish_original_vcf
    nanopolish_vcf_corrected = snakemake.output.nanopolish_vcf_corrected
    sample = snakemake.wildcards.sample
    fixer = FixNanopolishVCF()
    fixer.process_vcf(nanopolish_original_vcf, nanopolish_vcf_corrected, sample)