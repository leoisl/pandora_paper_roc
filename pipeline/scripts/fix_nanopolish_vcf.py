from pathlib import Path
import sys
sys.path.append(str(Path().absolute()))
if __name__=="__main__":
    from fix_snippy_vcf import FixSnippyVCF
else:
    from pipeline.scripts.fix_snippy_vcf import FixSnippyVCF


class FixNanopolishVCF(FixSnippyVCF):
    def process_vcf(self, original_vcf, corrected_vcf, sample):
        with open(original_vcf) as original_vcf_filehandler,\
             open(corrected_vcf, "w") as corrected_vcf_filehandler:
            headers, records = self.get_header_and_record_lines(original_vcf_filehandler)
            corrected_headers = self.correct_headers(headers, sample)
            print("\n".join(corrected_headers), file=corrected_vcf_filehandler)

            if len(records) > 0:
                corrected_records = self.correct_records(records)

                # remove repeated records
                corrected_records = list(set(corrected_records))

                # sort records by chrom and pos
                corrected_records = sorted(corrected_records,
                                           key=lambda record: (record.split()[0], int(record.split()[1])))

                print("\n".join(corrected_records), file=corrected_vcf_filehandler)


if __name__=="__main__":
    # setup
    nanopolish_original_vcf = snakemake.input.nanopolish_original_vcf
    nanopolish_vcf_corrected = snakemake.output.nanopolish_vcf_corrected
    sample = snakemake.wildcards.sample
    fixer = FixNanopolishVCF()
    fixer.process_vcf(nanopolish_original_vcf, nanopolish_vcf_corrected, sample)