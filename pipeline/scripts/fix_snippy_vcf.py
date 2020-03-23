from pathlib import Path
import sys
sys.path.append(str(Path().absolute()))
from typing import List
import copy
if __name__=="__main__":
    from fix_vcf_common import FixVCF
else:
    from pipeline.scripts.fix_vcf_common import FixVCF


class FixSnippyVCF(FixVCF):
    def correct_sample_names(self, line: str, sample: str) -> str:
        words = line.split("\t")
        header_has_ten_fields = len(words) == 10
        assert header_has_ten_fields, f"Snippy {line} should have 10 fields (only one sample)."
        corrected_words = words
        corrected_words[-1] = sample
        return "\t".join(corrected_words)


    def get_gt_confs(self, record: str) -> List[float]:
        record_split = record.split("\t")
        qual_field_index = 5
        all_gt_confs = [ float(record_split[qual_field_index]) ]
        return all_gt_confs


    def set_gt_confs(self, record: str, gt_confs: List[float]) -> str:
        record_split = record.split("\t")
        record_split_corrected = copy.deepcopy(record_split)
        qual_field_index = 5
        record_split_corrected[qual_field_index] = str(gt_confs.pop(0))
        record_corrected = "\t".join(record_split_corrected)
        return record_corrected


    def process_snippy_vcf(self, snippy_original_vcf, snippy_vcf_corrected, sample):
        with open(snippy_original_vcf) as snippy_original_vcf_filehandler,\
             open(snippy_vcf_corrected, "w") as snippy_vcf_corrected_filehandler:
            headers, records = self.get_header_and_record_lines(snippy_original_vcf_filehandler)
            corrected_headers = self.correct_headers(headers, sample)
            corrected_records = self.correct_records(records)
            print("\n".join(corrected_headers), file=snippy_vcf_corrected_filehandler)
            print("\n".join(corrected_records), file=snippy_vcf_corrected_filehandler)


if __name__=="__main__":
    # setup
    snippy_original_vcf = snakemake.input.snippy_original_vcf
    snippy_vcf_corrected = snakemake.output.snippy_vcf_corrected
    sample = snakemake.wildcards.sample
    fixer = FixSnippyVCF()
    fixer.process_snippy_vcf(snippy_original_vcf, snippy_vcf_corrected, sample)