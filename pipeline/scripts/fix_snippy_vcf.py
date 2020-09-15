from pathlib import Path
import sys
sys.path.append(str(Path().absolute()))
from typing import Deque
import copy
from collections import deque
if __name__=="__main__":
    from fix_vcf_common import FixVCF
else:
    from pipeline.scripts.fix_vcf_common import FixVCF


class FixSnippyVCF(FixVCF):
    def correct_sample_names(self, line: str, sample: str) -> str:
        words = line.split("\t")
        header_has_ten_fields = len(words) == 10
        assert header_has_ten_fields, f"{line} should have 10 fields (only one sample)."
        corrected_words = words
        corrected_words[-1] = sample
        return "\t".join(corrected_words)


    def get_gt_confs(self, record: str) -> Deque[float]:
        record_split = record.split("\t")
        qual_field_index = 5
        all_gt_confs = deque([ float(record_split[qual_field_index]) ])
        return all_gt_confs


    def set_gt_confs(self, record: str, gt_confs: Deque[float]) -> str:
        record_split = record.split("\t")
        record_split_corrected = copy.deepcopy(record_split)
        qual_field_index = 5
        record_split_corrected[qual_field_index] = str(gt_confs.popleft())
        record_corrected = "\t".join(record_split_corrected)
        return record_corrected


    def process_vcf(self, original_vcf, corrected_vcf, sample):
        with open(original_vcf) as original_vcf_filehandler,\
             open(corrected_vcf, "w") as corrected_vcf_filehandler:
            headers, records = self.get_header_and_record_lines(original_vcf_filehandler)
            corrected_headers = self.correct_headers(headers, sample)
            print("\n".join(corrected_headers), file=corrected_vcf_filehandler)

            if len(records) > 0:
                corrected_records = self.correct_records(records)
                print("\n".join(corrected_records), file=corrected_vcf_filehandler)


if __name__=="__main__":
    # setup
    snippy_original_vcf = snakemake.input.snippy_original_vcf
    snippy_vcf_corrected = snakemake.output.snippy_vcf_corrected
    sample = snakemake.wildcards.sample
    fixer = FixSnippyVCF()
    fixer.process_vcf(snippy_original_vcf, snippy_vcf_corrected, sample)