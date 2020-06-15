from pathlib import Path
import sys
sys.path.append(str(Path().absolute()))
from typing import List
import copy
if __name__=="__main__":
    from fix_vcf_common import FixVCF
else:
    from pipeline.scripts.fix_vcf_common import FixVCF


class FixPandoraVCF(FixVCF):
    def correct_sample_names(self, line: str, suffix: str) -> str:
        words = line.split("\t")
        corrected_words = [word.replace(suffix, "") for word in words]
        corrected_header = "\t".join(corrected_words)
        return corrected_header

    def get_gt_conf_percentiles(self, record: str) -> List[float]:
        record_split = record.split("\t")
        all_gt_conf_percentiles = []
        for index, word in enumerate(record_split):
            is_sample_info_field = index >= 9
            if is_sample_info_field:
                sample_info = word
                sample_info_split = sample_info.split(":")
                gt_conf_percentile = float(sample_info_split[-1])
                all_gt_conf_percentiles.append(gt_conf_percentile)
        return all_gt_conf_percentiles

    def get_gt_confs(self, record: str) -> List[float]:
        return self.get_gt_conf_percentiles(record)

    def set_gt_confs(self, record: str, gt_confs: List[float]) -> str:
        record_split = record.split("\t")
        record_split_corrected = []
        for index, word in enumerate(record_split):
            is_sample_info_field = index >= 9
            if is_sample_info_field:
                # correction
                sample_info = word
                sample_info_split = sample_info.split(":")
                sample_info_split_corrected = copy.deepcopy(sample_info_split)

                # assign gt_conf_percentile to gt_conf
                sample_info_split_corrected[-2] = str(gt_confs.pop(0))

                word = ":".join(sample_info_split_corrected)
            record_split_corrected.append(word)
        record_corrected = "\t".join(record_split_corrected)
        return record_corrected


    def correct_records(self, records: List[str]) -> List[str]:
        all_gt_conf_percentiles = self.get_all_gt_confs(records)
        corrected_records = self.correct_gt_confs(records, all_gt_conf_percentiles)
        return corrected_records


    def process_pandora_vcf(self, pandora_original_vcf, pandora_vcf_corrected, technology, coverage, subsampling):
        suffix = f".{coverage}.{subsampling}.{technology}"
        with open(pandora_original_vcf) as pandora_original_vcf_filehandler,\
             open(pandora_vcf_corrected, "w") as pandora_vcf_corrected_filehandler:
            headers, records = self.get_header_and_record_lines(pandora_original_vcf_filehandler)
            corrected_headers = self.correct_headers(headers, suffix)
            corrected_records = self.correct_records(records)
            print("\n".join(corrected_headers), file=pandora_vcf_corrected_filehandler)
            print("\n".join(corrected_records), file=pandora_vcf_corrected_filehandler)


if __name__=="__main__":
    # setup
    pandora_original_vcf = snakemake.input.pandora_original_vcf
    pandora_vcf_corrected = snakemake.output.pandora_vcf_corrected
    technology = snakemake.wildcards.technology
    coverage = snakemake.wildcards.coverage
    subsampling = snakemake.wildcards.subsampling
    fixer = FixPandoraVCF()
    fixer.process_pandora_vcf(pandora_original_vcf, pandora_vcf_corrected, technology, coverage, subsampling)
