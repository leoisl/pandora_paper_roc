from pathlib import Path
import sys
sys.path.append(str(Path().absolute()))
import math


def correct_sample_names(line: str, sample: str) -> str:
    words = line.split("\t")
    header_has_ten_fields = len(words) == 10
    assert header_has_ten_fields, f"Snippy {line} should have 10 fields (only one sample)."
    corrected_words = words
    corrected_words[-1] = sample
    return "\t".join(corrected_words)


def turn_gt_conf_in_sample_to_log(qual: str) -> str:
    gt_conf = float(qual)
    assert gt_conf >= 0.0, f"Error: gt_conf is negative: {gt_conf}"
    gt_conf += 1.0  # avoids calculating log of values between 0.0 and 1.0 (which can get exponentially small)
    log_gt_conf = math.log2(gt_conf)
    return str(log_gt_conf)


def turn_gt_conf_into_log(line: str):
    line_split = line.split("\t")
    vcf_line_with_gt_conf_as_log = []
    for index, word in enumerate(line_split):
        if index == 5:
            word = turn_gt_conf_in_sample_to_log(word)
        vcf_line_with_gt_conf_as_log.append(word)
    return "\t".join(vcf_line_with_gt_conf_as_log)


def process_snippy_vcf(snippy_original_vcf, snippy_vcf_corrected, sample):
    with open(snippy_original_vcf) as snippy_original_vcf_filehandler,\
         open(snippy_vcf_corrected, "w") as snippy_vcf_corrected_filehandler:
        for line in snippy_original_vcf_filehandler:
            line = line.strip()
            is_header = line.startswith("#")
            is_header_with_sample_names = line.startswith("#CHROM")

            if is_header:
                if is_header_with_sample_names:
                    corrected_header = correct_sample_names(line, sample)
                    print(corrected_header, file=snippy_vcf_corrected_filehandler)
                else:
                    print(line, file=snippy_vcf_corrected_filehandler)
            else:
                corrected_vcf = turn_gt_conf_into_log(line)
                print(corrected_vcf, file=snippy_vcf_corrected_filehandler)


if __name__=="__main__":
    # setup
    snippy_original_vcf = snakemake.input.snippy_original_vcf
    snippy_vcf_corrected = snakemake.output.snippy_vcf_corrected
    sample = snakemake.wildcards.sample
    process_snippy_vcf(snippy_original_vcf, snippy_vcf_corrected, sample)
