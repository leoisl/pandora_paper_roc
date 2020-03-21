from pathlib import Path
import sys
sys.path.append(str(Path().absolute()))
import math

def correct_sample_names(line: str, suffix: str) -> str:
    words = line.split("\t")
    corrected_words = [word.replace(suffix, "") for word in words]
    corrected_header = "\t".join(corrected_words)
    return corrected_header


def turn_gt_conf_in_sample_to_log(word: str) -> str:
    word_split = word.split(":")
    gt_conf = float(word_split[-1])
    assert gt_conf >= 0.0, f"Error: gt_conf is negative: {gt_conf}"
    gt_conf += 1.0  # avoids calculating log of values between 0.0 and 1.0 (which can get exponentially small)
    log_gt_conf = math.log2(gt_conf)
    sample_info_with_gt_conf_as_log = ":".join(word_split[:-1] + [str(log_gt_conf)])
    return sample_info_with_gt_conf_as_log


def turn_gt_conf_into_log(line: str):
    line_split = line.split("\t")
    vcf_line_with_gt_conf_as_log = []
    for index, word in enumerate(line_split):
        if index >= 9:
            word = turn_gt_conf_in_sample_to_log(word)
        vcf_line_with_gt_conf_as_log.append(word)
    return "\t".join(vcf_line_with_gt_conf_as_log)


def process_pandora_vcf(pandora_original_vcf, pandora_vcf_corrected, technology, coverage, subsampling):
    suffix = f".{coverage}.{subsampling}.{technology}"
    with open(pandora_original_vcf) as pandora_original_vcf_filehandler,\
         open(pandora_vcf_corrected, "w") as pandora_vcf_corrected_filehandler:
        for line in pandora_original_vcf_filehandler:
            line = line.strip()
            is_header = line.startswith("#")
            is_header_with_sample_names = line.startswith("#CHROM")

            if is_header:
                if is_header_with_sample_names:
                    corrected_header = correct_sample_names(line, suffix)
                    print(corrected_header, file=pandora_vcf_corrected_filehandler)
                else:
                    print(line, file=pandora_vcf_corrected_filehandler)
            else:
                corrected_vcf = turn_gt_conf_into_log(line)
                print(corrected_vcf, file=pandora_vcf_corrected_filehandler)



if __name__=="__main__":
    # setup
    pandora_original_vcf = snakemake.input.pandora_original_vcf
    pandora_vcf_corrected = snakemake.output.pandora_vcf_corrected
    technology = snakemake.wildcards.technology
    coverage = snakemake.wildcards.coverage
    subsampling = snakemake.wildcards.subsampling
    process_pandora_vcf(pandora_original_vcf, pandora_vcf_corrected, technology, coverage, subsampling)
