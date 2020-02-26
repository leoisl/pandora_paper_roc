from pathlib import Path
import sys
sys.path.append(str(Path().absolute()))
import logging
log_level = "INFO"
logging.basicConfig(
    filename=str(snakemake.log),
    filemode="w",
    level=log_level,
    format="[%(asctime)s]:%(levelname)s: %(message)s",
    datefmt="%d/%m/%Y %I:%M:%S %p",
)

def correct(word, suffix):
    return word.replace(suffix, "")

# setup
pandora_original_vcf = snakemake.input.pandora_original_vcf
pandora_vcf_with_sample_names_corrected = snakemake.output.pandora_vcf_with_sample_names_corrected
technology = snakemake.wildcards.technology
coverage = snakemake.wildcards.coverage
subsampling = snakemake.wildcards.subsampling


# processing
suffix = f".{coverage}.{subsampling}.{technology}"
with open(pandora_original_vcf) as pandora_original_vcf_filehandler,\
     open(pandora_vcf_with_sample_names_corrected, "w") as pandora_vcf_with_sample_names_corrected_filehandler:
    for line in pandora_original_vcf_filehandler:
        if line.startswith("#CHROM"):
            words = line.strip().split()
            corrected_words = [word.replace(suffix, "") for word in words]
            print("\t".join(corrected_words), file=pandora_vcf_with_sample_names_corrected_filehandler)
        else:
            pandora_vcf_with_sample_names_corrected_filehandler.write(line)


