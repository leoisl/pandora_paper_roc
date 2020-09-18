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
from evaluate.filtered_vcf_file import FilteredVCFFile
from evaluate.vcf_filters import VCF_Filters
from evaluate.vcf import VCFFactory
import pysam
import subprocess
def run_command(command):
    subprocess.check_call(command, shell=True)


# setup
singlesample_vcf_files_gt_conf_percentile_filtered = snakemake.input.singlesample_vcf_files_gt_conf_percentile_filtered
vcf_ref = snakemake.input.vcf_ref
empty_depth_file = snakemake.input.empty_depth_file
filtered_vcf_filepaths = snakemake.output.filtered_vcf_files
mutated_vcf_refs = snakemake.output.mutated_vcf_refs
gt_conf_percentiles = snakemake.params.gt_conf_percentiles
coverage_threshold = snakemake.wildcards.coverage_threshold
strand_bias_threshold = snakemake.wildcards.strand_bias_threshold
gaps_threshold = snakemake.wildcards.gaps_threshold
tool = snakemake.wildcards.tool

filters = VCF_Filters.get_all_VCF_Filters(
    coverage_threshold=coverage_threshold,
    strand_bias_threshold=strand_bias_threshold,
    gaps_threshold=gaps_threshold,
)

if tool.startswith("pandora"):
    VCF_creator_method = VCFFactory.create_Pandora_VCF_from_VariantRecord_and_Sample
elif tool.startswith("snippy"):
    VCF_creator_method = VCFFactory.create_Snippy_VCF_from_VariantRecord_and_Sample
elif tool.startswith("samtools"):
    VCF_creator_method = VCFFactory.create_Samtools_VCF_from_VariantRecord_and_Sample
elif tool.startswith("medaka"):
    VCF_creator_method = VCFFactory.create_Medaka_VCF_from_VariantRecord_and_Sample
elif tool.startswith("nanopolish"):
    VCF_creator_method = VCFFactory.create_Nanopolish_VCF_from_VariantRecord_and_Sample
else:
    raise RuntimeError("VCFs should be from either pandora or snippy or samtools or medaka or nanopolish (should start with either these values)")

for gt_conf_percentile, vcf_filepath, filtered_vcf_filepath, mutated_vcf_ref in \
        zip(gt_conf_percentiles, singlesample_vcf_files_gt_conf_percentile_filtered, filtered_vcf_filepaths, mutated_vcf_refs):
    logging.info(f"Applying filters to {vcf_filepath}")
    with pysam.VariantFile(vcf_filepath) as pysam_variant_file:
        filtered_vcf_file = FilteredVCFFile(pysam_variant_file=pysam_variant_file, filters=filters,
                                            VCF_creator_method=VCF_creator_method)
    with open(filtered_vcf_filepath, "w") as filtered_vcf_filehandler:
        filtered_vcf_file.write(filtered_vcf_filehandler)

    logging.info("Running vcf_consensus_builder")
    run_command(f"python vcf_consensus_builder/cli.py -v {filtered_vcf_filepath} -r {vcf_ref} -d {empty_depth_file} "
                f"-o {mutated_vcf_ref} --low-coverage 0 --no-coverage 0 -V")

    logging.info("bwa index")
    run_command(f"bwa index {mutated_vcf_ref}")
