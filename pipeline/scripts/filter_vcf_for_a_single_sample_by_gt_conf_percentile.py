import subprocess
def run_command(command):
    subprocess.check_call(command, shell=True)

for gt_conf_percentile, output_file in zip(snakemake.params.gt_conf_percentiles,
                                           snakemake.output.singlesample_vcf_files_gt_conf_percentile_filtered):
    run_command(f"bash {snakemake.params.filter_script} "
                f"{snakemake.input.gzipped_singlesample_vcf_file} "
                f"{gt_conf_percentile} "
                f"{output_file}")
