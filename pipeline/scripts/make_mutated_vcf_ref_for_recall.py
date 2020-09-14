import subprocess
def run_command(command):
    subprocess.check_call(command, shell=True)

for gt_conf_percentile, input_file, output_file in zip(snakemake.params.gt_conf_percentiles,
                                                       snakemake.input.singlesample_vcf_files_gt_conf_percentile_filtered,
                                                       snakemake.output.mutated_vcf_refs):
    run_command(f"python vcf_consensus_builder/cli.py -v {input_file} -r {snakemake.input.vcf_ref} -d {snakemake.input.empty_depth_file} "
                f"-o {output_file} --low-coverage 0 --no-coverage 0 -V")
    run_command(f"bwa index {output_file}")