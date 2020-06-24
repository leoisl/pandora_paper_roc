import subprocess
from pipeline.scripts.utils import get_sample_pairs_containing_given_sample

def run_command(command):
    subprocess.check_call(command, shell=True)

rule make_vcf_for_a_single_sample:
    input:
        gzipped_multisample_vcf_file = rules.gzip_vcf_file.output.gzipped_vcf_file,
        indexed_gzipped_vcf_file = rules.index_gzipped_vcf_file.output.indexed_gzipped_vcf_file
    output:
        singlesample_vcf_file = "{filename}.vcf.sample_{sample_id}.vcf",
    threads: 1
    resources:
        mem_mb = lambda wildcards, attempt: 2000 * attempt
    log:
        "logs/make_vcf_for_a_single_sample{filename}_sample_{sample_id}.log"
    shell:
        "bcftools view -s {wildcards.sample_id} {input.gzipped_multisample_vcf_file} > {output.singlesample_vcf_file}"


rule filter_vcf_for_a_single_sample_by_gt_conf_percentile_for_pandora:
    input:
        gzipped_singlesample_vcf_file = "{filename}.vcf.sample_{sample_id}.vcf.gz",
        indexed_gzipped_singlesample_vcf_file = "{filename}.vcf.sample_{sample_id}.vcf.gz.tbi"
    output:
        singlesample_vcf_files_gt_conf_percentile_filtered = expand("{{filename}}.vcf.sample_{{sample_id}}.gt_conf_percentile_{gt_conf_percentile}.vcf", gt_conf_percentile=gt_conf_percentiles)
    wildcard_constraints:
        filename=".*/pandora_multisample_genotyped_.*\.vcf\.\~\~vcf\~\~fixed\~\~"
    threads: 1
    resources:
        mem_mb = lambda wildcards, attempt: 2000 * attempt
    log:
        "logs/filter_vcf_for_a_single_sample_by_gt_conf_percentile_for_pandora{filename}_sample_{sample_id}.log"
    run:
        for gt_conf_percentile, output_file in zip(gt_conf_percentiles, output.singlesample_vcf_files_gt_conf_percentile_filtered):
            run_command(f"bash pipeline/scripts/filter_vcf_for_a_single_sample_by_gt_conf_percentile_for_pandora.sh "
                        f"{input.gzipped_singlesample_vcf_file} "
                        f"{gt_conf_percentile} "
                        f"{output_file}")
ruleorder: filter_vcf_for_a_single_sample_by_gt_conf_percentile_for_pandora > make_vcf_for_a_single_sample


rule filter_vcf_for_a_single_sample_by_gt_conf_percentile_for_snippy:
    input:
        gzipped_singlesample_vcf_file = "{filename}.vcf.sample_{sample_id}.vcf.gz",
        indexed_gzipped_singlesample_vcf_file = "{filename}.vcf.sample_{sample_id}.vcf.gz.tbi"
    output:
        singlesample_vcf_files_gt_conf_percentile_filtered = expand("{{filename}}.vcf.sample_{{sample_id}}.gt_conf_percentile_{gt_conf_percentile}.vcf", gt_conf_percentile=gt_conf_percentiles)
    wildcard_constraints:
        filename=".*/snippy_[^/]+\.vcf\.\~\~vcf\~\~fixed\~\~"
    threads: 1
    resources:
        mem_mb = lambda wildcards, attempt: 2000 * attempt
    log:
        "logs/filter_vcf_for_a_single_sample_by_gt_conf_percentile_for_snippy{filename}_sample_{sample_id}.log"
    run:
        for gt_conf_percentile, output_file in zip(gt_conf_percentiles, output.singlesample_vcf_files_gt_conf_percentile_filtered):
            run_command(f"bash pipeline/scripts/filter_vcf_for_a_single_sample_by_gt_conf_percentile_for_snippy.sh "
                        f"{input.gzipped_singlesample_vcf_file} "
                        f"{gt_conf_percentile} "
                        f"{output_file}")
ruleorder: filter_vcf_for_a_single_sample_by_gt_conf_percentile_for_snippy > make_vcf_for_a_single_sample


rule filter_vcf_for_a_single_sample_by_gt_conf_percentile_for_samtools:
    input:
        gzipped_singlesample_vcf_file = "{filename}.vcf.sample_{sample_id}.vcf.gz",
        indexed_gzipped_singlesample_vcf_file = "{filename}.vcf.sample_{sample_id}.vcf.gz.tbi"
    output:
        singlesample_vcf_files_gt_conf_percentile_filtered = expand("{{filename}}.vcf.sample_{{sample_id}}.gt_conf_percentile_{gt_conf_percentile}.vcf", gt_conf_percentile=gt_conf_percentiles)
    wildcard_constraints:
        filename=".*/samtools_[^/]+\.vcf\.\~\~vcf\~\~fixed\~\~"
    threads: 1
    resources:
        mem_mb = lambda wildcards, attempt: 2000 * attempt
    log:
        "logs/filter_vcf_for_a_single_sample_by_gt_conf_percentile_for_samtools{filename}_sample_{sample_id}.log"
    run:
        for gt_conf_percentile, output_file in zip(gt_conf_percentiles, output.singlesample_vcf_files_gt_conf_percentile_filtered):
            run_command(f"bash pipeline/scripts/filter_vcf_for_a_single_sample_by_gt_conf_percentile_for_samtools.sh "
                        f"{input.gzipped_singlesample_vcf_file} "
                        f"{gt_conf_percentile} "
                        f"{output_file}")
ruleorder: filter_vcf_for_a_single_sample_by_gt_conf_percentile_for_samtools > make_vcf_for_a_single_sample


rule make_mutated_vcf_ref_for_recall:
    input:
         singlesample_vcf_files_gt_conf_percentile_filtered = lambda wildcards: expand(f"{data.xs((wildcards.sample_id, wildcards.coverage, wildcards.tool))['vcf']}.sample_{wildcards.sample_id}.gt_conf_percentile_{{gt_conf_percentile}}.vcf", gt_conf_percentile=gt_conf_percentiles),
         vcf_ref = lambda wildcards: data.xs((wildcards.sample_id, wildcards.coverage, wildcards.tool))['vcf_reference'],
         empty_depth_file = lambda wildcards: f"{data.xs((wildcards.sample_id, wildcards.coverage, wildcards.tool))['vcf_reference']}.depth",
    output:
          mutated_vcf_refs = expand(output_folder + "/recall/mutated_refs/{{sample_id}}/{{coverage}}/{{tool}}/coverage_filter_{{coverage_threshold}}/strand_bias_filter_{{strand_bias_threshold}}/gaps_filter_{{gaps_threshold}}/gt_conf_percentile_{gt_conf_percentile}/mutated_ref.fa", gt_conf_percentile=gt_conf_percentiles),
          indexes = expand(output_folder + "/recall/mutated_refs/{{sample_id}}/{{coverage}}/{{tool}}/coverage_filter_{{coverage_threshold}}/strand_bias_filter_{{strand_bias_threshold}}/gaps_filter_{{gaps_threshold}}/gt_conf_percentile_{gt_conf_percentile}/mutated_ref.fa.amb", gt_conf_percentile=gt_conf_percentiles)
    threads: 1
    resources:
        mem_mb = lambda wildcards, attempt: 4000 * attempt
    log:
        "logs/make_mutated_vcf_ref_for_recall/{sample_id}/{coverage}/{tool}/coverage_filter_{coverage_threshold}/strand_bias_filter_{strand_bias_threshold}/gaps_filter_{gaps_threshold}/mutated_ref.log"
    run:
        for gt_conf_percentile, input_file, output_file in zip(gt_conf_percentiles, input.singlesample_vcf_files_gt_conf_percentile_filtered, output.mutated_vcf_refs):
            run_command(f"python vcf_consensus_builder/cli.py -v {input_file} -r {input.vcf_ref} -d {input.empty_depth_file} "
                        f"-o {output_file} --low-coverage 0 --no-coverage 0 -V")
            run_command(f"bwa index {output_file}")


rule map_recall_truth_probeset_to_mutated_vcf_ref:
    input:
         truth_probeset = deduplicated_variants_output_folder + "/truth_probesets/{sample_id}/{sample_pair}.truth_probeset.fa",
         mutated_vcf_refs = rules.make_mutated_vcf_ref_for_recall.output.mutated_vcf_refs,
    output:
         sams = expand(output_folder + "/recall/map_probes/{{sample_id}}/{{coverage}}/{{tool}}/coverage_filter_{{coverage_threshold}}/strand_bias_filter_{{strand_bias_threshold}}/gaps_filter_{{gaps_threshold}}/gt_conf_percentile_{gt_conf_percentile}/{{sample_pair}}.sam", gt_conf_percentile=gt_conf_percentiles)
    threads: 1
    resources:
        mem_mb = 4000
    log:
        "logs/map_recall_truth_probeset_to_mutated_vcf_ref/{sample_id}/{coverage}/{tool}/coverage_filter_{coverage_threshold}/strand_bias_filter_{strand_bias_threshold}/gaps_filter_{gaps_threshold}/{sample_pair}.log"
    script:
        "../scripts/map_recall_truth_variants_to_mutated_vcf_ref.py"


rule create_recall_report_for_truth_variants_mappings:
    input:
        sams = rules.map_recall_truth_probeset_to_mutated_vcf_ref.output.sams,
        mask = lambda wildcards: samples.xs(wildcards.sample_id)["mask"]
    output:
        reports = expand(output_folder + "/recall/reports/{{sample_id}}/{{coverage}}/{{tool}}/coverage_filter_{{coverage_threshold}}/strand_bias_filter_{{strand_bias_threshold}}/gaps_filter_{{gaps_threshold}}/gt_conf_percentile_{gt_conf_percentile}/{{sample_pair}}.report.tsv", gt_conf_percentile=gt_conf_percentiles)
    params:
        gt_conf_percentiles = gt_conf_percentiles
    threads: 1
    resources:
        mem_mb = lambda wildcards, attempt: 4000 * attempt
    log:
        "logs/create_recall_report_for_truth_variants_mappings/{sample_id}/{coverage}/{tool}/coverage_filter_{coverage_threshold}/strand_bias_filter_{strand_bias_threshold}/gaps_filter_{gaps_threshold}/{sample_pair}.report.log"
    script:
        "../scripts/create_recall_report_for_probe_mappings.py"


rule create_recall_report_per_sample_for_calculator:
    input:
         recall_report_files_for_one_sample_and_all_gt_conf_percentiles = lambda wildcards: sample_cov_tool_and_filters_to_recall_report_files[(wildcards.sample, wildcards.coverage, wildcards.tool, wildcards.coverage_threshold, wildcards.strand_bias_threshold, wildcards.gaps_threshold)]
    output:
         recall_report_per_sample_for_calculator = output_folder + "/recall/recall_report_per_sample_for_calculator/{sample}/{coverage}/{tool}/coverage_filter_{coverage_threshold}/strand_bias_filter_{strand_bias_threshold}/gaps_filter_{gaps_threshold}/recall_report_per_sample_for_calculator.tsv"
    params:
         gt_conf_percentiles = gt_conf_percentiles
    threads: 1
    resources:
        mem_mb = lambda wildcards, attempt: 8000 * attempt
    log:
        "logs/create_recall_report_per_sample_for_calculator/{sample}/{coverage}/{tool}/coverage_filter_{coverage_threshold}/strand_bias_filter_{strand_bias_threshold}/gaps_filter_{gaps_threshold}/create_recall_report_per_sample_for_calculator.log"
    script:
        "../scripts/create_recall_report_per_sample_for_calculator.py"



rule calculate_recall:
    input:
         recall_report_per_sample_for_calculator = expand(output_folder + "/recall/recall_report_per_sample_for_calculator/{sample}/{{coverage}}/{{tool}}/coverage_filter_{{coverage_threshold}}/strand_bias_filter_{{strand_bias_threshold}}/gaps_filter_{{gaps_threshold}}/recall_report_per_sample_for_calculator.tsv", sample=samples["sample_id"])
    output:
         recall_file_for_all_samples_and_all_gt_conf_percentile = output_folder + "/recall/recall_files/{coverage}/{tool}/coverage_filter_{coverage_threshold}/strand_bias_filter_{strand_bias_threshold}/gaps_filter_{gaps_threshold}/recall.tsv",
    params:
         gt_conf_percentiles = gt_conf_percentiles
    threads: 1
    resources:
        mem_mb = lambda wildcards, attempt: 8000 * attempt
    log:
        "logs/calculate_recall/{coverage}/{tool}/coverage_filter_{coverage_threshold}/strand_bias_filter_{strand_bias_threshold}/gaps_filter_{gaps_threshold}/recall.log"
    script:
        "../scripts/calculate_recall.py"


rule calculate_recall_per_sample_no_gt_conf_filter:
    input:
         all_recall_reports_for_one_sample_with_no_gt_conf_filter = lambda wildcards: expand(output_folder + "/recall/reports/{{sample_id}}/{{coverage}}/{{tool}}/coverage_filter_{{coverage_threshold}}/strand_bias_filter_{{strand_bias_threshold}}/gaps_filter_{{gaps_threshold}}/gt_conf_percentile_0/{sample_pairs}.report.tsv", sample_pairs = get_sample_pairs_containing_given_sample(sample_pairs, wildcards.sample_id))
    output:
         recall_file_for_one_sample_with_no_gt_conf_filter = output_folder + "/recall/recall_files_per_sample/{sample_id}/{coverage}/{tool}/coverage_filter_{coverage_threshold}/strand_bias_filter_{strand_bias_threshold}/gaps_filter_{gaps_threshold}/recall.tsv"
    threads: 1
    resources:
        mem_mb = lambda wildcards, attempt: 4000 * attempt
    log:
        "logs/calculate_recall_per_sample_no_gt_conf_filter/{sample_id}/{coverage}/{tool}/coverage_filter_{coverage_threshold}/strand_bias_filter_{strand_bias_threshold}/gaps_filter_{gaps_threshold}/recall.log"
    script:
        "../scripts/calculate_recall_per_sample_no_gt_conf_filter.py"


# TODO: put this back? removed because it yields 15k jobs, and it is not useful right now
# rule calculate_recall_per_sample_pair_no_gt_conf_filter:
#     input:
#          all_recall_reports_for_one_sample_pair_with_no_gt_conf_filter = output_folder + "/recall/reports/{sample_id}/{coverage}/{tool}/coverage_filter_{coverage_threshold}/strand_bias_filter_{strand_bias_threshold}/gaps_filter_{gaps_threshold}/gt_conf_percentile_0/{sample_pair}.report.tsv",
#          aligned_bases_percentage = deduplicated_variants_output_folder + "/dnadiff_reports/{sample_id}/{sample_pair}.aligned_bases_percentage",
#     output:
#          recall_file_for_one_sample_pair_with_no_gt_conf_filter = output_folder + "/recall/recall_files_per_sample_pair/{sample_id}/{coverage}/{tool}/coverage_filter_{coverage_threshold}/strand_bias_filter_{strand_bias_threshold}/gaps_filter_{gaps_threshold}/{sample_pair}.recall.tsv"
#     threads: 1
#     resources:
#         mem_mb = lambda wildcards, attempt: 4000 * attempt
#     log:
#         "logs/calculate_recall_per_sample_pair_no_gt_conf_filter/{sample_id}/{coverage}/{tool}/coverage_filter_{coverage_threshold}/strand_bias_filter_{strand_bias_threshold}/gaps_filter_{gaps_threshold}/{sample_pair}.recall.log"
#     script:
#         "../scripts/calculate_recall_per_sample_pair_no_gt_conf_filter.py"


rule calculate_recall_per_number_of_samples_no_gt_conf_filter:
    input:
         all_recall_reports_with_no_gt_conf_filter = lambda wildcards: cov_tool_and_filters_to_recall_reports_with_no_gt_conf_filter[(
             wildcards.coverage, wildcards.tool, wildcards.coverage_threshold, wildcards.strand_bias_threshold, wildcards.gaps_threshold
         )]
    output:
         recall_per_number_of_samples = output_folder + "/recall/recall_per_number_of_samples/{coverage}/{tool}/coverage_filter_{coverage_threshold}/strand_bias_filter_{strand_bias_threshold}/gaps_filter_{gaps_threshold}/recall_per_number_of_samples.csv"
    params:
         list_with_number_of_samples = list_with_number_of_samples
    threads: 1
    resources:
        mem_mb = lambda wildcards, attempt: 12000 * attempt
    log:
        "logs/calculate_recall_per_number_of_samples_no_gt_conf_filter/{coverage}/{tool}/coverage_filter_{coverage_threshold}/strand_bias_filter_{strand_bias_threshold}/gaps_filter_{gaps_threshold}/recall_per_number_of_samples.log"
    script:
        "../scripts/calculate_recall_per_number_of_samples_no_gt_conf_filter.py"

