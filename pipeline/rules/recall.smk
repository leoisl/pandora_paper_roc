rule make_empty_depth_file:
    input:
        file = "{file}"
    output:
        empty_depth_file = "{file}.depth"
    threads: 1
    resources:
        mem_mb = lambda wildcards, attempt: 100 * attempt
    log:
        "logs/make_empty_depth_file{file}.log"
    shell:
        "touch {output.empty_depth_file}"


rule gzip_vcf_file:
    input:
        vcf_file = "{filename}.vcf"
    output:
        gzipped_vcf_file = "{filename}.vcf.gz"
    threads: 1
    resources:
        mem_mb = lambda wildcards, attempt: 4000 * attempt
    log:
        "logs/gzip_vcf_file{filename}.log"
    shell:
        "bgzip -c {input.vcf_file} > {output.gzipped_vcf_file}"


rule index_gzipped_vcf_file:
    input:
        gzipped_vcf_file = rules.gzip_vcf_file.output.gzipped_vcf_file
    output:
        indexed_gzipped_vcf_file = "{filename}.vcf.gz.tbi"
    threads: 1
    resources:
        mem_mb = lambda wildcards, attempt: 4000 * attempt
    log:
        "logs/index_gzipped_vcf_file{filename}.log"
    shell:
        "tabix -p vcf {input.gzipped_vcf_file}"


rule make_vcf_for_a_single_sample:
    input:
        gzipped_multisample_vcf_file = rules.gzip_vcf_file.output.gzipped_vcf_file,
        indexed_gzipped_vcf_file = rules.index_gzipped_vcf_file.output.indexed_gzipped_vcf_file
    output:
        singlesample_vcf_file = "{filename}.vcf.sample_{sample_id}.vcf",
    threads: 1
    resources:
        mem_mb = lambda wildcards, attempt: 4000 * attempt
    log:
        "logs/make_vcf_for_a_single_sample{filename}_sample_{sample_id}.log"
    shell:
        "bcftools view -s {wildcards.sample_id} {input.gzipped_multisample_vcf_file} > {output.singlesample_vcf_file}"


rule filter_vcf_for_a_single_sample_by_gt_conf_percentile_for_pandora:
    input:
        gzipped_singlesample_vcf_file = "{filename}.vcf.sample_{sample_id}.vcf.gz",
        indexed_gzipped_singlesample_vcf_file = "{filename}.vcf.sample_{sample_id}.vcf.gz.tbi"
    output:
        singlesample_vcf_file_gt_conf_percentile_filtered = "{filename}.vcf.sample_{sample_id}.gt_conf_percentile_{gt_conf_percentile}.vcf"
    wildcard_constraints:
        filename=".*/pandora_multisample_genotyped"
    threads: 1
    resources:
        mem_mb = lambda wildcards, attempt: 4000 * attempt
    log:
        "logs/filter_vcf_for_a_single_sample_by_gt_conf_percentile_for_pandora{filename}_sample_{sample_id}.gt_conf_percentile_{gt_conf_percentile}.log"
    shell:
        "bash pipeline/scripts/filter_vcf_for_a_single_sample_by_gt_conf_percentile_for_pandora.sh {input.gzipped_singlesample_vcf_file} "
        "{wildcards.gt_conf_percentile} {output.singlesample_vcf_file_gt_conf_percentile_filtered}"
ruleorder: filter_vcf_for_a_single_sample_by_gt_conf_percentile_for_pandora > make_vcf_for_a_single_sample


rule filter_vcf_for_a_single_sample_by_gt_conf_percentile_for_snippy:
    input:
        gzipped_singlesample_vcf_file = "{filename}.vcf.sample_{sample_id}.vcf.gz",
        indexed_gzipped_singlesample_vcf_file = "{filename}.vcf.sample_{sample_id}.vcf.gz.tbi"
    output:
        singlesample_vcf_file_gt_conf_percentile_filtered = "{filename}.vcf.sample_{sample_id}.gt_conf_percentile_{gt_conf_percentile}.vcf"
    wildcard_constraints:
        filename=".*/snippy_[^/]+"
    threads: 1
    resources:
        mem_mb = lambda wildcards, attempt: 4000 * attempt
    log:
        "logs/filter_vcf_for_a_single_sample_by_gt_conf_percentile_for_snippy{filename}_sample_{sample_id}.gt_conf_percentile_{gt_conf_percentile}.log"
    shell:
         "bash pipeline/scripts/filter_vcf_for_a_single_sample_by_gt_conf_percentile_for_snippy.sh {input.gzipped_singlesample_vcf_file} "
         "{wildcards.gt_conf_percentile} {output.singlesample_vcf_file_gt_conf_percentile_filtered}"
ruleorder: filter_vcf_for_a_single_sample_by_gt_conf_percentile_for_snippy > make_vcf_for_a_single_sample


rule make_mutated_vcf_ref_for_recall:
    input:
         singlesample_vcf_file_gt_conf_percentile_filtered = lambda wildcards: f"{data.xs((wildcards.sample_id, wildcards.coverage, wildcards.tool))['vcf']}.sample_{wildcards.sample_id}.gt_conf_percentile_{wildcards.gt_conf_percentile}.vcf",
         vcf_ref = lambda wildcards: data.xs((wildcards.sample_id, wildcards.coverage, wildcards.tool))['vcf_reference'],
         empty_depth_file = lambda wildcards: f"{data.xs((wildcards.sample_id, wildcards.coverage, wildcards.tool))['vcf_reference']}.depth",
    output:
          mutated_vcf_ref = output_folder + "/recall/mutated_refs/{sample_id}/{coverage}/{tool}/coverage_filter_{coverage_threshold}/strand_bias_filter_{strand_bias_threshold}/gaps_filter_{gaps_threshold}/gt_conf_percentile_{gt_conf_percentile}/mutated_ref.fa"
    threads: 1
    resources:
        mem_mb = lambda wildcards, attempt: 4000 * attempt
    log:
        "logs/make_mutated_vcf_ref_for_recall/{sample_id}/{coverage}/{tool}/coverage_filter_{coverage_threshold}/strand_bias_filter_{strand_bias_threshold}/gaps_filter_{gaps_threshold}/gt_conf_percentile_{gt_conf_percentile}/mutated_ref.log"
    shell:
        "python vcf_consensus_builder/cli.py -v {input.singlesample_vcf_file_gt_conf_percentile_filtered} -r {input.vcf_ref} -d {input.empty_depth_file} "
        "-o {output.mutated_vcf_ref} --low-coverage 0 --no-coverage 0 -VVVV"


rule make_recall_truth_probeset:
    input:
        truth1 = lambda wildcards: samples.xs(wildcards.sample1)["reference_assembly"],
        truth2 = lambda wildcards: samples.xs(wildcards.sample2)["reference_assembly"],
    output:
        probeset1 = output_folder + "/recall/truth_probesets/{sample1}/{sample1}_and_{sample2}.truth_probeset.fa",
        probeset2 = output_folder + "/recall/truth_probesets/{sample2}/{sample1}_and_{sample2}.truth_probeset.fa",
    params:
         flank_length = config["truth_probes_flank_length_for_recall"]
    shadow:
        "shallow"
    threads: 1
    resources:
        mem_mb = lambda wildcards, attempt: 2000 * attempt
    log:
        "logs/make_recall_truth_probeset/{sample1}_and_{sample2}.log"
    script:
        "../scripts/make_recall_truth_probeset.py"


rule map_recall_truth_probeset_to_mutated_vcf_ref:
    input:
         truth_probeset = output_folder + "/recall/truth_probesets/{sample_id}/{sample_pair}.truth_probeset.fa",
         mutated_vcf_ref = rules.make_mutated_vcf_ref_for_recall.output.mutated_vcf_ref,
         mutated_vcf_ref_index = rules.make_mutated_vcf_ref_for_recall.output.mutated_vcf_ref + ".amb"
    output:
         sam = output_folder + "/recall/map_probes/{sample_id}/{coverage}/{tool}/coverage_filter_{coverage_threshold}/strand_bias_filter_{strand_bias_threshold}/gaps_filter_{gaps_threshold}/gt_conf_percentile_{gt_conf_percentile}/{sample_pair}.sam"
    threads: 4
    resources:
        mem_mb = lambda wildcards, attempt: 4000 * attempt
    log:
        "logs/map_recall_truth_probeset_to_mutated_vcf_ref/{sample_id}/{coverage}/{tool}/coverage_filter_{coverage_threshold}/strand_bias_filter_{strand_bias_threshold}/gaps_filter_{gaps_threshold}/gt_conf_percentile_{gt_conf_percentile}/{sample_pair}.log"
    script:
        "../scripts/map_recall_truth_variants_to_mutated_vcf_ref.py"


rule create_recall_report_for_truth_variants_mappings:
    input:
        sam = rules.map_recall_truth_probeset_to_mutated_vcf_ref.output.sam,
        mask = lambda wildcards: samples.xs(wildcards.sample_id)["mask"]
    output:
        report = output_folder + "/recall/reports/{sample_id}/{coverage}/{tool}/coverage_filter_{coverage_threshold}/strand_bias_filter_{strand_bias_threshold}/gaps_filter_{gaps_threshold}/gt_conf_percentile_{gt_conf_percentile}/{sample_pair}.report.tsv"
    threads: 1
    resources:
        mem_mb = lambda wildcards, attempt: 2000 * attempt
    log:
        "logs/create_recall_report_for_truth_variants_mappings/{sample_id}/{coverage}/{tool}/coverage_filter_{coverage_threshold}/strand_bias_filter_{strand_bias_threshold}/gaps_filter_{gaps_threshold}/gt_conf_percentile_{gt_conf_percentile}/{sample_pair}.report.log"
    script:
        "../scripts/create_recall_report_for_probe_mappings.py"


rule calculate_recall:
    input:
         recall_report_files_for_all_samples_and_all_gt_conf_percentile = lambda wildcards: cov_tool_and_filters_to_recall_report_files[(wildcards.coverage, wildcards.tool, wildcards.coverage_threshold, wildcards.strand_bias_threshold, wildcards.gaps_threshold)]
    output:
         recall_file_for_all_samples_and_all_gt_conf_percentile = output_folder + "/recall/recall_files/{coverage}/{tool}/coverage_filter_{coverage_threshold}/strand_bias_filter_{strand_bias_threshold}/gaps_filter_{gaps_threshold}/recall.tsv"
    params:
         number_of_points_in_ROC_curve = 100 # from 0 to 100
    threads: 1
    resources:
        mem_mb = lambda wildcards, attempt: 16000 * attempt
    log:
        "logs/calculate_recall/{coverage}/{tool}/coverage_filter_{coverage_threshold}/strand_bias_filter_{strand_bias_threshold}/gaps_filter_{gaps_threshold}/recall.log"
    script:
        "../scripts/calculate_recall.py"
