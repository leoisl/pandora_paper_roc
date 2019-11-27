rule make_recall_truth_probeset:
    input:
        truth1 = lambda wildcards: samples.xs(wildcards.sample1)["reference_assembly"],
        truth2 = lambda wildcards: samples.xs(wildcards.sample2)["reference_assembly"],
    output:
        probeset1 = output_folder + "/recall/truth_probesets/{sample1}/{sample1}_and_{sample2}.truth_probeset.fa",
        probeset2 = output_folder + "/recall/truth_probesets/{sample2}/{sample1}_and_{sample2}.truth_probeset.fa",
    params:
         flank_length = config["truth_probes_flank_length"]
    shadow:
        "shallow"
    threads: 1
    resources:
        mem_mb = lambda wildcards, attempt: 2000 * attempt
    log:
        "logs/make_recall_truth_probeset/{sample1}_and_{sample2}.log"
    script:
        "../scripts/make_recall_truth_probeset.py"


rule map_recall_truth_probes_to_variant_call_probes:
    input:
         truth_probeset = output_folder + "/recall/truth_probesets/{sample_id}/{filename_prefix}.truth_probeset.fa",
         variant_calls_probeset = output_folder + "/variant_calls_probesets/{sample_id}/{coverage}/{tool}/coverage_filter_{coverage_threshold}/strand_bias_filter_{strand_bias_threshold}/gaps_filter_{gaps_threshold}/variant_calls_probeset.fa",
         variant_calls_probeset_index = output_folder + "/variant_calls_probesets/{sample_id}/{coverage}/{tool}/coverage_filter_{coverage_threshold}/strand_bias_filter_{strand_bias_threshold}/gaps_filter_{gaps_threshold}/variant_calls_probeset.fa.amb",
    output:
         sam = output_folder + "/recall/map_probes/{sample_id}/{coverage}/{tool}/coverage_filter_{coverage_threshold}/strand_bias_filter_{strand_bias_threshold}/gaps_filter_{gaps_threshold}/{filename_prefix}.sam"
    threads: 4
    resources:
        mem_mb = lambda wildcards, attempt: 2000 * attempt
    log:
        "logs/map_recall_truth_probes_to_variant_call_probes/{sample_id}/{coverage}/{tool}/coverage_filter_{coverage_threshold}/strand_bias_filter_{strand_bias_threshold}/gaps_filter_{gaps_threshold}/{filename_prefix}.log"
    script:
        "../scripts/map_recall_truth_probes_to_variant_call_probes.py"

rule create_recall_report_for_probe_mappings:
    input:
        sam = rules.map_recall_truth_probes_to_variant_call_probes.output.sam,
        mask = lambda wildcards: samples.xs(wildcards.sample_id)["mask"]
    output:
        report = output_folder + "/recall/reports/{sample_id}/{coverage}/{tool}/coverage_filter_{coverage_threshold}/strand_bias_filter_{strand_bias_threshold}/gaps_filter_{gaps_threshold}/{filename_prefix}.report.tsv"
    threads: 1
    resources:
        mem_mb = lambda wildcards, attempt: 2000 * attempt
    log:
        "logs/create_recall_report_for_probe_mappings/{sample_id}/{coverage}/{tool}/coverage_filter_{coverage_threshold}/strand_bias_filter_{strand_bias_threshold}/gaps_filter_{gaps_threshold}/{filename_prefix}.report.log"
    script:
        "../scripts/create_recall_report_for_probe_mappings.py"

rule calculate_recall:
    input:
         recall_report_files_for_all_samples = lambda wildcards: cov_tool_and_filters_to_recall_report_files[wildcards.coverage, wildcards.tool, wildcards.coverage_threshold, wildcards.strand_bias_threshold, wildcards.gaps_threshold]
    output:
         recall_file_for_all_samples = output_folder + "/recall/recall_files/{coverage}/{tool}/coverage_filter_{coverage_threshold}/strand_bias_filter_{strand_bias_threshold}/gaps_filter_{gaps_threshold}/recall.tsv"
    params:
         number_of_points_in_ROC_curve = number_of_points_in_ROC_curve
    threads: 1
    resources:
        mem_mb = lambda wildcards, attempt: 4000 * attempt
    log:
        "logs/calculate_recall/{coverage}/{tool}/coverage_filter_{coverage_threshold}/strand_bias_filter_{strand_bias_threshold}/gaps_filter_{gaps_threshold}/recall.log"
    script:
        "../scripts/calculate_recall.py"
