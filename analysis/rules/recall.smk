rule make_recall_truth_probeset:
    input:
        truth1 = lambda wildcards: samples.xs(wildcards.sample1)["reference_assembly"],
        truth2 = lambda wildcards: samples.xs(wildcards.sample2)["reference_assembly"],
    output:
        probeset1 = "analysis/recall/truth_probesets/{sample1}/{sample1}_and_{sample2}.truth_probeset.fa",
        probeset2 = "analysis/recall/truth_probesets/{sample2}/{sample1}_and_{sample2}.truth_probeset.fa",
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
         truth_probeset = "analysis/recall/truth_probesets/{sample_id}/{filename_prefix}.truth_probeset.fa",
         variant_calls_probeset = "analysis/variant_calls_probesets/{sample_id}/{coverage}/{tool}.variant_calls_probeset.fa",
         variant_calls_probeset_index = "analysis/variant_calls_probesets/{sample_id}/{coverage}/{tool}.variant_calls_probeset.fa.amb"
    output:
         sam = "analysis/recall/map_probes/{sample_id}/{coverage}/{tool}/coverage_filter_{coverage_threshold}/strand_bias_filter_{strand_bias_threshold}/gaps_filter_{gaps_threshold}/{filename_prefix}.sam"
    threads: 2
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
        report = "analysis/recall/reports/{sample_id}/{coverage}/{tool}/coverage_filter_{coverage_threshold}/strand_bias_filter_{strand_bias_threshold}/gaps_filter_{gaps_threshold}/{filename_prefix}.report.tsv"
    threads: 1
    resources:
        mem_mb = lambda wildcards, attempt: 2000 * attempt
    log:
        "logs/create_recall_report_for_probe_mappings/{sample_id}/{coverage}/{tool}/coverage_filter_{coverage_threshold}/strand_bias_filter_{strand_bias_threshold}/gaps_filter_{gaps_threshold}/{filename_prefix}.report.log"
    script:
        "../scripts/create_recall_report_for_probe_mappings.py"

rule calculate_recall:
    input:
         recall_report_files_for_tool_and_coverage = lambda wildcards: tool_and_coverage_to_recall_report_files[wildcards.tool_and_coverage]
    output:
         recall_file_for_tool_and_coverage = "analysis/recall/recall_{tool_and_coverage}_gt_min_{min_gt}_step_{step_gt}_max_{max_gt}.tsv"
    threads: 1
    resources:
        mem_mb = lambda wildcards, attempt: 4000 * attempt
    log:
        "logs/calculate_recall/{tool_and_coverage}_gt_min_{min_gt}_step_{step_gt}_max_{max_gt}.log"
    script:
        "../scripts/calculate_recall.py"
