rule map_variant_call_probeset_to_reference_assembly:
    input:
        variant_call_probeset = rules.make_variant_calls_probeset.output.probeset,
        reference_assembly = lambda wildcards: data.xs((wildcards.sample_id, wildcards.coverage, wildcards.tool))["reference_assembly"]
    output:
          variant_call_probeset_mapped_to_ref = "analysis/precision/variant_calls_probesets_mapped_to_refs/{sample_id}/{coverage}/{tool}/coverage_filter_{coverage_threshold}/strand_bias_filter_{strand_bias_threshold}/gaps_filter_{gaps_threshold}/variant_calls_probeset_mapped.sam"
    threads: 1
    resources:
        mem_mb = lambda wildcards, attempt: 4000 * attempt
    log:
        "logs/map_variant_call_probeset_to_reference_assembly/{sample_id}_{coverage}_{tool}_coverage_filter_{coverage_threshold}_strand_bias_filter_{strand_bias_threshold}_gaps_filter_{gaps_threshold}.log"
    script:
        "../scripts/map_variant_call_probeset_to_reference_assembly.py"


rule create_precision_report_from_probe_mappings:
    input:
        variant_call_probeset_mapped_to_ref = rules.map_variant_call_probeset_to_reference_assembly.output.variant_call_probeset_mapped_to_ref,
        mask = lambda wildcards: data.xs((wildcards.sample_id, wildcards.coverage, wildcards.tool))["mask"]
    output:
        variant_call_precision_report = "analysis/precision/reports_from_probe_mappings/{sample_id}/{coverage}/{tool}/variant_calls_probeset_report.tsv"
    threads: 1
    resources:
        mem_mb = lambda wildcards, attempt: 1000 * attempt
    log:
        "logs/create_precision_report_from_probe_mappings/{sample_id}_{coverage}_{tool}.log"
    script:
        "../scripts/create_precision_report_from_probe_mappings.py"


rule calculate_precision:
    input:
         precision_report_files_for_tool_and_coverage = lambda wildcards: tool_and_coverage_to_precision_report_files[wildcards.tool_and_coverage]
    output:
         precision_file_for_tool_and_coverage = "analysis/precision/precision_{tool_and_coverage}_gt_min_{min_gt}_step_{step_gt}_max_{max_gt}.tsv"
    threads: 1
    resources:
        mem_mb = lambda wildcards, attempt: 4000 * attempt
    log:
        "logs/calculate_precision/{tool_and_coverage}_gt_min_{min_gt}_step_{step_gt}_max_{max_gt}.log"
    script:
        "../scripts/calculate_precision.py"
