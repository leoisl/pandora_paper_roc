rule map_variant_call_probeset_to_reference_assembly:
    input:
        variant_call_probeset = rules.make_variant_calls_probeset.output.probeset,
        reference_assembly = lambda wildcards: data.xs((wildcards.sample_id, wildcards.coverage, wildcards.tool))["reference_assembly"]
    output:
          variant_call_probeset_mapped_to_ref = "analysis/precision/variant_calls_probesets_mapped_to_refs/{sample_id}/{coverage}/{tool}/variant_calls_probeset_mapped.sam"
    threads: 1
    resources:
        mem_mb = lambda wildcards, attempt: 4000 * attempt
    log:
        "logs/map_variant_call_probeset_to_reference_assembly/{sample_id}_{coverage}_{tool}.log"
    script:
        "../scripts/map_variant_call_probeset_to_reference_assembly.py"


rule create_precision_report_from_probe_mappings:
    input:
        variant_call_probeset_mapped_to_ref = rules.map_variant_call_probeset_to_reference_assembly.output.variant_call_probeset_mapped_to_ref,
        mask = lambda wildcards: data.xs((wildcards.sample_id, wildcards.coverage, wildcards.tool))["mask"]
    output:
        variant_call_precision_report = "analysis/precision/reports_from_probe_mappings/{sample_id}/{coverage}/{tool}/variant_calls_probeset_report"
    threads: 1
    resources:
        mem_mb = lambda wildcards, attempt: 1000 * attempt
    log:
        "logs/create_precision_report_from_probe_mappings/{sample_id}_{coverage}_{tool}.log"
    script:
        "../scripts/create_precision_report_from_probe_mappings.py"
