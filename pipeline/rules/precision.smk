rule make_variant_calls_probeset_for_precision:
    input:
         vcf = lambda wildcards: data.xs((wildcards.sample_id, wildcards.coverage, wildcards.tool))["vcf"],
         vcf_ref = lambda wildcards: data.xs((wildcards.sample_id, wildcards.coverage, wildcards.tool))["vcf_reference"]
    output:
          probeset = output_folder + "/precision/variant_calls_probesets/{sample_id}/{coverage}/{tool}/coverage_filter_{coverage_threshold}/strand_bias_filter_{strand_bias_threshold}/gaps_filter_{gaps_threshold}/variant_calls_probeset.fa"
    params:
          flank_length = config["variant_calls_flank_length_for_precision"]
    threads: 1
    resources:
        mem_mb = lambda wildcards, attempt: 1000 * attempt
    log:
        "logs/make_variant_calls_probeset_for_precision/{sample_id}/{coverage}/{tool}/coverage_filter_{coverage_threshold}/strand_bias_filter_{strand_bias_threshold}/gaps_filter_{gaps_threshold}/variant_calls_probeset.log"
    script:
        "../scripts/make_variant_calls_probeset.py"


rule map_variant_call_probeset_to_reference_assembly:
    input:
        variant_call_probeset = rules.make_variant_calls_probeset_for_precision.output.probeset,
        reference_assembly = lambda wildcards: data.xs((wildcards.sample_id, wildcards.coverage, wildcards.tool))["reference_assembly"],
        reference_assembly_index = lambda wildcards: data.xs((wildcards.sample_id, wildcards.coverage, wildcards.tool))["reference_assembly"]+".amb"
    output:
          variant_call_probeset_mapped_to_ref = output_folder + "/precision/variant_calls_probesets_mapped_to_refs/{sample_id}/{coverage}/{tool}/coverage_filter_{coverage_threshold}/strand_bias_filter_{strand_bias_threshold}/gaps_filter_{gaps_threshold}/variant_calls_probeset_mapped.sam"
    threads: 4
    resources:
        mem_mb = lambda wildcards, attempt: 4000 * attempt
    log:
        "logs/map_variant_call_probeset_to_reference_assembly/{sample_id}/{coverage}/{tool}/coverage_filter_{coverage_threshold}/strand_bias_filter_{strand_bias_threshold}/gaps_filter_{gaps_threshold}/variant_calls_probeset_mapped.log"
    script:
        "../scripts/map_variant_call_probeset_to_reference_assembly.py"


rule create_precision_report_from_probe_mappings:
    input:
        variant_call_probeset_mapped_to_ref = rules.map_variant_call_probeset_to_reference_assembly.output.variant_call_probeset_mapped_to_ref,
        mask = lambda wildcards: data.xs((wildcards.sample_id, wildcards.coverage, wildcards.tool))["mask"]
    output:
        variant_call_precision_report = output_folder + "/precision/reports_from_probe_mappings/{sample_id}/{coverage}/{tool}/coverage_filter_{coverage_threshold}/strand_bias_filter_{strand_bias_threshold}/gaps_filter_{gaps_threshold}/variant_calls_probeset_report.tsv",
        nb_of_records_removed_with_mapq_sam_records_filter_filepath = output_folder + "/precision/reports_from_probe_mappings/{sample_id}/{coverage}/{tool}/coverage_filter_{coverage_threshold}/strand_bias_filter_{strand_bias_threshold}/gaps_filter_{gaps_threshold}/nb_of_records_removed_with_mapq_sam_records_filter.csv"
    threads: 1
    resources:
        mem_mb = lambda wildcards, attempt: 1000 * attempt
    log:
        "logs/create_precision_report_from_probe_mappings/{sample_id}/{coverage}/{tool}/coverage_filter_{coverage_threshold}/strand_bias_filter_{strand_bias_threshold}/gaps_filter_{gaps_threshold}/variant_calls_probeset_report.log"
    script:
        "../scripts/create_precision_report_from_probe_mappings.py"


rule calculate_precision:
    input:
         precision_report_files_for_all_samples = lambda wildcards: cov_tool_and_filters_to_precision_report_files[wildcards.coverage, wildcards.tool, wildcards.coverage_threshold, wildcards.strand_bias_threshold, wildcards.gaps_threshold]
    output:
         precision_file_for_all_samples = output_folder + "/precision/precision_files/{coverage}/{tool}/coverage_filter_{coverage_threshold}/strand_bias_filter_{strand_bias_threshold}/gaps_filter_{gaps_threshold}/precision.tsv"
    params:
         number_of_points_in_ROC_curve = number_of_points_in_ROC_curve,
         tool = lambda wildcards: wildcards.tool
    threads: 1
    resources:
        mem_mb = lambda wildcards, attempt: 4000 * attempt
    log:
        "logs/calculate_precision/{coverage}/{tool}/coverage_filter_{coverage_threshold}/strand_bias_filter_{strand_bias_threshold}/gaps_filter_{gaps_threshold}/precision.log"
    script:
        "../scripts/calculate_precision.py"

def get_precision_report_files_for_all_samples_for_all_snippy(wildcards):
    return list(itertools.chain.from_iterable(
        cov_tool_and_filters_to_precision_report_files[wildcards.coverage, snippy_tool, wildcards.coverage_threshold, \
            wildcards.strand_bias_threshold, wildcards.gaps_threshold] \
            for snippy_tool in all_tools if snippy_tool.startswith("snippy")))

rule calculate_precision_for_all_snippy:
    input:
         precision_report_files_for_all_samples = lambda wildcards: get_precision_report_files_for_all_samples_for_all_snippy(wildcards)
    output:
         precision_file_for_all_samples = output_folder + "/precision/precision_files/{coverage}/snippy_all_curves/coverage_filter_{coverage_threshold}/strand_bias_filter_{strand_bias_threshold}/gaps_filter_{gaps_threshold}/precision.tsv"
    params:
         number_of_points_in_ROC_curve = number_of_points_in_ROC_curve,
         tool = "snippy_all_curves"
    threads: 1
    resources:
        mem_mb = lambda wildcards, attempt: {1: 40000, 2: 60000, 3: 80000, 4: 120000}
    log:
        "logs/calculate_precision/{coverage}/snippy_all_curves/coverage_filter_{coverage_threshold}/strand_bias_filter_{strand_bias_threshold}/gaps_filter_{gaps_threshold}/precision.log"
    script:
        "../scripts/calculate_precision.py"


ruleorder: calculate_precision_for_all_snippy > calculate_precision