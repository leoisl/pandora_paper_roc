rule make_variant_calls_probeset_for_recall:
    input:
         vcf = lambda wildcards: data.xs((wildcards.sample_id, wildcards.coverage, wildcards.tool))["vcf"],
         vcf_ref = lambda wildcards: data.xs((wildcards.sample_id, wildcards.coverage, wildcards.tool))["vcf_reference"]
    output:
          probeset = output_folder + "/recall/variant_calls_probesets/{sample_id}/{coverage}/{tool}/coverage_filter_{coverage_threshold}/strand_bias_filter_{strand_bias_threshold}/gaps_filter_{gaps_threshold}/variant_calls_probeset.fa"
    params:
          flank_length = config["variant_calls_flank_length_for_recall"]
    threads: 1
    resources:
        mem_mb = lambda wildcards, attempt: 1000 * attempt
    log:
        "logs/make_variant_calls_probeset_for_recall/{sample_id}/{coverage}/{tool}/coverage_filter_{coverage_threshold}/strand_bias_filter_{strand_bias_threshold}/gaps_filter_{gaps_threshold}/variant_calls_probeset.log"
    script:
        "../scripts/make_variant_calls_probeset.py"


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


rule restrict_recall_truth_probeset_to_probes_that_map_uniquely_to_the_truth:
    input:
        truth = lambda wildcards: samples.xs(wildcards.sample_id)["reference_assembly"],
        unrestricted_probeset = output_folder + "/recall/truth_probesets/{sample_id}/{sample_pair}.truth_probeset.fa"
    output:
        unrestricted_probeset_mapped_to_truth_sam = output_folder + "/recall/restricted_truth_probesets/{sample_id}/{sample_pair}.unrestricted_truth_probeset.mapped_to_truth.sam",
        restricted_probeset = output_folder + "/recall/restricted_truth_probesets/{sample_id}/{sample_pair}.restricted_truth_probeset.fa",
        nb_of_truth_probes_removed_with_unique_sam_records_filter_filepath = output_folder + "/recall/restricted_truth_probesets/{sample_id}/{sample_pair}.nb_of_truth_probes_removed_with_unique_sam_records_filter.csv"
    threads: 4
    resources:
        mem_mb = lambda wildcards, attempt: 2000 * attempt
    log:
        "logs/restrict_recall_truth_probeset_to_probes_that_map_uniquely_to_the_truth/{sample_id}_{sample_pair}.log"
    script:
        "../scripts/restrict_recall_truth_probeset_to_probes_that_map_uniquely_to_the_truth.py"


rule map_recall_truth_probes_to_variant_call_probes:
    input:
         truth_probeset = output_folder + "/recall/restricted_truth_probesets/{sample_id}/{sample_pair}.restricted_truth_probeset.fa",
         variant_calls_probeset = output_folder + "/recall/variant_calls_probesets/{sample_id}/{coverage}/{tool}/coverage_filter_{coverage_threshold}/strand_bias_filter_{strand_bias_threshold}/gaps_filter_{gaps_threshold}/variant_calls_probeset.fa",
         variant_calls_probeset_index = output_folder + "/recall/variant_calls_probesets/{sample_id}/{coverage}/{tool}/coverage_filter_{coverage_threshold}/strand_bias_filter_{strand_bias_threshold}/gaps_filter_{gaps_threshold}/variant_calls_probeset.fa.amb",
    output:
         sam = output_folder + "/recall/map_probes/{sample_id}/{coverage}/{tool}/coverage_filter_{coverage_threshold}/strand_bias_filter_{strand_bias_threshold}/gaps_filter_{gaps_threshold}/{sample_pair}.sam"
    threads: 4
    resources:
        mem_mb = lambda wildcards, attempt: 2000 * attempt
    log:
        "logs/map_recall_truth_probes_to_variant_call_probes/{sample_id}/{coverage}/{tool}/coverage_filter_{coverage_threshold}/strand_bias_filter_{strand_bias_threshold}/gaps_filter_{gaps_threshold}/{sample_pair}.log"
    script:
        "../scripts/map_recall_truth_probes_to_variant_call_probes.py"

rule create_recall_report_for_probe_mappings:
    input:
        sam = rules.map_recall_truth_probes_to_variant_call_probes.output.sam,
        mask = lambda wildcards: samples.xs(wildcards.sample_id)["mask"]
    output:
        report = output_folder + "/recall/reports/{sample_id}/{coverage}/{tool}/coverage_filter_{coverage_threshold}/strand_bias_filter_{strand_bias_threshold}/gaps_filter_{gaps_threshold}/{sample_pair}.report.tsv"
    threads: 1
    resources:
        mem_mb = lambda wildcards, attempt: 2000 * attempt
    log:
        "logs/create_recall_report_for_probe_mappings/{sample_id}/{coverage}/{tool}/coverage_filter_{coverage_threshold}/strand_bias_filter_{strand_bias_threshold}/gaps_filter_{gaps_threshold}/{sample_pair}.report.log"
    script:
        "../scripts/create_recall_report_for_probe_mappings.py"

rule calculate_recall:
    input:
         recall_report_files_for_all_samples = lambda wildcards: cov_tool_and_filters_to_recall_report_files[wildcards.coverage, wildcards.tool, wildcards.coverage_threshold, wildcards.strand_bias_threshold, wildcards.gaps_threshold]
    output:
         recall_file_for_all_samples = output_folder + "/recall/recall_files/{coverage}/{tool}/coverage_filter_{coverage_threshold}/strand_bias_filter_{strand_bias_threshold}/gaps_filter_{gaps_threshold}/recall.tsv"
    params:
         number_of_points_in_ROC_curve = number_of_points_in_ROC_curve,
         tool = lambda wildcards: wildcards.tool
    threads: 1
    resources:
        mem_mb = lambda wildcards, attempt: 4000 * attempt
    log:
        "logs/calculate_recall/{coverage}/{tool}/coverage_filter_{coverage_threshold}/strand_bias_filter_{strand_bias_threshold}/gaps_filter_{gaps_threshold}/recall.log"
    script:
        "../scripts/calculate_recall.py"


def get_recall_report_files_for_all_samples_for_all_snippy(wildcards):
    return list(itertools.chain.from_iterable(
        cov_tool_and_filters_to_recall_report_files[wildcards.coverage, snippy_tool, wildcards.coverage_threshold, \
            wildcards.strand_bias_threshold, wildcards.gaps_threshold] \
            for snippy_tool in all_tools if snippy_tool.startswith("snippy")))

rule calculate_recall_for_all_snippy:
    input:
         recall_report_files_for_all_samples = \
            lambda wildcards: get_recall_report_files_for_all_samples_for_all_snippy(wildcards)
    output:
         recall_file_for_all_samples = output_folder + "/recall/recall_files/{coverage}/snippy_all_curves/coverage_filter_{coverage_threshold}/strand_bias_filter_{strand_bias_threshold}/gaps_filter_{gaps_threshold}/recall.tsv"
    params:
         number_of_points_in_ROC_curve = number_of_points_in_ROC_curve,
         tool = "snippy_all_curves"
    threads: 1
    resources:
        mem_mb = lambda wildcards, attempt: {1: 40000, 2: 60000, 3: 80000}.get(attempt, 120000)
    log:
        "logs/calculate_recall/{coverage}/snippy_all_curves/coverage_filter_{coverage_threshold}/strand_bias_filter_{strand_bias_threshold}/gaps_filter_{gaps_threshold}/recall.log"
    script:
        "../scripts/calculate_recall.py"


ruleorder: calculate_recall_for_all_snippy > calculate_recall