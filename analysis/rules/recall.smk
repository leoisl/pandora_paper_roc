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
    # singularity:
    #     config["singularity_image"]
    script:
        "../scripts/make_recall_truth_probeset.py"


rule map_recall_truth_probes_to_variant_call_probes:
    input:
         truth_probeset = "analysis/recall/truth_probesets/{sample_id}/{filename_prefix}.truth_probeset.fa",
         variant_calls_probeset = "analysis/variant_calls_probesets/{sample_id}/{coverage}/{tool}.variant_calls_probeset.fa",
    output:
         sam = "analysis/recall/map_probes/{sample_id}/{coverage}/{tool}/{filename_prefix}.sam"
    threads: 2
    resources:
        mem_mb = lambda wildcards, attempt: 2000 * attempt
    log:
        "logs/map_recall_truth_probes_to_variant_call_probes/{sample_id}/{coverage}/{tool}/{filename_prefix}.log"
    # singularity:
    #     config["singularity_image"]
    script:
        "../scripts/map_recall_truth_probes_to_variant_call_probes.py"

