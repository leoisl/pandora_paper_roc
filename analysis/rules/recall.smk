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

#
# rule map_recall_truth_probes_to_pandora_probes:
#     input:
#          truth_probeset_fasta_file = rules.make_recall_truth_probeset.output.truth_probeset_fasta_file,
#          pandora_probes_from_sample_1 = output_dir_path / "{sample_1}.pandora_probeset",
#          pandora_probes_from_sample_2 = output_dir_path / "{sample_2}.pandora_probeset"
#     output:
#          truth_probes_mapped_to_pandora_probes_from_sample_1 = rules.make_recall_truth_probeset.output.truth_probeset_fasta_file + ".mapped_to_pandora_probes_{sample_1}",
#          truth_probes_mapped_to_pandora_probes_from_sample_2 = rules.make_recall_truth_probeset.output.truth_probeset_fasta_file + ".mapped_to_pandora_probes_{sample_2}"
#     threads: 1
#     resources:
#         mem_mb = lambda wildcards, attempt: config["mem_mb"][attempt-1]
#     log:
#         "logs/map_recall_truth_probes_to_pandora_probes_{sample_1}---{sample_2}.log"
#     singularity:
#         config["singularity_image"]
#     script:
#         "../scripts/map_recall_truth_probes_to_pandora_probes.py"
#
