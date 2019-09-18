rule make_variant_calls_probeset:
    input:
         vcf = lambda wildcards: data.xs((wildcards.sample_id, wildcards.coverage, wildcards.tool))["vcf"],
         vcf_ref = lambda wildcards: data.xs((wildcards.sample_id, wildcards.coverage, wildcards.tool))["vcf_reference"]
    output:
          probeset = "analysis/truth_probesets/{sample_id}/{coverage}/{tool}.truth_probeset.fa"
    params:
          flank_length = config["truth_calls_probe_length"]
    threads: 1
    resources:
        mem_mb = lambda wildcards, attempt: 1000 * attempt
    log:
        "logs/make_variant_calls_probeset_{sample_id}_{coverage}_{tool}.log"
    script:
        "../scripts/make_variant_calls_probeset.py"


