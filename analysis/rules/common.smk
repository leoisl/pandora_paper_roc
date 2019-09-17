rule make_pandora_probeset:
    input:
         vcf = lambda wildcards: data.query(f"sample_id == '{wildcards.sample_id}'")['vcf'][0],
         vcf_ref = config["pandora_vcf_ref"]
    output:
          probeset = "analysis/truth_probesets/{sample_id}/{coverage}/{tool}.truth_probeset.fa"
    params:
          flank_length = config["pandora_probes_flank_length"]
    threads: 1
    resources:
        mem_mb = lambda wildcards, attempt: config["mem_mb"][attempt-1]
    log:
        "logs/make_pandora_probeset_{sample}.log"
    singularity:
        config["singularity_image"]
    script:
        "../scripts/make_pandora_probeset.py"


