rule make_pandora_probeset:
    input:
         vcf = lambda wildcards: data.query(f"sample_id == '{wildcards.sample_id}'")['vcf'][0],
         vcf_ref = lambda wildcards: data.query(f"sample_id == '{wildcards.sample_id}'")['vcf_reference'][0]
    output:
          probeset = "analysis/truth_probesets/{sample_id}/{coverage}/{tool}.truth_probeset.fa"
    params:
          flank_length = config["truth_calls_probe_length"]
    threads: 1
    resources:
        mem_mb = lambda wildcards, attempt: config["mem_mb"][attempt-1]
    log:
        "logs/make_pandora_probeset_{sample_id}_{coverage}_{tool}.log"
    # singularity:
    #     config["singularity_image"]
    script:
        "../scripts/make_pandora_probeset.py"


