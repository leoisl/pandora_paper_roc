rule index_vcf:
    input:
        vcf = "{vcf}"
    output:
        indexed_vcf = "{vcf}.gz.tbi"
    threads: 1
    resources:
        mem_mb = lambda wildcards, attempt: 1000 * attempt
    run:
        import pysam
        pysam.tabix_index(input.vcf, preset="vcf", keep_original=True)


rule bwa_index:
    input:
        fasta = "{fasta}"
    output:
        indexed_fasta = "{fasta}.amb"
    threads: 1
    resources:
        mem_mb = lambda wildcards, attempt: 1000 * attempt
    shell: "bwa index {input.fasta}"



rule make_variant_calls_probeset:
    input:
         vcf = lambda wildcards: data.xs((wildcards.sample_id, wildcards.coverage, wildcards.tool))["vcf"],
         vcf_index = lambda wildcards: data.xs((wildcards.sample_id, wildcards.coverage, wildcards.tool))["vcf"]+".gz.tbi",
         vcf_ref = lambda wildcards: data.xs((wildcards.sample_id, wildcards.coverage, wildcards.tool))["vcf_reference"]
    output:
          probeset = "analysis/variant_calls_probesets/{sample_id}/{coverage}/{tool}/coverage_filter_{coverage_threshold}/strand_bias_filter_{strand_bias_threshold}/gaps_filter_{gaps_threshold}/variant_calls_probeset.fa"
    params:
          flank_length = config["variant_calls_probe_length"]
    threads: 1
    resources:
        mem_mb = lambda wildcards, attempt: 1000 * attempt
    log:
        "logs/make_variant_calls_probeset/{sample_id}/{coverage}/{tool}/coverage_filter_{coverage_threshold}/strand_bias_filter_{strand_bias_threshold}/gaps_filter_{gaps_threshold}/variant_calls_probeset.log"
    script:
        "../scripts/make_variant_calls_probeset.py"
