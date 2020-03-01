rule bwa_index:
    input:
        fasta = "{fasta}"
    output:
        indexed_fasta = "{fasta}.amb"
    threads: 1
    log: "{fasta}.bwa_index.log"
    resources:
        mem_mb = lambda wildcards, attempt: 4000 * attempt
    shell: "bwa index {input.fasta} > {log} 2>&1"


rule fix_pandora_sample_name:
    input:
         pandora_original_vcf = "{pandora_results_dir}/{technology}/{coverage}/{subsampling}/compare_{mode}/pandora_multisample_genotyped.vcf"
    output:
         pandora_vcf_with_sample_names_corrected = "{pandora_results_dir}/{technology}/{coverage}/{subsampling}/compare_{mode}/pandora_multisample_genotyped.vcf.~~sample~~names~~fixed.vcf"
    threads: 1
    resources:
        mem_mb = lambda wildcards, attempt: 4000 * attempt
    log:
        "logs/fix_pandora_sample_name{pandora_results_dir}/{technology}/{coverage}/{subsampling}/compare_{mode}/pandora_multisample_genotyped.vcf.log"
    script:
        "../scripts/fix_pandora_sample_name.py"


rule fix_snippy_sample_name:
    # TODO: for now we are assuming all snippy vcf were corrected with scripts/make_snippy_variant_calls.csv.py
    input:
         snippy_original_vcf = "{filename}.vcf"
    output:
         snippy_vcf_with_sample_names_corrected = "{filename}.vcf.~~sample~~names~~fixed.vcf"
    wildcard_constraints:
        filename=".*/snippy_[^/]+"
    threads: 1
    resources:
        mem_mb = lambda wildcards, attempt: 4000 * attempt
    log:
        "logs/fix_snippy_sample_name{filename}.log"
    shell:
        "cp {input.snippy_original_vcf} {output.snippy_vcf_with_sample_names_corrected}"