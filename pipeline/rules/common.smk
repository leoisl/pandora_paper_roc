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


rule fix_pandora_vcf_for_pipeline:
    input:
         pandora_original_vcf = "{pandora_results_dir}/{technology}/{coverage}/{subsampling}/compare_{mode}_{genotyping_mode}_genotyping/pandora_multisample_genotyped_{genotyping_mode}.vcf"
    output:
         pandora_vcf_corrected = "{pandora_results_dir}/{technology}/{coverage}/{subsampling}/compare_{mode}_{genotyping_mode}_genotyping/pandora_multisample_genotyped_{genotyping_mode}.vcf.~~vcf~~fixed~~.vcf"
    threads: 1
    resources:
        mem_mb = lambda wildcards, attempt: 4000 * attempt
    log:
        "logs/fix_pandora_vcf_for_pipeline{pandora_results_dir}/{technology}/{coverage}/{subsampling}/compare_{mode}_{genotyping_mode}_genotyping/pandora_multisample_genotyped_{genotyping_mode}.vcf.log"
    script:
        "../scripts/fix_pandora_vcf.py"


rule fix_snippy_vcf_for_pipeline:
    input:
         snippy_original_vcf =  "{directory}/snippy_{sample}_AND_{ref}.vcf"
    output:
         snippy_vcf_corrected = "{directory}/snippy_{sample}_AND_{ref}.vcf.~~vcf~~fixed~~.vcf"
    threads: 1
    resources:
        mem_mb = lambda wildcards, attempt: 4000 * attempt
    log:
        "logs/fix_snippy_vcf_for_pipeline{directory}/snippy_{sample}_AND_{ref}.vcf.log"
    shell:
        "../scripts/fix_snippy_vcf.py"


rule remove_plasmids_from_truth_sample:
    input:
         truth_sample = "{truth_sample}"
    output:
         truth_sample_with_chrom_only = "{truth_sample}.chrom_only"
    threads: 1
    resources:
        mem_mb = lambda wildcards, attempt: 4000 * attempt
    log:
        "logs/remove_plasmids_from_truth_sample{truth_sample}.log"
    script:
        "../scripts/remove_plasmids_from_truth_sample.py"