rule bwa_index:
    input:
        fasta = "{fasta}"
    output:
        indexed_fasta = "{fasta}.amb"
    threads: 1
    log: "{fasta}.bwa_index.log"
    resources:
        mem_mb = lambda wildcards, attempt: 2000 * attempt
    shell: "bwa index {input.fasta} > {log} 2>&1"


rule fix_pandora_vcf_for_pipeline:
    input:
         pandora_original_vcf = "{pandora_results_dir}/{technology}/{coverage}/{subsampling}/compare_{mode}_{genotyping_mode}_genotyping/pandora_multisample_genotyped_{genotyping_mode}.vcf"
    output:
         pandora_vcf_corrected = "{pandora_results_dir}/{technology}/{coverage}/{subsampling}/compare_{mode}_{genotyping_mode}_genotyping/pandora_multisample_genotyped_{genotyping_mode}.vcf.~~vcf~~fixed~~.vcf"
    threads: 1
    resources:
        mem_mb = lambda wildcards, attempt: 2000 * attempt
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
        mem_mb = lambda wildcards, attempt: 2000 * attempt
    log:
        "logs/fix_snippy_vcf_for_pipeline{directory}/snippy_{sample}_AND_{ref}.vcf.log"
    script:
        "../scripts/fix_snippy_vcf.py"


rule make_empty_depth_file:
    input:
        file = "{file}"
    output:
        empty_depth_file = "{file}.depth"
    threads: 1
    resources:
        mem_mb = lambda wildcards, attempt: 100 * attempt
    log:
        "logs/make_empty_depth_file{file}.log"
    shell:
        "touch {output.empty_depth_file}"


rule gzip_vcf_file:
    input:
        vcf_file = "{filename}.vcf"
    output:
        gzipped_vcf_file = "{filename}.vcf.gz"
    threads: 1
    resources:
        mem_mb = lambda wildcards, attempt: 2000 * attempt
    log:
        "logs/gzip_vcf_file{filename}.log"
    shell:
        "bgzip -c {input.vcf_file} > {output.gzipped_vcf_file}"


rule index_gzipped_vcf_file:
    input:
        gzipped_vcf_file = rules.gzip_vcf_file.output.gzipped_vcf_file
    output:
        indexed_gzipped_vcf_file = "{filename}.vcf.gz.tbi"
    threads: 1
    resources:
        mem_mb = lambda wildcards, attempt: 2000 * attempt
    log:
        "logs/index_gzipped_vcf_file{filename}.log"
    shell:
        "tabix -p vcf {input.gzipped_vcf_file}"
