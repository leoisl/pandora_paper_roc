rule bwa_index:
    input:
        fasta = "{fasta}"
    output:
        indexed_fasta = "{fasta}.amb"
    threads: 1
    log: "{fasta}.bwa_index.log"
    resources:
        mem_mb = lambda wildcards, attempt: 2000 * attempt
    singularity:
        "docker://leandroishilima/pandora1_paper_basic_tools:pandora_paper_tag1"
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
        "../scripts/fix_pandora_vcf_before_GCP.py"


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


rule fix_samtools_vcf_for_pipeline:
    input:
         samtools_original_vcf =  "{directory}/samtools_{sample}_AND_{ref}.vcf"
    output:
         samtools_vcf_corrected = "{directory}/samtools_{sample}_AND_{ref}.vcf.~~vcf~~fixed~~.vcf"
    threads: 1
    resources:
        mem_mb = lambda wildcards, attempt: 2000 * attempt
    log:
        "logs/fix_samtools_vcf_for_pipeline{directory}/samtools_{sample}_AND_{ref}.vcf.log"
    script:
        "../scripts/fix_samtools_vcf.py"

rule fix_medaka_vcf_for_pipeline:
    input:
         medaka_original_vcf =  "{directory}/medaka_{sample}_AND_{ref}.vcf"
    output:
         medaka_vcf_corrected = "{directory}/medaka_{sample}_AND_{ref}.vcf.~~vcf~~fixed~~.vcf"
    threads: 1
    resources:
        mem_mb = lambda wildcards, attempt: 2000 * attempt
    log:
        "logs/fix_medaka_vcf_for_pipeline{directory}/medaka_{sample}_AND_{ref}.vcf.log"
    script:
        "../scripts/fix_medaka_vcf.py"



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
    singularity:
        "docker://leandroishilima/pandora1_paper_basic_tools:pandora_paper_tag1"
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
    singularity:
        "docker://leandroishilima/pandora1_paper_basic_tools:pandora_paper_tag1"
    shell:
        "tabix -p vcf {input.gzipped_vcf_file}"
