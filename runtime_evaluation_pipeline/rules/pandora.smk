rule copy_PRG_to_output_folder:
    input:
        PRG_original_path = PRG_folder / "{PRG_name}"
    output:
        PRG_output_path = "analysis/PRGs/{PRG_name}"
    threads: 1
    resources:
        mem_mb = lambda wildcards, attempt: 1000 * attempt
    log:
        "analysis/logs/copy_PRG_to_output_folder/{PRG_name}.log"
    shell: "cp {input.PRG_original_path} {output.PRG_output_path} 2>{log}"


rule index_PRG:
    input:
         PRG = rules.copy_PRG_to_output_folder.output.PRG_output_path,
    output:
         touch("analysis/index/{PRG_name}---threads_{threads}.index_done_flag")
    resources:
        mem_mb = lambda wildcards, attempt: 1000 * attempt
    threads:
        lambda wildcards: int(wildcards.threads)
    log:
        "analysis/logs/index_PRG/{PRG_name}---threads_{threads}.log"
    singularity:
        "shub://rmcolq/pandora:pandora"
    benchmark:
        repeat("analysis/benchmarks/index_PRG/{PRG_name}---threads_{threads}.txt", benchmark_repeat_times)
    shell: f"pandora index -w {window_size} -k {kmer_size} -t {{wildcards.threads}} {{input.PRG}} 2>{{log}}"
