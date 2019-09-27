rule copy_PRG_to_output_folder:
    input:
        PRG_original_path = PRG_folder / "{PRG_name}"
    output:
        PRG_output_path = "analysis/PRGs/{PRG_name}"
    threads: 1
    resources:
        mem_mb = lambda wildcards, attempt: 1000 * attempt
    log:
        "logs/copy_PRG_to_output_folder/{PRG_name}.log"
    shell: "cp {input.PRG_original_path} {output.PRG_output_path}"

