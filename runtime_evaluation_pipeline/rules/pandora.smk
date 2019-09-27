rule copy_PRGs_to_output_folder:
    input:
        PRGs_original_paths = PRGs_original_paths
    output:
        PRGs_output_paths = PRGs_output_paths
    threads: 1
    resources:
        mem_mb = lambda wildcards, attempt: 1000 * attempt
    log:
        "logs/copy_PRG_to_output.log"
    run:
        for PRG_original_path, PRG_output_path in zip(input.PRGs_original_paths, output.PRGs_output_paths):
            shell(f"cp {PRG_original_path} {PRG_output_path}")