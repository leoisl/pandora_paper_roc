import argparse
import pandas as pd
from pathlib import Path
import subprocess
import os

def run_command(command):
    print(f"Running {command}")
    subprocess.check_call(command, shell=True)

def get_args():
    parser = argparse.ArgumentParser(description="Run snippy using Rachel's nextflow pipeline.")
    parser.add_argument('--reference_directory', type=str, help='Directory of different references to use', required=True)
    parser.add_argument('--samples_csv', type=str, help='CSV with two columns. 1st is the sample id. 2nd is the dir with sample reads.', required=True)
    parser.add_argument('--coverage', type=str, help='Subsampling coverage', required=True)
    parser.add_argument('--output_dir', type=str, help='Output dir', required=True)
    args = parser.parse_args()
    return args

def main():
    args = get_args()
    reference_directory = args.reference_directory
    samples = pd.read_csv(args.samples_csv)
    output_dir = Path(args.output_dir)
    coverage = args.coverage
    script_dir = Path(os.path.dirname(os.path.realpath(__file__))).absolute()

    for sample_id, sample_path in zip(samples["sample_id"], samples["sample_path"]):
        print(f"Processing {sample_id}...")
        snippy_output_dir = output_dir / sample_id
        snippy_output_dir.mkdir(parents=True, exist_ok=True)
        run_command(f"cd {snippy_output_dir} && "
                    f"bsub.py 2 run_genotype_snippy_{sample_id} "
                    f"nextflow run {script_dir / 'genotype_snippy.nf'} -w work -c {script_dir / 'nextflow.config'} "
                    f"--reference_directory {reference_directory} --reads_directory {sample_path} "
                    f"--sample_id {sample_path} --coverage {coverage} -resume")


if __name__ == "__main__":
    main()