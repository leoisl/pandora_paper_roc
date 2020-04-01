import argparse
import pandas as pd
from pathlib import Path
import subprocess
import os

def run_command(command, dry_run):
    print(f"Running {command}")

    if not dry_run:
        subprocess.check_call(command, shell=True)

def get_args():
    parser = argparse.ArgumentParser(description="Run snippy using Rachel's nextflow pipeline.")
    parser.add_argument('--reference_directory', type=str, help='Directory of different references to use', required=True)
    parser.add_argument('--samples_csv', type=str, help='CSV with two columns. 1st is the sample id. 2nd is the dir with sample reads.', required=True)
    parser.add_argument('--coverage', type=str, help='Subsampling coverage', required=True)
    parser.add_argument('--output_dir', type=str, help='Output dir', required=True)
    parser.add_argument('--dry_run', action="store_true", help='Do not run, but just outputs the command lines', required=False, default=False)
    args = parser.parse_args()
    return args

def main():
    args = get_args()
    reference_directory = Path(args.reference_directory).absolute()
    samples = pd.read_csv(args.samples_csv)
    output_dir = Path(args.output_dir).absolute()
    coverage = args.coverage
    script_dir = Path(os.path.dirname(os.path.realpath(__file__))).absolute()
    dry_run = args.dry_run

    for sample_id, sample_path in zip(samples["sample_id"], samples["sample_path"]):
        print(f"Processing {sample_id}...")
        sample_path = Path(sample_path).absolute()
        snippy_output_dir = output_dir / sample_id
        snippy_output_dir.mkdir(parents=True, exist_ok=True)
        run_command(f"cd {snippy_output_dir} && "
                    f"bsub.py 2 run_genotype_snippy_{sample_id} "
                    f"nextflow run {script_dir / 'genotype_snippy.nf'} -w work -c {script_dir / 'nextflow.config'} "
                    f"--reference_directory {reference_directory} --reads_directory {sample_path} "
                    f"--sample_id {sample_id} --coverage {coverage} -resume",
                    dry_run)


if __name__ == "__main__":
    main()