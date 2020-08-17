# configs
references_folder = "references_subset"
suffix = ".fna.gz"

from pathlib import Path
from glob import glob
import subprocess
import pandas as pd

def run_command(command):
    subprocess.check_call(command, shell=True)

def get_reference_id(uncompressed_file):
    with open(uncompressed_file) as uncompressed_filehandler:
        line = uncompressed_filehandler.readline()
        reference_id = line.split()[0][1:]
        return reference_id


def main():
    reference_ids = []
    compressed_files = []
    uncompressed_files = []
    for file in glob(f"{references_folder}/*{suffix}"):
        compressed_file = Path(file).absolute()
        uncompressed_file = compressed_file.with_suffix("")
        run_command(f"gunzip -c {compressed_file} > {uncompressed_file}")

        reference_id = get_reference_id(uncompressed_file)

        reference_ids.append(reference_id)
        compressed_files.append(compressed_file)
        uncompressed_files.append(uncompressed_file)

    df = pd.DataFrame(data={"reference_id": reference_ids,
                            "compressed_file": compressed_files,
                            "uncompressed_file": uncompressed_files})
    df.to_csv("references.csv", index=False)


main()