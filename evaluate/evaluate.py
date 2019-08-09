from evaluate.cli import cli
from evaluate.mummer import Nucmer
from io import StringIO
from pathlib import Path

def generate_mummer_snps(reference: Path, query: Path) -> StringIO:
    nucmer_params = "--maxmatch"
    nucmer = Nucmer(reference, query, prefix, extra_params=nucmer_params)

def main():
    args = cli()
    reference: Path = args.query1
    query: Path = args.query2
    prefix: Path = args.temp / reference.name
    mummer_snps = generate_mummer_snps(reference, query)



if __name__ == "__main__":
    main()
