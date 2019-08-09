from evaluate.cli import cli
from evaluate.mummer import Nucmer
from evaluate.utils import strip_extensions
from io import StringIO
from pathlib import Path


def generate_mummer_snps(
    reference: Path, query: Path, prefix: Path = Path("out")
) -> StringIO:
    nucmer_params = "--maxmatch"
    nucmer = Nucmer(reference, query, str(prefix), extra_params=nucmer_params)


def main():
    args = cli()
    reference: Path = args.query1
    reference_name: str = strip_extensions(reference).name
    query: Path = args.query2
    query_name: str = strip_extensions(query).name
    prefix: Path = args.temp / f"{reference_name}_{query_name}"
    mummer_snps = generate_mummer_snps(reference, query, prefix)


if __name__ == "__main__":
    main()
