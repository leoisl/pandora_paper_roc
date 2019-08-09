from evaluate.cli import cli
from evaluate.mummer import Nucmer, DeltaFilter, ShowSnps
from evaluate.utils import strip_extensions
from io import StringIO
from pathlib import Path
import logging


def generate_mummer_snps(
    reference: Path, query: Path, prefix: Path = Path("out"), flank_width: int = 0
) -> StringIO:
    logging.info("Generating MUMmer SNPs file.")

    nucmer_params = "--maxmatch"
    nucmer = Nucmer(reference, query, str(prefix), extra_params=nucmer_params)
    nucmer_result = nucmer.run()
    nucmer_result.check_returncode()

    deltafile = prefix.with_suffix(".delta")
    deltafilter_params = "-1"
    deltafilter = DeltaFilter(deltafile, extra_params=deltafilter_params)
    deltafilter_result = deltafilter.run()
    deltafilter_result.check_returncode()

    filtered_deltafile = prefix.with_suffix(".delta1")
    _ = filtered_deltafile.write_text(deltafilter_result.stdout.decode())

    showsnps_params = "-rlTC"
    showsnps = ShowSnps(
        filtered_deltafile, context=flank_width, extra_params=showsnps_params
    )
    showsnps_result = showsnps.run()
    showsnps_result.check_returncode()
    showsnps_content = showsnps_result.stdout.decode()

    snpsfile = prefix.with_suffix(".snps")
    _ = snpsfile.write_text(showsnps_content)

    logging.info("Finished generating MUMmer SNPs file.")

    return StringIO(showsnps_content)


def main():
    args = cli()
    reference: Path = args.query1
    reference_name: str = strip_extensions(reference).name
    query: Path = args.query2
    query_name: str = strip_extensions(query).name
    prefix: Path = args.temp / f"{reference_name}_{query_name}"
    mummer_snps: StringIO = generate_mummer_snps(
        reference, query, prefix, args.truth_flank
    )
    snps_df = ShowSnps.to_dataframe(mummer_snps)


if __name__ == "__main__":
    main()
