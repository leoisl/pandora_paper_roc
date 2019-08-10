from mummer import Nucmer, DeltaFilter, ShowSnps
from cli import cli
from utils import strip_extensions
from query import Query
from io import StringIO
from pathlib import Path
import logging
import pandas as pd
from typing import Tuple


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


def make_truth_panels_from_snps_dataframe(snps_df: pd.DataFrame) -> Tuple[str, str]:

    ref_probes = ""
    query_probes = ""

    for index, row in snps_df.iterrows():
        ref_name = f">{row.ref_chrom}_POS={row.ref_pos}_SUB={row.ref_sub}"
        ref_probe = row.ref_context.replace(".", "").replace("-", "")
        ref_probes += f"{ref_name}\n{ref_probe}\n"
        query_name = f">{row.query_chrom}_POS={row.query_pos}_SUB={row.query_sub}"
        query_probe = row.query_context.replace(".", "").replace("-", "")
        query_probes += f"{query_name}\n{query_probe}\n"

    return ref_probes, query_probes


def write_vcf_probes_to_file(vcf: Path, probes: str, tempdir: Path):
    vcf_name: str = strip_extensions(vcf).name
    vcf_probes_path: Path = tempdir / f"{vcf_name}.probes.fa"
    vcf_probes_path.write_text(probes)
    logging.info(f"VCF probes written to file: {vcf_probes_path}")


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
    # todo: merge consecutive positions
    logging.info("Making truth probesets.")
    ref_truth_probes, query_truth_probes = make_truth_panels_from_snps_dataframe(
        snps_df
    )

    ref_truth_probes_path: Path = args.temp / f"{reference_name}.truth_probes.fa"
    query_truth_probes_path: Path = args.temp / f"{query_name}.truth_probes.fa"
    ref_truth_probes_path.write_text(ref_truth_probes)
    logging.info(
        f"{reference_name} truth probes written to: {str(ref_truth_probes_path)}"
    )
    query_truth_probes_path.write_text(query_truth_probes)
    logging.info(
        f"{query_name} truth probes written to: {str(query_truth_probes_path)}"
    )

    logging.info("Making probes for VCF")
    query_vcf = Query(args.vcf, args.vcf_ref, flank_width=args.query_flank)
    vcf_probes: str = query_vcf.make_probes()
    write_vcf_probes_to_file(args.vcf, vcf_probes, args.temp)


if __name__ == "__main__":
    main()
