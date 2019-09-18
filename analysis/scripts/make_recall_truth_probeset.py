from pathlib import Path
import sys

sys.path.append(str(Path().absolute()))
import logging

log_level = "DEBUG"
logging.basicConfig(
    filename=str(snakemake.log),
    filemode="w",
    level=log_level,
    format="[%(asctime)s]:%(levelname)s: %(message)s",
    datefmt="%d/%m/%Y %I:%M:%S %p",
)


from pathlib import Path
from io import StringIO
from evaluate.mummer import ShowSnps, Nucmer, DeltaFilter, NucmerError
from evaluate.utils import strip_extensions


def generate_mummer_snps(
    reference: Path,
    query: Path,
    prefix: Path = Path("out"),
    flank_width: int = 0,
    indels: bool = True,
    print_header: bool = True,
) -> StringIO:
    nucmer_params = "--maxmatch"
    nucmer = Nucmer(reference, query, str(prefix), extra_params=nucmer_params)
    nucmer_result = nucmer.run()
    if nucmer_result.returncode != 0:
        raise NucmerError(nucmer_result.stderr.decode())

    deltafile = Path(str(prefix) + ".delta")
    deltafilter_params = "-1"
    deltafilter = DeltaFilter(deltafile, extra_params=deltafilter_params)
    deltafilter_result = deltafilter.run()
    deltafilter_result.check_returncode()

    filtered_deltafile = prefix.with_suffix(".delta1")
    _ = filtered_deltafile.write_text(deltafilter_result.stdout.decode())

    showsnps_params = "-rlTC"
    showsnps = ShowSnps(
        filtered_deltafile,
        context=flank_width,
        extra_params=showsnps_params,
        indels=indels,
        print_header=print_header,
    )
    showsnps_result = showsnps.run()
    showsnps_result.check_returncode()
    showsnps_content = showsnps_result.stdout.decode()

    snpsfile = prefix.with_suffix(".snps")
    _ = snpsfile.write_text(showsnps_content)

    return StringIO(showsnps_content)


def write_truth_probeset_to_temp_file(
    truth_probes: str, query: Path, temp_folder: Path
) -> Path:
    query_name: str = strip_extensions(query).name
    truth_probes_temp_filepath: Path = temp_folder / f"{query_name}.truth_probes.fa"
    truth_probes_temp_filepath.write_text(truth_probes)

    return truth_probes_temp_filepath


# ==================================================================================
# MAIN
# ==================================================================================
query1: Path = Path(snakemake.input.truth1)
query2: Path = Path(snakemake.input.truth2)
query1_name: str = snakemake.wildcards.sample1
query2_name: str = snakemake.wildcards.sample2
prefix: Path = Path(f"{query1_name}_and_{query2_name}")
flank_width: int = snakemake.params.flank_length
mummer_snps: StringIO = generate_mummer_snps(
    reference=query1, query=query2, prefix=prefix, flank_width=flank_width
)
logging.info("Converting show-snps output to dataframe")
snps_df = ShowSnps.to_dataframe(mummer_snps)
logging.info("Creating probes from dataframe")
query1_truth_probes, query2_truth_probes = snps_df.get_probes()

query1_truth_probes_path: Path = Path(snakemake.output.probeset1)
query2_truth_probes_path: Path = Path(snakemake.output.probeset2)
logging.info("Writing output files")
query1_truth_probes_path.write_text(query1_truth_probes)
query2_truth_probes_path.write_text(query2_truth_probes)
