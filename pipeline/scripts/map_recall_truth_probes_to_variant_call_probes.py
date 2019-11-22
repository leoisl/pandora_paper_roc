from pathlib import Path
import sys
sys.path.append(str(Path().absolute()))
import logging
log_level = "INFO"
logging.basicConfig(
    filename=str(snakemake.log),
    filemode="w",
    level=log_level,
    format="[%(asctime)s]:%(levelname)s: %(message)s",
    datefmt="%d/%m/%Y %I:%M:%S %p",
)



from evaluate.bwa import BWA


# setup
query = Path(snakemake.input.truth_probeset)
ref = Path(snakemake.input.variant_calls_probeset)
output = Path(snakemake.output.sam)
threads = int(snakemake.threads)


# API usage
logging.info(f"Mapping {query} to {ref}")
BWA.map_query_to_ref(
    query=query,
    ref=ref,
    output=output,
    threads=threads,
)


logging.info(f"Done")