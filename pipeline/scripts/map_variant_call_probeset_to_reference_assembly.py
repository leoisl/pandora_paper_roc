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
query = Path(snakemake.input.variant_call_probeset)
ref = Path(snakemake.input.reference_assembly)
output = Path(snakemake.output.variant_call_probeset_mapped_to_ref)
threads = snakemake.threads



# API usage
logging.info(f"Mapping {query} to {ref}")
BWA.map_query_to_ref(
    query=query,
    ref=ref,
    output=output,
    threads=threads,
)
logging.info(f"Done")