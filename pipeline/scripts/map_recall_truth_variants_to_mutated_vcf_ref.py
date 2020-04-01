from pathlib import Path
import sys
sys.path.append(str(Path().absolute()))
import logging
import subprocess
import time


log_level = "INFO"
logging.basicConfig(
    filename=str(snakemake.log),
    filemode="w",
    level=log_level,
    format="[%(asctime)s]:%(levelname)s: %(message)s",
    datefmt="%d/%m/%Y %I:%M:%S %p",
)
from evaluate.bwa import BWA


def run_bwa_mem(query, ref):
    nb_tries = 10
    success = False
    for _ in range(nb_tries):
        try:
            logging.info(f"Mapping {query} to {ref}")
            logging.info(f"bwa index {ref}")
            subprocess.check_call(f"bwa index {ref}", shell=True)
            time.sleep(60.0)  # sleeps 60 seconds due to cluster file latency

            logging.info(f"bwa mem {ref} {query}")
            BWA.map_query_to_ref(
                query=query,
                ref=ref,
                output=output,
                threads=threads,
            )

            # all good
            success = True
            break
        except subprocess.CalledProcessError as error:
            last_error = error
            logging.info("CalledProcessError captured, retrying...")

        if not success:
            raise last_error


# setup
query = Path(snakemake.input.truth_probeset)
refs = [Path(ref) for ref in snakemake.input.mutated_vcf_refs]
outputs = [Path(sam) for sam in snakemake.output.sams]
threads = int(snakemake.threads)


# API usage
for ref, output in zip(refs, outputs):
    run_bwa_mem(query, ref)

logging.info(f"Done")