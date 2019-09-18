from pathlib import Path
import sys

sys.path.append(str(Path().absolute()))

from evaluate.bwa import BWA

BWA.map_query_to_ref(
    Path(snakemake.input.variant_call_probeset),
    Path(snakemake.input.reference_assembly),
    Path(snakemake.output.variant_call_probeset_mapped_to_ref),
    snakemake.threads,
)
