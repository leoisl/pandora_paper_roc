from pathlib import Path
import sys

sys.path.append(str(Path().absolute()))
from typing import Dict
from evaluate.query import Query

sample_id = snakemake.wildcards.sample_id
query_vcf = Query(
    snakemake.input.vcf,
    snakemake.input.vcf_ref,
    samples=[sample_id],
    flank_width=snakemake.params.flank_length,
)
vcf_probes: Dict[str, str] = query_vcf.make_probes()
sample_vcf_probes = vcf_probes[sample_id]
output = Path(snakemake.output.probeset)
output.write_text(sample_vcf_probes)
