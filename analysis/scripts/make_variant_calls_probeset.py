from pathlib import Path
import sys

sys.path.append(str(Path().absolute()))
from typing import Dict
from evaluate.query import Query
from evaluate.filtered_vcf_file import FilteredVCFFile
from evaluate.vcf_filters import VCF_Filters

# setup
sample_id = snakemake.wildcards.sample_id
vcf_filepath = snakemake.input.vcf
coverage_threshold = float(snakemake.wildcards.coverage_threshold)
strand_bias_threshold = float(snakemake.wildcards.strand_bias_threshold)
gaps_threshold = float(snakemake.wildcards.gaps_threshold)
vcf_ref = snakemake.input.vcf_ref
flank_width = snakemake.params.flank_length
output = Path(snakemake.output.probeset)


# API usage
filters = VCF_Filters.get_all_VCF_Filters(
    coverage_threshold=coverage_threshold,
    strand_bias_threshold=strand_bias_threshold,
    gaps_threshold=gaps_threshold,
)
filtered_vcf_file = FilteredVCFFile(vcf_filepath=vcf_filepath, filters=filters)

query_vcf = Query(
    filtered_vcf_file,
    vcf_ref,
    samples=[sample_id],
    flank_width=flank_width,
)
vcf_probes: Dict[str, str] = query_vcf.make_probes()
sample_vcf_probes = vcf_probes[sample_id]


# output
output.write_text(sample_vcf_probes)
