# configs
# pandora configs
pandora_technologies = ["illumina", "nanopore"]
pandora_modes = ["nodenovo", "withdenovo"]
pandora_genotyping = ["global"]  # "local" can also be included
pandora_base_path="../pandora_workflow/pandora_output_pandora_paper_tag1"

# other tools configs
tools=["snippy", "samtools", "medaka", "nanopolish"]
# tools=["snippy", "samtools", "medaka"]  # no nanopolish
tool_to_technology = {
    "snippy": "illumina",
    "samtools": "illumina",
    "medaka": "nanopore",
    "nanopolish": "nanopore",
}
references=['NC_011993.1', 'NC_022648.1', 'CP010121.1', 'NZ_LM995446.1', 'NZ_LT632320.1', 'NZ_CP018109.1', 'NC_017646.1', 'CP018206.1', 'NZ_CP011134.1', 'NZ_CP008697.1', 'CP010230.1', 'NC_010498.1', 'NZ_CP015228.1', 'NZ_CP009859.1', 'NZ_HG941718.1', 'CU928163.2', 'NC_011742.1', 'CP010171.1', 'NC_007779.1', 'NZ_CP016007.1', 'NZ_CP013483.1', 'CP010116.1', 'CP010226.1', 'CP010170.1']
other_tools_base_path="../variant_callers_pipeline/output_other_variant_callers_pandora_paper_tag1"

# common configs
coverages=["100x"]
subsamplings=["random"]
samples=['063_STEC', 'CFT073', 'Escherichia_coli_MINF_1A', 'Escherichia_coli_MINF_1D', 'Escherichia_coli_MINF_7C', 'Escherichia_coli_MINF_8D', 'Escherichia_coli_MINF_9A', 'Escherichia_coli_MSB1_1A', 'Escherichia_coli_MSB1_3B', 'Escherichia_coli_MSB1_4E', 'Escherichia_coli_MSB1_4I', 'Escherichia_coli_MSB1_6C', 'Escherichia_coli_MSB1_7A', 'Escherichia_coli_MSB1_7C', 'Escherichia_coli_MSB1_8B', 'Escherichia_coli_MSB1_8G', 'Escherichia_coli_MSB1_9D', 'Escherichia_coli_MSB2_1A', 'H131800734', 'ST38']
output_file="variant_calls_pandora_paper_tag1.csv"
# output_file="variant_calls_pandora_paper_tag1.no_nanopolish.csv"  # no nanopolish




# processing
from itertools import product
with open(output_file, "w") as output_filehandler:
    print("sample_id,tool,coverage,reference,vcf", file=output_filehandler)
    for sample, coverage, technology, subsampling, mode, genotyping in product(
            samples, coverages, pandora_technologies, subsamplings, pandora_modes, pandora_genotyping):
        print(",".join([
            sample,
            f"pandora_{technology}_{mode}_{genotyping}_genotyping",
            coverage,
            f"{pandora_base_path}/{technology}/{coverage}/{subsampling}/compare_{mode}_{genotyping}_genotyping/pandora_multisample.vcf_ref.fa",
            f"{pandora_base_path}/{technology}/{coverage}/{subsampling}/compare_{mode}_{genotyping}_genotyping/pandora_multisample_genotyped_{genotyping}.vcf"]),
        file=output_filehandler)

    for tool, sample, reference, coverage, subsampling in product(
            tools, samples, references, coverages, subsamplings):
        technology=tool_to_technology[tool]
        print(",".join([
            sample,
            f"{tool}_{reference}",
            coverage,
            f"{other_tools_base_path}/{tool}/{technology}/{coverage}/{subsampling}/{sample}/{tool}_{sample}_AND_{reference}.ref.fa",
            f"{other_tools_base_path}/{tool}/{technology}/{coverage}/{subsampling}/{sample}/{tool}_{sample}_AND_{reference}.vcf"]),
        file=output_filehandler)
