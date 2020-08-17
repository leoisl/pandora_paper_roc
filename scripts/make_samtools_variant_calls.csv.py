# configs
technologies=["illumina"]
coverages=["100x"]
subsamplings=["random"]
samples=["063_STEC", "CFT073", "Escherichia_coli_MINF_1A", "Escherichia_coli_MINF_1D", "Escherichia_coli_MINF_7C", "Escherichia_coli_MINF_8D", "Escherichia_coli_MINF_9A", "Escherichia_coli_MSB1_1A", "Escherichia_coli_MSB1_3B", "Escherichia_coli_MSB1_4E", "Escherichia_coli_MSB1_4I", "Escherichia_coli_MSB1_6C", "Escherichia_coli_MSB1_7A", "Escherichia_coli_MSB1_7C", "Escherichia_coli_MSB1_8B", "Escherichia_coli_MSB1_8G", "Escherichia_coli_MSB1_9D", "Escherichia_coli_MSB2_1A", "H131800734", "ST38"]
references=['NC_011993.1', 'NC_022648.1', 'CP010121.1', 'NZ_LM995446.1', 'NZ_LT632320.1', 'NZ_CP018109.1', 'NC_017646.1', 'CP018206.1', 'NZ_CP011134.1', 'NZ_CP008697.1', 'NC_004431.1', 'CP010230.1', 'NC_010498.1', 'NZ_CP015228.1', 'NZ_CP009859.1', 'NZ_HG941718.1', 'CU928163.2', 'NC_011742.1', 'CP010171.1', 'NC_007779.1', 'NZ_CP016007.1', 'NZ_CP013483.1', 'CP010116.1', 'CP010226.1', 'CP010170.1']
base_path="/hps/nobackup/iqbal/leandro/samtools_analysis_pipeline"

from itertools import product

print("sample_id,tool,coverage,reference,vcf")
for sample, reference, coverage, subsampling, technology in product(samples, references, coverages, subsamplings, technologies):
    print(",".join([sample,
                    f"samtools_{reference}",
                    coverage,
                    f"{base_path}/{technology}/{coverage}/{subsampling}/{sample}/samtools_{sample}_AND_{reference}.ref.fa",
                    f"{base_path}/{technology}/{coverage}/{subsampling}/{sample}/samtools_{sample}_AND_{reference}.vcf"]))
