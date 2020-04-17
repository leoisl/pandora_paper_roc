from itertools import product
samples = ["063_STEC", "CFT073", "Escherichia_coli_MINF_1A", "Escherichia_coli_MINF_1D", "Escherichia_coli_MINF_2E", "Escherichia_coli_MINF_7C", "Escherichia_coli_MINF_8D", "Escherichia_coli_MINF_9A", "Escherichia_coli_MSB1_1A", "Escherichia_coli_MSB1_3B", "Escherichia_coli_MSB1_3I", "Escherichia_coli_MSB1_4E", "Escherichia_coli_MSB1_4I", "Escherichia_coli_MSB1_6C", "Escherichia_coli_MSB1_6J", "Escherichia_coli_MSB1_7A", "Escherichia_coli_MSB1_7C", "Escherichia_coli_MSB1_8B", "Escherichia_coli_MSB1_8G", "Escherichia_coli_MSB1_9D", "Escherichia_coli_MSB1_9I", "Escherichia_coli_MSB2_1A", "H131800734", "ST38"]
coverages = ["100x"]
technologies = ["illumina"]
modes = ["nodenovo", "withdenovo"]
genotyping_modes = ["global"]
base_path = "/hps/nobackup/iqbal/leandro/pandora_analysis_pipeline/analysis_24_way"

print("sample_id,tool,coverage,reference,vcf")
for sample, coverage, technology, mode, genotyping_mode in product(samples, coverages, technologies, modes, genotyping_modes):
    print(",".join([sample,
                    f"pandora_{technology}_{mode}_{genotyping_mode}_genotyping",
                    coverage,
                    f"{base_path}/{technology}/{coverage}/random/compare_{mode}_{genotyping_mode}_genotyping/pandora_multisample.vcf_ref.fa",
                    f"{base_path}/{technology}/{coverage}/random/compare_{mode}_{genotyping_mode}_genotyping/pandora_multisample_genotyped_{genotyping_mode}.vcf"]))
