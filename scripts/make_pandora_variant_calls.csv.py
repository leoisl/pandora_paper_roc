from itertools import product
samples = ["CFT073", "H131800734", "ST38", "063_STEC"]
coverages = ["30x", "60x", "100x"]
technologies = ["illumina", "nanopore"]
modes = ["no_denovo", "with_denovo"]
base_path = "/hps/nobackup/research/zi/projects/pandora_paper_leandro/analysis"

for sample, coverage, technology, mode in product(samples, coverages, technologies, modes):
    print(",".join([sample, f"pandora_{technology}_{mode}", coverage, f"{base_path}/{technology}/{coverage}/random/compare_{mode}/pandora_multisample_genotyped.vcf"]))
