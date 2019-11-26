import pysam


length_filter = 500

long_enough_genes = set()
with pysam.FastxFile("pandora_multisample.vcf_ref.fa") as genes:
    for gene in genes:
        gene_is_long_enough = len(gene.sequence) > length_filter
        if gene_is_long_enough:
            long_enough_genes.add(gene.name)

# with pysam.VariantFile(pysam.tabix_index("pandora_multisample_genotyped.vcf", preset="vcf", keep_original=True, force=True)) as input_vcf:
#     with pysam.VariantFile(f"pandora_multisample_genotyped.length_filtered_{length_filter}.vcf", "w", template=input_vcf) as filtered_vcf:
#         for gene in long_enough_genes:
#             try:
#                 for record in input_vcf.fetch(contig=gene):
#                     filtered_vcf.write(record)
#             except ValueError:
#                 pass


with open("pandora_multisample_genotyped.vcf") as input_vcf:
    with open(f"pandora_multisample_genotyped.length_filtered_{length_filter}.vcf", "w") as filtered_vcf:
        for line in input_vcf:
            line_split = line.split()
            gene = line_split[0]
            if line.startswith("#"):
                filtered_vcf.write(line)
            elif gene in long_enough_genes:
                filtered_vcf.write(line)
            else:
                print(f"Filtering out a record from {gene}")
