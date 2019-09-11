# from evaluate.Positioned_SNPs import *
#
# length_of_longest_contig = 5572075
#
#
# def test_PositionedSNP_merge():
#     positioned_SNP_1 = PositionedSNP.from_dict({'alleles': ('C', 'T'), 'positions': {Position(genome='genome.1', chrom='NC_011353.1', pos=2372140), Position(genome='genome.2', chrom='NZ_CP009644.1', pos=1846009)}})
#     positioned_SNP_2 = PositionedSNP.from_dict({'alleles': ('C', 'T'), 'positions': {Position(genome='genome.3', chrom='NZ_CP016628.1', pos=2720609), Position(genome='genome.1', chrom='NC_011353.1', pos=2157387)}})
#
#     actual = positioned_SNP_1.merge(positioned_SNP_2)
#     expected = PositionedSNP.from_dict({'alleles': ('C', 'T'), 'positions': {Position(genome='genome.3', chrom='NZ_CP016628.1', pos=2720609), Position(genome='genome.2', chrom='NZ_CP009644.1', pos=1846009), Position(genome='genome.1', chrom='NC_011353.1', pos=2157387), Position(genome='genome.1', chrom='NC_011353.1', pos=2372140)}})
#
#     assert actual == expected
# 
#
# def test_PositionedSNPsIndex_add_PositionedSNP_no_previous_PositionedSNP():
#     positionedSNPsIndex = PositionedSNPsIndex(length_of_longest_contig)
#     position_1 = Position(genome='genome.1', chrom='NC_011350.1', pos=75224)
#     position_2 = Position(genome='genome.2', chrom='NZ_CP009644.1', pos=4166248)
#     allele_1 = "G"
#     allele_2 = "T"
#     positionedSNPsIndex.add_PositionedSNP(position_1, position_2, allele_1, allele_2)
#
#     actual = positionedSNPsIndex.to_dict()
#     expected = {75224: {('genome.1', 'NC_011350.1', 'G', 'T'): {'alleles': ('G', 'T'), 'positions': {Position(genome='genome.2', chrom='NZ_CP009644.1', pos=4166248), Position(genome='genome.1', chrom='NC_011350.1', pos=75224)}}}, 4166248: {('genome.2', 'NZ_CP009644.1', 'G', 'T'): {'alleles': ('G', 'T'), 'positions': {Position(genome='genome.2', chrom='NZ_CP009644.1', pos=4166248), Position(genome='genome.1', chrom='NC_011350.1', pos=75224)}}}}
#     assert actual == expected
#
# def test_PositionedSNPsIndex_add_PositionedSNP_only_previous_PositionedSNP_1_exists():
#     position_1 = Position(genome='genome.1', chrom='NC_011353.1', pos=3772004)
#     position_2 = Position(genome='genome.2', chrom='NZ_CP009644.1', pos=2929307)
#     position_3 = Position(genome='genome.1', chrom='NC_011353.1', pos=1743676)
#     allele_1 = "A"
#     allele_2 = "G"
#
#     positionedSNPsIndex = PositionedSNPsIndex(length_of_longest_contig)
#     positionedSNPsIndex.add_PositionedSNP(position_2, position_3, allele_1, allele_2)
#     positionedSNPsIndex.add_PositionedSNP(position_1, position_2, allele_1, allele_2)
#
#     actual = positionedSNPsIndex.to_dict()
#     expected = {1743676: {('genome.1', 'NC_011353.1', 'A', 'G'): {'alleles': ('A', 'G'), 'positions': {Position(genome='genome.1', chrom='NC_011353.1', pos=1743676), Position(genome='genome.2', chrom='NZ_CP009644.1', pos=2929307), Position(genome='genome.1', chrom='NC_011353.1', pos=3772004)}}}, 2929307: {('genome.2', 'NZ_CP009644.1', 'A', 'G'): {'alleles': ('A', 'G'), 'positions': {Position(genome='genome.1', chrom='NC_011353.1', pos=1743676), Position(genome='genome.2', chrom='NZ_CP009644.1', pos=2929307), Position(genome='genome.1', chrom='NC_011353.1', pos=3772004)}}}, 3772004: {('genome.1', 'NC_011353.1', 'A', 'G'): {'alleles': ('A', 'G'), 'positions': {Position(genome='genome.1', chrom='NC_011353.1', pos=1743676), Position(genome='genome.2', chrom='NZ_CP009644.1', pos=2929307), Position(genome='genome.1', chrom='NC_011353.1', pos=3772004)}}}}
#
#     assert actual == expected
#
# def test_PositionedSNPsIndex_add_PositionedSNP_only_previous_PositionedSNP_2_exists():
#     position_1 = Position(genome='genome.1', chrom='NC_011353.1', pos=3772004)
#     position_2 = Position(genome='genome.2', chrom='NZ_CP009644.1', pos=2929307)
#     position_3 = Position(genome='genome.1', chrom='NC_011353.1', pos=1743676)
#     allele_1 = "A"
#     allele_2 = "G"
#
#     positionedSNPsIndex = PositionedSNPsIndex(length_of_longest_contig)
#     positionedSNPsIndex.add_PositionedSNP(position_2, position_3, allele_1, allele_2)
#     positionedSNPsIndex.add_PositionedSNP(position_2, position_1, allele_1, allele_2)
#
#     actual = positionedSNPsIndex.to_dict()
#     expected = {1743676: {('genome.1', 'NC_011353.1', 'A', 'G'): {'alleles': ('A', 'G'), 'positions': {Position(genome='genome.1', chrom='NC_011353.1', pos=1743676), Position(genome='genome.2', chrom='NZ_CP009644.1', pos=2929307), Position(genome='genome.1', chrom='NC_011353.1', pos=3772004)}}}, 2929307: {('genome.2', 'NZ_CP009644.1', 'A', 'G'): {'alleles': ('A', 'G'), 'positions': {Position(genome='genome.1', chrom='NC_011353.1', pos=1743676), Position(genome='genome.2', chrom='NZ_CP009644.1', pos=2929307), Position(genome='genome.1', chrom='NC_011353.1', pos=3772004)}}}, 3772004: {('genome.1', 'NC_011353.1', 'A', 'G'): {'alleles': ('A', 'G'), 'positions': {Position(genome='genome.1', chrom='NC_011353.1', pos=1743676), Position(genome='genome.2', chrom='NZ_CP009644.1', pos=2929307), Position(genome='genome.1', chrom='NC_011353.1', pos=3772004)}}}}
#
#     assert actual == expected
#
# def test_PositionedSNPsIndex_add_PositionedSNP_both_previous_PositionedSNPs_exists_but_are_not_equal():
#     position_1 = Position(genome='genome.1', chrom='NC_011353.1', pos=2372140)
#     position_2 = Position(genome='genome.3', chrom='NZ_CP016628.1', pos=2720609)
#     position_3 = Position(genome='genome.2', chrom='NZ_CP009644.1', pos=1846009)
#     position_4 = Position(genome='genome.1', chrom='NC_011353.1', pos=2157387)
#     allele_1 = 'C'
#     allele_2 = 'T'
#     positionedSNPsIndex = PositionedSNPsIndex(length_of_longest_contig)
#     positionedSNPsIndex.add_PositionedSNP(position_1, position_3, allele_1, allele_2)
#     positionedSNPsIndex.add_PositionedSNP(position_2, position_4, allele_1, allele_2)
#     positionedSNPsIndex.add_PositionedSNP(position_1, position_2, allele_1, allele_2)
#
#     actual = positionedSNPsIndex.to_dict()
#     expected = {1846009: {('genome.2', 'NZ_CP009644.1', 'C', 'T'): {'alleles': ('C', 'T'), 'positions': {Position(genome='genome.3', chrom='NZ_CP016628.1', pos=2720609), Position(genome='genome.1', chrom='NC_011353.1', pos=2157387), Position(genome='genome.1', chrom='NC_011353.1', pos=2372140), Position(genome='genome.2', chrom='NZ_CP009644.1', pos=1846009)}}}, 2157387: {('genome.1', 'NC_011353.1', 'C', 'T'): {'alleles': ('C', 'T'), 'positions': {Position(genome='genome.3', chrom='NZ_CP016628.1', pos=2720609), Position(genome='genome.1', chrom='NC_011353.1', pos=2157387), Position(genome='genome.1', chrom='NC_011353.1', pos=2372140), Position(genome='genome.2', chrom='NZ_CP009644.1', pos=1846009)}}}, 2372140: {('genome.1', 'NC_011353.1', 'C', 'T'): {'alleles': ('C', 'T'), 'positions': {Position(genome='genome.3', chrom='NZ_CP016628.1', pos=2720609), Position(genome='genome.1', chrom='NC_011353.1', pos=2157387), Position(genome='genome.1', chrom='NC_011353.1', pos=2372140), Position(genome='genome.2', chrom='NZ_CP009644.1', pos=1846009)}}}, 2720609: {('genome.3', 'NZ_CP016628.1', 'C', 'T'): {'alleles': ('C', 'T'), 'positions': {Position(genome='genome.3', chrom='NZ_CP016628.1', pos=2720609), Position(genome='genome.1', chrom='NC_011353.1', pos=2157387), Position(genome='genome.1', chrom='NC_011353.1', pos=2372140), Position(genome='genome.2', chrom='NZ_CP009644.1', pos=1846009)}}}}
#
#     assert actual == expected
#
#
#
# def test_PositionedSNPsIndex_add_PositionedSNP_both_previous_PositionedSNPs_exists_and_are_equal():
#     position_1 = Position(genome='genome.1', chrom='NC_011353.1', pos=2372140)
#     position_2 = Position(genome='genome.3', chrom='NZ_CP016628.1', pos=2720609)
#     position_3 = Position(genome='genome.2', chrom='NZ_CP009644.1', pos=1846009)
#     position_4 = Position(genome='genome.1', chrom='NC_011353.1', pos=2157387)
#     allele_1 = 'C'
#     allele_2 = 'T'
#     positionedSNPsIndex = PositionedSNPsIndex(length_of_longest_contig)
#     positionedSNPsIndex.add_PositionedSNP(position_1, position_3, allele_1, allele_2)
#     positionedSNPsIndex.add_PositionedSNP(position_2, position_4, allele_1, allele_2)
#     positionedSNPsIndex.add_PositionedSNP(position_1, position_2, allele_1, allele_2)
#     positionedSNPsIndex.add_PositionedSNP(position_3, position_4, allele_1, allele_2)
#
#     actual = positionedSNPsIndex.to_dict()
#     expected = {1846009: {('genome.2', 'NZ_CP009644.1', 'C', 'T'): {'alleles': ('C', 'T'), 'positions': {Position(genome='genome.3', chrom='NZ_CP016628.1', pos=2720609), Position(genome='genome.1', chrom='NC_011353.1', pos=2157387), Position(genome='genome.1', chrom='NC_011353.1', pos=2372140), Position(genome='genome.2', chrom='NZ_CP009644.1', pos=1846009)}}}, 2157387: {('genome.1', 'NC_011353.1', 'C', 'T'): {'alleles': ('C', 'T'), 'positions': {Position(genome='genome.3', chrom='NZ_CP016628.1', pos=2720609), Position(genome='genome.1', chrom='NC_011353.1', pos=2157387), Position(genome='genome.1', chrom='NC_011353.1', pos=2372140), Position(genome='genome.2', chrom='NZ_CP009644.1', pos=1846009)}}}, 2372140: {('genome.1', 'NC_011353.1', 'C', 'T'): {'alleles': ('C', 'T'), 'positions': {Position(genome='genome.3', chrom='NZ_CP016628.1', pos=2720609), Position(genome='genome.1', chrom='NC_011353.1', pos=2157387), Position(genome='genome.1', chrom='NC_011353.1', pos=2372140), Position(genome='genome.2', chrom='NZ_CP009644.1', pos=1846009)}}}, 2720609: {('genome.3', 'NZ_CP016628.1', 'C', 'T'): {'alleles': ('C', 'T'), 'positions': {Position(genome='genome.3', chrom='NZ_CP016628.1', pos=2720609), Position(genome='genome.1', chrom='NC_011353.1', pos=2157387), Position(genome='genome.1', chrom='NC_011353.1', pos=2372140), Position(genome='genome.2', chrom='NZ_CP009644.1', pos=1846009)}}}}
#
#     assert actual == expected
#
#
# def loadPositionedSNPsIndexFromCSVs():
#     positionedSNPsIndex = PositionedSNPsIndex(length_of_longest_contig)
#     positionedSNPsIndex.add_SNPs_from_csv('test_cases/Positioned_SNPs/genome.0-SEP-genome.1.mummer.filtered.csv', 'genome.0', 'genome.1')
#     positionedSNPsIndex.add_SNPs_from_csv('test_cases/Positioned_SNPs/genome.0-SEP-genome.2.mummer.filtered.csv', 'genome.0', 'genome.2')
#     positionedSNPsIndex.add_SNPs_from_csv('test_cases/Positioned_SNPs/genome.0-SEP-genome.3.mummer.filtered.csv', 'genome.0', 'genome.3')
#     positionedSNPsIndex.add_SNPs_from_csv('test_cases/Positioned_SNPs/genome.1-SEP-genome.2.mummer.filtered.csv', 'genome.1', 'genome.2')
#     positionedSNPsIndex.add_SNPs_from_csv('test_cases/Positioned_SNPs/genome.1-SEP-genome.3.mummer.filtered.csv', 'genome.1', 'genome.3')
#     positionedSNPsIndex.add_SNPs_from_csv('test_cases/Positioned_SNPs/genome.2-SEP-genome.3.mummer.filtered.csv', 'genome.2', 'genome.3')
#     return positionedSNPsIndex
#
# def test_PositionedSNPsIndex_check_if_there_are_no_PositionedSNP_copies():
#     '''
#     Ensures that there are no PositionedSNP copies.
#     TODO: improve this, this is a very long test
#     '''
#     positionedSNPsIndex = loadPositionedSNPsIndexFromCSVs()
#     all_positioned_SNPs = []
#     for dict_in_pos in positionedSNPsIndex._PositionedSNPIndexKey_to_PositionedSNP:
#         if dict_in_pos is not None:
#             for positioned_snp in dict_in_pos.values():
#                 all_positioned_SNPs.append(positioned_snp)
#
#     for positioned_snp_1 in all_positioned_SNPs:
#         for positioned_snp_2 in all_positioned_SNPs:
#             assert (positioned_snp_1 == positioned_snp_2 and id(positioned_snp_1) == id(positioned_snp_2)) or \
#                    (positioned_snp_1 != positioned_snp_2 and id(positioned_snp_1) != id(positioned_snp_2))
#
#
# def test_PositionedSNPsIndex_load_save():
#     positionedSNPsIndex = loadPositionedSNPsIndexFromCSVs()
#     positionedSNPsIndex.save('test_cases/Positioned_SNPs/positionedSNPsIndex')
#     loaded_positionedSNPsIndex = PositionedSNPsIndex.load('test_cases/Positioned_SNPs/positionedSNPsIndex')
#     assert positionedSNPsIndex == loaded_positionedSNPsIndex
#
# def test_PositionedSNPsIndex_get_all_unique_Positioned_SNPs():
#     '''
#     TODO: improve this, this is a very long test
#     '''
#     positionedSNPsIndex = loadPositionedSNPsIndexFromCSVs()
#     all_unique_positioned_SNPs = positionedSNPsIndex.get_all_unique_Positioned_SNPs()
#     for i, positioned_snp_1 in enumerate(all_unique_positioned_SNPs):
#         for j, positioned_snp_2 in enumerate(all_unique_positioned_SNPs):
#             if i != j:
#                 assert positioned_snp_1 != positioned_snp_2 and id(positioned_snp_1) != id(positioned_snp_2)
#
#
# def test_get_nb_SNPs_that_can_be_found_with_the_given_genomes():
#     pass  # TODO