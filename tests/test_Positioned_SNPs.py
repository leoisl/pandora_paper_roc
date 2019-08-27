import pytest

from evaluate.Positioned_SNPs import *


def test_PositionedSNP_merge():
    positioned_SNP_1 = PositionedSNP.from_dict({'alleles': ['C', 'T'], 'positions': {Position(genome='genome.1', chrom='NC_011353.1', pos=2372140), Position(genome='genome.2', chrom='NZ_CP009644.1', pos=1846009)}})
    positioned_SNP_2 = PositionedSNP.from_dict({'alleles': ['C', 'T'], 'positions': {Position(genome='genome.3', chrom='NZ_CP016628.1', pos=2720609), Position(genome='genome.1', chrom='NC_011353.1', pos=2157387)}})

    actual = positioned_SNP_1.merge(positioned_SNP_2)
    expected = PositionedSNP.from_dict({'alleles': ['C', 'T'], 'positions': {Position(genome='genome.3', chrom='NZ_CP016628.1', pos=2720609), Position(genome='genome.2', chrom='NZ_CP009644.1', pos=1846009), Position(genome='genome.1', chrom='NC_011353.1', pos=2157387), Position(genome='genome.1', chrom='NC_011353.1', pos=2372140)}})

    assert actual == expected


def test_PositionedSNPsIndex_add_PositionedSNP_no_previous_PositionedSNP():
    positionedSNPsIndex = PositionedSNPsIndex()
    position_1 = Position(genome='genome.1', chrom='NC_011350.1', pos=75224)
    position_2 = Position(genome='genome.2', chrom='NZ_CP009644.1', pos=4166248)
    allele_1 = "G"
    allele_2 = "T"
    positionedSNPsIndex.add_PositionedSNP(position_1, position_2, allele_1, allele_2)

    positionedSNP = PositionedSNP.from_dict({'alleles': ['G', 'T'], 'positions': {position_2, position_1}})
    expected = PositionedSNPsIndex.from_dict({'PositionedSNPIndexKey_to_PositionedSNP':{
        PositionedSNPIndexKeyType(position=position_1, allele_1='G', allele_2='T'): positionedSNP,
        PositionedSNPIndexKeyType(position=position_2, allele_1='G', allele_2='T'): positionedSNP}})

    assert positionedSNPsIndex == expected

def test_PositionedSNPsIndex_add_PositionedSNP_only_previous_PositionedSNP_1_exists():
    position_1 = Position(genome='genome.1', chrom='NC_011353.1', pos=3772004)
    position_2 = Position(genome='genome.2', chrom='NZ_CP009644.1', pos=2929307)
    position_3 = Position(genome='genome.1', chrom='NC_011353.1', pos=1743676)
    allele_1 = "A"
    allele_2 = "G"
    positionedSNPIndexKey_1 = PositionedSNPIndexKeyType(position_1, allele_1, allele_2)
    positionedSNPIndexKey_2 = PositionedSNPIndexKeyType(position_2, allele_1, allele_2)
    positionedSNPIndexKey_3 = PositionedSNPIndexKeyType(position_3, allele_1, allele_2)
    positionedSNP = PositionedSNP.from_dict({'alleles': ['A', 'G'], 'positions': {position_2, position_3}})

    positionedSNPsIndex = PositionedSNPsIndex.from_dict({
        "PositionedSNPIndexKey_to_PositionedSNP": {
            positionedSNPIndexKey_2: positionedSNP,
            positionedSNPIndexKey_3: positionedSNP,
        }
    })
    positionedSNPsIndex.add_PositionedSNP(position_2, position_1, allele_1, allele_2)

    positionedSNP_expected = PositionedSNP.from_dict({'alleles': ['A', 'G'], 'positions': {position_1, position_2, position_3}})
    expected = PositionedSNPsIndex.from_dict({
        "PositionedSNPIndexKey_to_PositionedSNP": {
            positionedSNPIndexKey_1: positionedSNP_expected,
            positionedSNPIndexKey_2: positionedSNP_expected,
            positionedSNPIndexKey_3: positionedSNP_expected,
        }
    })

    assert positionedSNPsIndex == expected

def test_PositionedSNPsIndex_add_PositionedSNP_only_previous_PositionedSNP_2_exists():
    position_1 = Position(genome='genome.1', chrom='NC_011353.1', pos=3772004)
    position_2 = Position(genome='genome.2', chrom='NZ_CP009644.1', pos=2929307)
    position_3 = Position(genome='genome.1', chrom='NC_011353.1', pos=1743676)
    allele_1 = "A"
    allele_2 = "G"
    positionedSNPIndexKey_1 = PositionedSNPIndexKeyType(position_1, allele_1, allele_2)
    positionedSNPIndexKey_2 = PositionedSNPIndexKeyType(position_2, allele_1, allele_2)
    positionedSNPIndexKey_3 = PositionedSNPIndexKeyType(position_3, allele_1, allele_2)
    positionedSNP = PositionedSNP.from_dict({'alleles': ['A', 'G'], 'positions': {position_2, position_3}})

    positionedSNPsIndex = PositionedSNPsIndex.from_dict({
        "PositionedSNPIndexKey_to_PositionedSNP": {
            positionedSNPIndexKey_2: positionedSNP,
            positionedSNPIndexKey_3: positionedSNP,
        }
    })
    positionedSNPsIndex.add_PositionedSNP(position_1, position_2, allele_1, allele_2)

    positionedSNP_expected = PositionedSNP.from_dict({'alleles': ['A', 'G'], 'positions': {position_1, position_2, position_3}})
    expected = PositionedSNPsIndex.from_dict({
        "PositionedSNPIndexKey_to_PositionedSNP": {
            positionedSNPIndexKey_1: positionedSNP_expected,
            positionedSNPIndexKey_2: positionedSNP_expected,
            positionedSNPIndexKey_3: positionedSNP_expected,
        }
    })

    assert positionedSNPsIndex == expected

def test_PositionedSNPsIndex_add_PositionedSNP_both_previous_PositionedSNPs_exists_but_are_not_equal():
    position_1 = Position(genome='genome.1', chrom='NC_011353.1', pos=2372140)
    position_2 = Position(genome='genome.3', chrom='NZ_CP016628.1', pos=2720609)
    position_3 = Position(genome='genome.2', chrom='NZ_CP009644.1', pos=1846009)
    position_4 = Position(genome='genome.1', chrom='NC_011353.1', pos=2157387)
    allele_1 = 'C'
    allele_2 = 'T'
    positionedSNPIndexKey_1 = PositionedSNPIndexKeyType(position_1, allele_1, allele_2)
    positionedSNPIndexKey_2 = PositionedSNPIndexKeyType(position_2, allele_1, allele_2)
    positionedSNPIndexKey_3 = PositionedSNPIndexKeyType(position_3, allele_1, allele_2)
    positionedSNPIndexKey_4 = PositionedSNPIndexKeyType(position_4, allele_1, allele_2)
    positionedSNP_1 = PositionedSNP.from_dict(
        {'alleles': ['C', 'T'],
         'positions': {position_1, position_3}})
    positionedSNP_2 = PositionedSNP.from_dict(
        {'alleles': ['C', 'T'],
         'positions': {position_2, position_4}})
    positionedSNP_all = PositionedSNP.from_dict(
        {'alleles': ['C', 'T'],
         'positions': {position_1, position_2, position_3, position_4}})
    positionedSNPsIndex = PositionedSNPsIndex.from_dict({
        "PositionedSNPIndexKey_to_PositionedSNP": {
            positionedSNPIndexKey_1: positionedSNP_1,
            positionedSNPIndexKey_2: positionedSNP_2,
            positionedSNPIndexKey_3: positionedSNP_1,
            positionedSNPIndexKey_4: positionedSNP_2
        }
    })
    positionedSNPsIndex.add_PositionedSNP(position_1, position_2, allele_1, allele_2)

    expected = PositionedSNPsIndex.from_dict({
        "PositionedSNPIndexKey_to_PositionedSNP": {
            positionedSNPIndexKey_1: positionedSNP_all,
            positionedSNPIndexKey_2: positionedSNP_all,
            positionedSNPIndexKey_3: positionedSNP_all,
            positionedSNPIndexKey_4: positionedSNP_all,
        }
    })

    assert positionedSNPsIndex == expected



def test_PositionedSNPsIndex_add_PositionedSNP_both_previous_PositionedSNPs_exists_and_are_equal():
    position_1 = Position(genome='genome.1', chrom='NC_011353.1', pos=2372140)
    position_2 = Position(genome='genome.3', chrom='NZ_CP016628.1', pos=2720609)
    position_3 = Position(genome='genome.2', chrom='NZ_CP009644.1', pos=1846009)
    position_4 = Position(genome='genome.1', chrom='NC_011353.1', pos=2157387)
    allele_1 = 'C'
    allele_2 = 'T'
    positionedSNPIndexKey_1 = PositionedSNPIndexKeyType(position_1, allele_1, allele_2)
    positionedSNPIndexKey_2 = PositionedSNPIndexKeyType(position_2, allele_1, allele_2)
    positionedSNPIndexKey_3 = PositionedSNPIndexKeyType(position_3, allele_1, allele_2)
    positionedSNPIndexKey_4 = PositionedSNPIndexKeyType(position_4, allele_1, allele_2)
    positionedSNP_all = PositionedSNP.from_dict(
        {'alleles': ['C', 'T'],
         'positions': {position_1, position_2, position_3, position_4}})
    positionedSNPsIndex = PositionedSNPsIndex.from_dict({
        "PositionedSNPIndexKey_to_PositionedSNP": {
            positionedSNPIndexKey_1: positionedSNP_all,
            positionedSNPIndexKey_2: positionedSNP_all,
            positionedSNPIndexKey_3: positionedSNP_all,
            positionedSNPIndexKey_4: positionedSNP_all
        }
    })
    positionedSNPsIndex.add_PositionedSNP(position_1, position_2, allele_1, allele_2)

    expected = PositionedSNPsIndex.from_dict({
        "PositionedSNPIndexKey_to_PositionedSNP": {
            positionedSNPIndexKey_1: positionedSNP_all,
            positionedSNPIndexKey_2: positionedSNP_all,
            positionedSNPIndexKey_3: positionedSNP_all,
            positionedSNPIndexKey_4: positionedSNP_all,
        }
    })

    assert positionedSNPsIndex == expected

def test_PositionedSNPsIndex_get_nb_SNPs_in_pangenome():
    position_1 = Position(genome='genome.1', chrom='NC_011353.1', pos=2372140)
    position_2 = Position(genome='genome.3', chrom='NZ_CP016628.1', pos=2720609)
    position_3 = Position(genome='genome.2', chrom='NZ_CP009644.1', pos=1846009)
    position_4 = Position(genome='genome.1', chrom='NC_011353.1', pos=2157387)
    allele_1 = 'C'
    allele_2 = 'T'
    positionedSNPsIndex = PositionedSNPsIndex()
    positionedSNPsIndex.add_PositionedSNP(position_1, position_2, allele_1, allele_2)
    positionedSNPsIndex.add_PositionedSNP(position_3, position_4, allele_1, allele_2)

    assert positionedSNPsIndex.get_nb_SNPs_in_pangenome() == 2
    pass


def test_PositionedSNPsIndex_check_if_there_are_no_PositionedSNP_copies():
    '''
    Ensures that there are no PositionedSNP copies.
    This does not need to be called in the true run
    '''
    positionedSNPsIndex = PositionedSNPsIndex()
    positionedSNPsIndex.add_SNPs_from_csv('test_cases/Positioned_SNPs/genome.0-SEP-genome.1.mummer.filtered.csv', 'genome.0', 'genome.1')
    positionedSNPsIndex.add_SNPs_from_csv('test_cases/Positioned_SNPs/genome.0-SEP-genome.2.mummer.filtered.csv', 'genome.0', 'genome.2')
    positionedSNPsIndex.add_SNPs_from_csv('test_cases/Positioned_SNPs/genome.0-SEP-genome.3.mummer.filtered.csv', 'genome.0', 'genome.3')
    positionedSNPsIndex.add_SNPs_from_csv('test_cases/Positioned_SNPs/genome.1-SEP-genome.2.mummer.filtered.csv', 'genome.1', 'genome.2')
    positionedSNPsIndex.add_SNPs_from_csv('test_cases/Positioned_SNPs/genome.1-SEP-genome.3.mummer.filtered.csv', 'genome.1', 'genome.3')
    positionedSNPsIndex.add_SNPs_from_csv('test_cases/Positioned_SNPs/genome.2-SEP-genome.3.mummer.filtered.csv', 'genome.2', 'genome.3')
    for positioned_snp_1 in positionedSNPsIndex.PositionedSNPIndexKey_to_PositionedSNP.values():
        for positioned_snp_2 in positionedSNPsIndex.PositionedSNPIndexKey_to_PositionedSNP.values():
            assert (positioned_snp_1 == positioned_snp_2 and positioned_snp_1 is positioned_snp_2) or \
                   (positioned_snp_1 != positioned_snp_2 and positioned_snp_1 is not positioned_snp_2)



def test_PositionedSNPsIndex_load_save():
    positionedSNPsIndex = PositionedSNPsIndex()
    positionedSNPsIndex.add_SNPs_from_csv('test_cases/Positioned_SNPs/genome.1-SEP-genome.2.mummer.30000.csv', 'genome.1', 'genome.2')
    positionedSNPsIndex.add_SNPs_from_csv('test_cases/Positioned_SNPs/genome.0-SEP-genome.1.mummer.30000.csv', 'genome.0', 'genome.1')
    positionedSNPsIndex.add_SNPs_from_csv('test_cases/Positioned_SNPs/genome.0-SEP-genome.2.mummer.30000.csv', 'genome.0', 'genome.2')
    positionedSNPsIndex.save('test_cases/Positioned_SNPs/positionedSNPsIndex')
    loaded_positionedSNPsIndex = PositionedSNPsIndex.load('test_cases/Positioned_SNPs/positionedSNPsIndex')
    assert positionedSNPsIndex == loaded_positionedSNPsIndex