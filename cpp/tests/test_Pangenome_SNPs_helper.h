#include <tuple>
#include "Pangenome_SNPs.h"
#include <boost/smart_ptr/make_shared.hpp>
using namespace std;

std::tuple<PangenomeSNP, PangenomeSNP> make_two_PangenomeSNPs_with_two_different_alleles() {
    PangenomeSNP pangenomeSnp_1('A', 'C');
    PangenomeSNP pangenomeSnp_2('G', 'T');
    return make_tuple(pangenomeSnp_1, pangenomeSnp_2);
}

std::tuple<PangenomeSNP, PangenomeSNP> make_two_PangenomeSNPs_with_one_different_allele() {
    PangenomeSNP pangenomeSnp_1('A', 'C');
    PangenomeSNP pangenomeSnp_2('A', 'T');
    return make_tuple(pangenomeSnp_1, pangenomeSnp_2);
}

std::tuple<PangenomeSNP, PangenomeSNP> make_two_PangenomeSNPs_with_no_shared_Genomic_Coordinates() {
    GenomicCoordinate genomic_coordinate_1("genome_1", "chrom_1", 1);
    GenomicCoordinate genomic_coordinate_2("genome_2", "chrom_2", 2);
    GenomicCoordinate genomic_coordinate_3("genome_3", "chrom_3", 3);
    GenomicCoordinate genomic_coordinate_4("genome_4", "chrom_4", 4);

    PangenomeSNP pangenomeSnp_1('A', 'G');
    pangenomeSnp_1.add_genomic_coordinate(genomic_coordinate_1);
    pangenomeSnp_1.add_genomic_coordinate(genomic_coordinate_2);
    PangenomeSNP pangenomeSnp_2('A', 'G');
    pangenomeSnp_2.add_genomic_coordinate(genomic_coordinate_3);
    pangenomeSnp_2.add_genomic_coordinate(genomic_coordinate_4);

    return make_tuple(pangenomeSnp_1, pangenomeSnp_2);
}

std::tuple<PangenomeSNP, PangenomeSNP, PangenomeSNP> make_two_PangenomeSNPs_with_one_shared_Genomic_Coordinates_and_expected() {
    GenomicCoordinate genomic_coordinate_1("genome_1", "chrom_1", 1);
    GenomicCoordinate genomic_coordinate_2("genome_2", "chrom_2", 2);
    GenomicCoordinate genomic_coordinate_3("genome_3", "chrom_3", 3);

    PangenomeSNP pangenomeSnp_1('A', 'G');
    pangenomeSnp_1.add_genomic_coordinate(genomic_coordinate_1);
    pangenomeSnp_1.add_genomic_coordinate(genomic_coordinate_2);
    PangenomeSNP pangenomeSnp_2('A', 'G');
    pangenomeSnp_2.add_genomic_coordinate(genomic_coordinate_2);
    pangenomeSnp_2.add_genomic_coordinate(genomic_coordinate_3);

    PangenomeSNP expected('A', 'G');
    expected.add_genomic_coordinate(genomic_coordinate_1);
    expected.add_genomic_coordinate(genomic_coordinate_2);
    expected.add_genomic_coordinate(genomic_coordinate_3);

    return make_tuple(pangenomeSnp_1, pangenomeSnp_2, expected);
}

std::tuple<PangenomeSNP, PangenomeSNP> make_two_PangenomeSNPs_with_all_four_shared_Genomic_Coordinates() {
    GenomicCoordinate genomic_coordinate_1("genome_1", "chrom_1", 1);
    GenomicCoordinate genomic_coordinate_2("genome_2", "chrom_2", 2);
    GenomicCoordinate genomic_coordinate_3("genome_3", "chrom_3", 3);
    GenomicCoordinate genomic_coordinate_4("genome_4", "chrom_4", 4);

    PangenomeSNP pangenomeSnp_1('A', 'G');
    pangenomeSnp_1.add_genomic_coordinate(genomic_coordinate_1);
    pangenomeSnp_1.add_genomic_coordinate(genomic_coordinate_2);
    pangenomeSnp_1.add_genomic_coordinate(genomic_coordinate_3);
    pangenomeSnp_1.add_genomic_coordinate(genomic_coordinate_4);
    PangenomeSNP pangenomeSnp_2('A', 'G');
    pangenomeSnp_2.add_genomic_coordinate(genomic_coordinate_1);
    pangenomeSnp_2.add_genomic_coordinate(genomic_coordinate_2);
    pangenomeSnp_2.add_genomic_coordinate(genomic_coordinate_3);
    pangenomeSnp_2.add_genomic_coordinate(genomic_coordinate_4);

    return make_tuple(pangenomeSnp_1, pangenomeSnp_2);
}

PangenomeSNP make_one_PangenomeSNP_with_one_Genomic_Coordinate() {
    GenomicCoordinate genomic_coordinate_1("genome_1", "chrom_1", 1);
    PangenomeSNP pangenomeSnp_1('A', 'G');
    pangenomeSnp_1.add_genomic_coordinate(genomic_coordinate_1);
    return pangenomeSnp_1;
}


PangenomeSNP make_one_PangenomeSNP_with_two_Genomic_Coordinates() {
    GenomicCoordinate genomic_coordinate_1("genome_1", "chrom_1", 1);
    GenomicCoordinate genomic_coordinate_2("genome_2", "chrom_2", 2);
    PangenomeSNP pangenomeSnp_1('A', 'G');
    pangenomeSnp_1.add_genomic_coordinate(genomic_coordinate_1);
    pangenomeSnp_1.add_genomic_coordinate(genomic_coordinate_2);
    return pangenomeSnp_1;
}

std::tuple<PangenomeSNPsIndexKeyType, PangenomeSNP_Pointer> create_key_and_pangenome_SNP_pointer(char allele_1='A', char allele_2 = 'T',
        const std::string &genome = "genome_1", const std::string &chrom = "chrom_1", uint32_t position=1) {
    GenomicCoordinate genomic_coordinate(genome, chrom, position);
    PangenomeSNP pangenomeSnp(allele_1, allele_2);
    pangenomeSnp.add_genomic_coordinate(genomic_coordinate);
    PangenomeSNP_Pointer pangenomeSnpPointer = boost::make_shared<PangenomeSNP>(pangenomeSnp);
    PangenomeSNPsIndexKeyType key(genomic_coordinate, allele_1, allele_2);
    return make_tuple(key, pangenomeSnpPointer);
}

std::tuple<PangenomeSNPsMapWithFastPositionAccess, PangenomeSNPsIndexKeyType, PangenomeSNP_Pointer> create_PangenomeSNPsMapWithFastPositionAccess_with_one_key_and_return_key_and_pangenomeSNP() {
    char allele_1 = 'A';
    char allele_2 = 'T';

    PangenomeSNPsMapWithFastPositionAccess pangenomeSNPsMap(10);
    GenomicCoordinate genomic_coordinate_5("genome_5", "chrom_5", 5);
    PangenomeSNP pangenomeSnp(allele_1, allele_2);
    pangenomeSnp.add_genomic_coordinate(genomic_coordinate_5);
    PangenomeSNPsIndexKeyType key(genomic_coordinate_5, 'A', 'T');
    PangenomeSNP_Pointer pangenomeSnpPointer = boost::make_shared<PangenomeSNP>(pangenomeSnp);
    pangenomeSNPsMap.set_PangenomeSNP(key, pangenomeSnpPointer);

    return make_tuple(pangenomeSNPsMap, key, pangenomeSnpPointer);
}