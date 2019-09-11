#include <tuple>
#include "Pangenome_SNPs.h"
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