#include "gtest/gtest.h"
#include "Pangenome_SNPs.h"
#include "test_Pangenome_SNPs_helper.h"
#include <tuple>

TEST(PangenomeSNP, create_SNP_alleles_not_ordered_expect_death) {
    EXPECT_DEATH(PangenomeSNP('C', 'A'), "");
}

TEST(PangenomeSNP, merge_two_SNPs_with_two_different_alleles_expect_death) {
    PangenomeSNP pangenomeSnp_1, pangenomeSnp_2;
    std::tie(pangenomeSnp_1, pangenomeSnp_2) = make_two_PangenomeSNPs_with_two_different_alleles();

    EXPECT_DEATH(pangenomeSnp_1.merge(pangenomeSnp_2), "");
}

TEST(PangenomeSNP, merge_two_PangenomeSNPs_with_one_different_allele_expect_death) {
    PangenomeSNP pangenomeSnp_1, pangenomeSnp_2;
    std::tie(pangenomeSnp_1, pangenomeSnp_2) = make_two_PangenomeSNPs_with_one_different_allele();

    EXPECT_DEATH(pangenomeSnp_1.merge(pangenomeSnp_2), "");
}


TEST(PangenomeSNP, merge_two_PangenomeSNPs_with_no_shared_Genomic_Coordinates_expect_death) {
    PangenomeSNP pangenomeSnp_1, pangenomeSnp_2;
    std::tie(pangenomeSnp_1, pangenomeSnp_2) = make_two_PangenomeSNPs_with_no_shared_Genomic_Coordinates();
    EXPECT_DEATH(pangenomeSnp_1.merge(pangenomeSnp_2), "");
}

TEST(PangenomeSNP, merge_two_PangenomeSNPs_with_one_shared_Genomic_Coordinates) {
    PangenomeSNP pangenomeSnp_1, pangenomeSnp_2, expected;
    std::tie(pangenomeSnp_1, pangenomeSnp_2, expected) = make_two_PangenomeSNPs_with_one_shared_Genomic_Coordinates_and_expected();

    PangenomeSNP pangenomeSnp_merged = pangenomeSnp_1.merge(pangenomeSnp_2);

    EXPECT_EQ(expected, pangenomeSnp_merged);
}


TEST(PangenomeSNP, merge_two_PangenomeSNPs_with_all_four_shared_Genomic_Coordinates) {
    PangenomeSNP pangenomeSnp_1, pangenomeSnp_2;
    std::tie(pangenomeSnp_1, pangenomeSnp_2) = make_two_PangenomeSNPs_with_all_four_shared_Genomic_Coordinates();

    PangenomeSNP pangenomeSnp_merged = pangenomeSnp_1.merge(pangenomeSnp_2);

    EXPECT_EQ(pangenomeSnp_1, pangenomeSnp_merged);
    EXPECT_EQ(pangenomeSnp_2, pangenomeSnp_merged);
}

TEST(PangenomeSNP, get_all_genomes_this_pangenome_snp_is_present_no_genomes) {
    PangenomeSNP pangenomeSnp;
    std::set<std::string> expected;

    auto actual = pangenomeSnp.get_all_genomes_this_pangenome_snp_is_present();

    EXPECT_EQ(expected, actual);
}

TEST(PangenomeSNP, get_all_genomes_this_pangenome_snp_is_present_one_genome) {
    PangenomeSNP pangenomeSnp = make_one_PangenomeSNP_with_one_Genomic_Coordinate();
    std::set<std::string> expected{"genome_1"};

    auto actual = pangenomeSnp.get_all_genomes_this_pangenome_snp_is_present();

    EXPECT_EQ(expected, actual);
}

TEST(PangenomeSNP, get_all_genomes_this_pangenome_snp_is_present_two_genomes) {
    PangenomeSNP pangenomeSnp = make_one_PangenomeSNP_with_two_Genomic_Coordinates();
    std::set<std::string> expected{"genome_1", "genome_2"};

    auto actual = pangenomeSnp.get_all_genomes_this_pangenome_snp_is_present();

    EXPECT_EQ(expected, actual);
}