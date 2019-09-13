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












TEST(PangenomeSNPsMapWithFastPositionAccess, set_PangenomeSNP_sets_not_indexed_position_expects_out_of_range_exception) {
    PangenomeSNPsMapWithFastPositionAccess pangenomeSNPsMap;
    PangenomeSNP_Pointer pangenomeSnpPointer;
    PangenomeSNPsIndexKeyType key;
    std::tie(key, pangenomeSnpPointer) = create_key_and_pangenome_SNP_pointer();

    try {
        pangenomeSNPsMap.set_PangenomeSNP(key, pangenomeSnpPointer);
        FAIL() << "Expected std::out_of_range";
    }
    catch(std::out_of_range const & err) {}
    catch(...) {
        FAIL() << "Expected std::out_of_range";
    }
}


TEST(PangenomeSNPsMapWithFastPositionAccess, set_PangenomeSNP_sets_indexed_position) {
    PangenomeSNPsMapWithFastPositionAccess pangenomeSNPsMap(10);
    PangenomeSNP_Pointer pangenomeSnpPointer;
    PangenomeSNPsIndexKeyType key;
    std::tie(key, pangenomeSnpPointer) = create_key_and_pangenome_SNP_pointer();

    pangenomeSNPsMap.set_PangenomeSNP(key, pangenomeSnpPointer);

    auto expected = pangenomeSnpPointer;
    auto actual = pangenomeSNPsMap.get_PangenomeSNP(key);

    EXPECT_EQ(expected, actual);
}


TEST(PangenomeSNPsMapWithFastPositionAccess, set_PangenomeSNP_sets_indexed_position_which_already_has_one_SNP) {
    PangenomeSNPsMapWithFastPositionAccess pangenomeSNPsMap(10);
    PangenomeSNP_Pointer pangenomeSnpPointer_1;
    PangenomeSNPsIndexKeyType key_1;
    std::tie(key_1, pangenomeSnpPointer_1) = create_key_and_pangenome_SNP_pointer();
    PangenomeSNP_Pointer pangenomeSnpPointer_2;
    PangenomeSNPsIndexKeyType key_2;
    std::tie(key_2, pangenomeSnpPointer_2) = create_key_and_pangenome_SNP_pointer('A', 'T', "genome_2");


    pangenomeSNPsMap.set_PangenomeSNP(key_1, pangenomeSnpPointer_1);
    pangenomeSNPsMap.set_PangenomeSNP(key_2, pangenomeSnpPointer_2);

    auto expected = pangenomeSnpPointer_2;
    auto actual = pangenomeSNPsMap.get_PangenomeSNP(key_2);

    EXPECT_EQ(expected, actual);
}

TEST(PangenomeSNPsMapWithFastPositionAccess, set_PangenomeSNP_resets_existing_position_with_a_different_SNP) {
    PangenomeSNPsMapWithFastPositionAccess pangenomeSNPsMap(10);
    PangenomeSNP_Pointer pangenomeSnpPointer;
    PangenomeSNPsIndexKeyType key;
    std::tie(key, pangenomeSnpPointer) = create_key_and_pangenome_SNP_pointer();
    PangenomeSNP_Pointer same_pangenomeSnpPointer_with_different_alleles_but_same_position;
    std::tie(std::ignore, same_pangenomeSnpPointer_with_different_alleles_but_same_position) = create_key_and_pangenome_SNP_pointer('C', 'G', "genome_1", "chrom_1", 1);

    pangenomeSNPsMap.set_PangenomeSNP(key, pangenomeSnpPointer);
    pangenomeSNPsMap.set_PangenomeSNP(key, same_pangenomeSnpPointer_with_different_alleles_but_same_position);

    auto expected = same_pangenomeSnpPointer_with_different_alleles_but_same_position;
    auto actual = pangenomeSNPsMap.get_PangenomeSNP(key);

    EXPECT_EQ(expected, actual);
}

TEST(PangenomeSNPsMapWithFastPositionAccess, set_PangenomeSNP_sets_two_different_positions) {
    PangenomeSNPsMapWithFastPositionAccess pangenomeSNPsMap(10);
    PangenomeSNP_Pointer pangenomeSnpPointer_1;
    PangenomeSNPsIndexKeyType key_1;
    std::tie(key_1, pangenomeSnpPointer_1) = create_key_and_pangenome_SNP_pointer();
    PangenomeSNP_Pointer pangenomeSnpPointer_2;
    PangenomeSNPsIndexKeyType key_2;
    std::tie(key_2, pangenomeSnpPointer_2) = create_key_and_pangenome_SNP_pointer('C', 'T', "genome_2", "chrom_2", 2);


    pangenomeSNPsMap.set_PangenomeSNP(key_1, pangenomeSnpPointer_1);
    pangenomeSNPsMap.set_PangenomeSNP(key_2, pangenomeSnpPointer_2);

    auto expected_1 = pangenomeSnpPointer_1;
    auto actual_1 = pangenomeSNPsMap.get_PangenomeSNP(key_1);
    EXPECT_EQ(expected_1, actual_1);

    auto expected_2 = pangenomeSnpPointer_2;
    auto actual_2 = pangenomeSNPsMap.get_PangenomeSNP(key_2);
    EXPECT_EQ(expected_2, actual_2);
}






TEST(PangenomeSNPsMapWithFastPositionAccess, get_PangenomeSNP_gets_not_indexed_position_expects_out_of_range_exception) {
    PangenomeSNPsMapWithFastPositionAccess pangenomeSNPsMap;
    PangenomeSNPsIndexKeyType key;
    std::tie(key, std::ignore) = create_key_and_pangenome_SNP_pointer();

    try {
        pangenomeSNPsMap.get_PangenomeSNP(key);
        FAIL() << "Expected std::out_of_range";
    }
    catch(std::out_of_range const & err) {}
    catch(...) {
        FAIL() << "Expected std::out_of_range";
    }
}


TEST(PangenomeSNPsMapWithFastPositionAccess, get_PangenomeSNP_position_exists_key_does_not_exist) {
    PangenomeSNPsMapWithFastPositionAccess pangenomeSNPsMap(10);
    PangenomeSNPsIndexKeyType key;
    std::tie(key, std::ignore) = create_key_and_pangenome_SNP_pointer();

    PangenomeSNP_Pointer expected = nullptr;
    PangenomeSNP_Pointer actual = pangenomeSNPsMap.get_PangenomeSNP(key);

    EXPECT_EQ(expected, actual);
}


TEST(PangenomeSNPsMapWithFastPositionAccess, merge_two_PangenomeSNPs) {
    PangenomeSNPsMapWithFastPositionAccess pangenomeSNPsMap(10);
    PangenomeSNPsIndexKeyType key_1;
    PangenomeSNP_Pointer pangenomeSnpPointer_1;
    std::tie(key_1, pangenomeSnpPointer_1) = create_key_and_pangenome_SNP_pointer();
    PangenomeSNPsIndexKeyType key_2;
    PangenomeSNP_Pointer pangenomeSnpPointer_2;
    std::tie(key_2, pangenomeSnpPointer_2) = create_key_and_pangenome_SNP_pointer('A', 'C', "genome_2", "chrom_2", 2);

    pangenomeSNPsMap.merge_two_PangenomeSNPs(pangenomeSnpPointer_1, pangenomeSnpPointer_2);
}