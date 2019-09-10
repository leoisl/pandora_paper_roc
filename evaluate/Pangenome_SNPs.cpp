#include "Pangenome_SNPs.h"

PangenomeSNP PangenomeSNP::merge(const PangenomeSNP &other) const {
    assert(this->allele_1 == other.allele_1 && this->allele_2 == other.allele_2); // we should not merge PangenomeSNPs that do not have the same alleles
    PangenomeSNP new_PangenomeSNP(this->allele_1, this->allele_2);
    new_PangenomeSNP.genomic_coordinates.insert(this->genomic_coordinates.begin(), this->genomic_coordinates.end());
    new_PangenomeSNP.genomic_coordinates.insert(other.genomic_coordinates.begin(), other.genomic_coordinates.end());
    return new_PangenomeSNP;
}

std::set<std::string> PangenomeSNP::get_all_genomes_this_pangenome_snp_is_present() const {
    std::set<std::string> all_genomes_this_pangenome_snp_is_present;
    for (const GenomicCoordinate &genomic_coordinate : genomic_coordinates) {
        all_genomes_this_pangenome_snp_is_present.insert(genomic_coordinate.get_genome());
    }
    return all_genomes_this_pangenome_snp_is_present;
}



PangenomeSNP_Pointer PangenomeSNPsIndex::get_PangenomeSNP(const PangenomeSNPsIndexKeyType &key) const {
    /** @return: the PangenomeSNP on the coordinates given by pangenomeSNPsIndexKey, or nullptr if it does not exist */
    uint32_t position = key.get_genomic_coordinate().get_pos();
    Key_To_PangenomeSNP_Pointer key_To_PangenomeSNP_Pointer =
            position_To_Key_To_PangenomeSNP[position];

    if (key_To_PangenomeSNP_Pointer != nullptr)
        return key_To_PangenomeSNP_Pointer->at(key);
    else
        return nullptr;
}


void PangenomeSNPsIndex::set_PangenomeSNP(const PangenomeSNPsIndexKeyType &key, const PangenomeSNP_Pointer &pangenomeSNP) {
    uint32_t position = key.get_genomic_coordinate().get_pos();
    Key_To_PangenomeSNP_Pointer key_To_PangenomeSNP_Pointer =
            position_To_Key_To_PangenomeSNP[position];

    if (key_To_PangenomeSNP_Pointer == nullptr) {
        key_To_PangenomeSNP_Pointer.reset(new Key_To_PangenomeSNP());
    }
    (*key_To_PangenomeSNP_Pointer)[key] = pangenomeSNP;
}


PangenomeSNP_Pointer PangenomeSNPsIndex::merge_two_PangenomeSNPs(const PangenomeSNP_Pointer &previous_PangenomeSNP_in_key_of_snp_1,
                                             const PangenomeSNP_Pointer &previous_PangenomeSNP_in_key_of_snp_2,
                                             const char allele_1, const char allele_2) {
    PangenomeSNP_Pointer merged_SNP = boost::make_shared<PangenomeSNP>(previous_PangenomeSNP_in_key_of_snp_1->merge(*previous_PangenomeSNP_in_key_of_snp_2));

    for (const GenomicCoordinate &genomic_coordinate : merged_SNP->get_genomic_coordinates())
        set_PangenomeSNP(PangenomeSNPsIndexKeyType(genomic_coordinate, allele_1, allele_2), merged_SNP);

    return merged_SNP;
}

//TODO: there is just too much logic in here, maybe a candidate for refactoring
void PangenomeSNPsIndex::make_the_two_keys_point_to_the_same_pangenome_SNP (
        const PangenomeSNPsIndexKeyType &key_of_snp_1,
        const PangenomeSNPsIndexKeyType &key_of_snp_2,
        const char allele_1, const char allele_2) {
    /**
     * Comment of shame because I could not express myself well with code.
     * This function gets two keys of two SNPs and makes them both point to the same pangenome SNP.
     * It might need to create this pangenome SNP, or make one of the keys point to the an existing one, or merge
     * two existing PangenomeSNPs.
     * The main post-condition is that the index is consistent in regard to these two keys (i.e. the two keys will
     * point to the same PangenomeSNP).
     */

    PangenomeSNP_Pointer pangenomeSNP_pointed_by_both_keys(nullptr);

    PangenomeSNP_Pointer previous_PangenomeSNP_in_key_of_snp_1 = get_PangenomeSNP(key_of_snp_1);
    PangenomeSNP_Pointer previous_PangenomeSNP_in_key_of_snp_2 = get_PangenomeSNP(key_of_snp_2);

    bool need_to_create_a_new_Pangenome_SNP = previous_PangenomeSNP_in_key_of_snp_1 == nullptr && previous_PangenomeSNP_in_key_of_snp_2 == nullptr;
    bool PangenomeSNP_exists_for_key_1 = previous_PangenomeSNP_in_key_of_snp_1 != nullptr;
    bool PangenomeSNP_exists_for_key_2 = previous_PangenomeSNP_in_key_of_snp_2 != nullptr;
    bool Pangenome_SNP_exists_for_both_keys = PangenomeSNP_exists_for_key_1 and PangenomeSNP_exists_for_key_2;
    bool keys_point_to_the_same_object = previous_PangenomeSNP_in_key_of_snp_1 == previous_PangenomeSNP_in_key_of_snp_2;

    if (need_to_create_a_new_Pangenome_SNP) {
        pangenomeSNP_pointed_by_both_keys = boost::make_shared<PangenomeSNP>(allele_1, allele_2);
    }
    else {
        if (Pangenome_SNP_exists_for_both_keys) {
            if (keys_point_to_the_same_object) {
                pangenomeSNP_pointed_by_both_keys = previous_PangenomeSNP_in_key_of_snp_1;
            }else {
                // Both keys exist but they do not point to the same object
                // This is an inconsistent state: merge both these pangenome SNPs so that they point to the same object now
                pangenomeSNP_pointed_by_both_keys = merge_two_PangenomeSNPs(previous_PangenomeSNP_in_key_of_snp_1, previous_PangenomeSNP_in_key_of_snp_2,
                                                                            allele_1, allele_2);
            }
        }
        else if (PangenomeSNP_exists_for_key_1) {
            pangenomeSNP_pointed_by_both_keys = previous_PangenomeSNP_in_key_of_snp_1;
        }else if (PangenomeSNP_exists_for_key_2) {
            pangenomeSNP_pointed_by_both_keys = previous_PangenomeSNP_in_key_of_snp_2;
        }
    }

    assert(pangenomeSNP_pointed_by_both_keys != nullptr);
    set_PangenomeSNP(key_of_snp_1, pangenomeSNP_pointed_by_both_keys);
    set_PangenomeSNP(key_of_snp_2, pangenomeSNP_pointed_by_both_keys);
    pangenomeSNP_pointed_by_both_keys->add_genomic_coordinate(key_of_snp_1.get_genomic_coordinate());
    pangenomeSNP_pointed_by_both_keys->add_genomic_coordinate(key_of_snp_2.get_genomic_coordinate());

}


std::set<PangenomeSNP_Raw_Pointer> PangenomeSNPsIndex::get_all_unique_Pangenome_SNPs_as_raw_pointers() const {
    std::set<PangenomeSNP_Raw_Pointer> all_unique_Pangenome_SNPs_as_raw_pointers;
    for (const auto &key_To_PangenomeSNP_Pointer : position_To_Key_To_PangenomeSNP) {
        for (const auto& [key, pangenomeSNP_Pointer] : *key_To_PangenomeSNP_Pointer) {
            all_unique_Pangenome_SNPs_as_raw_pointers.insert(pangenomeSNP_Pointer.get());
        }
    }
    return all_unique_Pangenome_SNPs_as_raw_pointers;
}


void PangenomeSNPsIndex::add_Positioned_SNP_from_genomic_coordinates(const GenomicCoordinate &genomic_coordinate_1,
                                                 const GenomicCoordinate &genomic_coordinate_2,
                                                 const char allele_1, const char allele_2) {
    PangenomeSNPsIndexKeyType key_of_snp_1(genomic_coordinate_1, allele_1, allele_2);
    PangenomeSNPsIndexKeyType key_of_snp_2(genomic_coordinate_2, allele_1, allele_2);
    make_the_two_keys_point_to_the_same_pangenome_SNP(key_of_snp_1, key_of_snp_2, allele_1, allele_2);
}



void PangenomeSNPsIndex::add_Positioned_SNP(
        const std::string &genome_1,
        const std::string &chrom_1,
        const uint32_t pos_1,
        const std::string &genome_2,
        const std::string &chrom_2,
        const uint32_t pos_2,
        const char allele_1,
        const char allele_2) {
    GenomicCoordinate genomic_coordinate_1 (genome_1, chrom_1, pos_1);
    GenomicCoordinate genomic_coordinate_2 (genome_2, chrom_2, pos_2);
    add_Positioned_SNP_from_genomic_coordinates(genomic_coordinate_1, genomic_coordinate_2, allele_1, allele_2);
}

std::map<std::string, uint32_t> PangenomeSNPsIndex::get_nb_SNPs_that_can_be_found_with_the_given_genomes (const std::vector<std::string> &genomes) const {
    /**
     * @return: A map with the genome names and the nb of SNPs that can be found in each. Additionally, an entry with "all" is added, with the nb of SNPs in the pangenome
     */
    std::map<std::string, uint32_t> genomes_to_nb_of_SNPs;
    std::set<PangenomeSNP_Raw_Pointer> all_unique_Pangenome_SNPs_as_raw_pointers = get_all_unique_Pangenome_SNPs_as_raw_pointers();

    for (const PangenomeSNP_Raw_Pointer &pangenomeSNP : all_unique_Pangenome_SNPs_as_raw_pointers) {
        std::set<std::string> all_genomes_this_pangenome_snp_is_present = pangenomeSNP->get_all_genomes_this_pangenome_snp_is_present();
        for (const std::string &genome : all_genomes_this_pangenome_snp_is_present) {
            bool genome_should_be_counted = std::find(genomes.begin(), genomes.end(), genome) != genomes.end();
            if (genome_should_be_counted)
                genomes_to_nb_of_SNPs[genome]++;
        }
    }
    genomes_to_nb_of_SNPs["all"] = all_unique_Pangenome_SNPs_as_raw_pointers.size();

    return genomes_to_nb_of_SNPs;
}


PangenomeSNPsIndex PangenomeSNPsIndex::load (const std::string &filename) {
    std::ifstream ifs(filename);
    boost::archive::text_iarchive ia(ifs);
    PangenomeSNPsIndex pangenomeSNPsIndex;
    ia >> pangenomeSNPsIndex;
    return pangenomeSNPsIndex;
}

PangenomeSNPsIndex PangenomeSNPsIndex::save (const std::string &filename) const {
    std::ofstream ofs(filename);
    boost::archive::text_oarchive oa(ofs);
    oa << *this;
}