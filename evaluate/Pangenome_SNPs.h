#include <string>
#include <utility>
#include <set>
#include <cassert>
#include <map>
#include <vector>
#include <algorithm>
#include <fstream>
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#include <boost/serialization/string.hpp>
#include <boost/serialization/set.hpp>
#include <boost/serialization/vector.hpp>
#include <boost/serialization/map.hpp>
#include <boost/serialization/shared_ptr.hpp>
#include <boost/smart_ptr/shared_ptr.hpp>
#include <boost/smart_ptr/make_shared.hpp>


class GenomicCoordinate {
public: // it was a choice to put all attributes as public and const - this class just aggregates info (like a struct)
    const std::string genome;
    const std::string chrom;
    const uint32_t pos;
    GenomicCoordinate (const std::string &genome, const std::string &chrom, uint32_t pos) :
        genome(genome), chrom(chrom), pos(pos) {}


private:
    friend class boost::serialization::access;
    template<class Archive>
    void serialize(Archive & ar, const unsigned int version)
    {
        ar & genome;
        ar & chrom;
        ar & pos;
    }
};


// Represents a pangenome SNP: two alleles and several positions in the genomes where this SNP appears
class PangenomeSNP {
private:
    const char allele_1;
    const char allele_2;
    std::set<GenomicCoordinate> genomic_coordinates;

    friend class boost::serialization::access;
    template<class Archive>
    void serialize(Archive & ar, const unsigned int version)
    {
        ar & allele_1;
        ar & allele_2;
        ar & genomic_coordinates;
    }
public:
    PangenomeSNP(const char allele_1, const char allele_2) : allele_1(allele_1), allele_2(allele_2), genomic_coordinates() {
        assert(allele_1 < allele_2);
    }

    void add_genomic_coordinate (const GenomicCoordinate &genomic_coordinate) {
        genomic_coordinates.insert(genomic_coordinate);
    }

    // Creates a new PangenomeSNP with the same alleles and positions merged
    PangenomeSNP merge(const PangenomeSNP &other) const {
        assert(this->allele_1 == other.allele_1 && this->allele_2 == other.allele_2); // we should not merge PangenomeSNPs that do not have the same alleles
        PangenomeSNP new_PangenomeSNP(this->allele_1, this->allele_2);
        new_PangenomeSNP.genomic_coordinates.insert(this->genomic_coordinates.begin(), this->genomic_coordinates.end());
        new_PangenomeSNP.genomic_coordinates.insert(other.genomic_coordinates.begin(), other.genomic_coordinates.end());
        return new_PangenomeSNP;
    }

    const std::set<GenomicCoordinate>& get_genomic_coordinates () const {
        return this->genomic_coordinates;
    }

    std::set<std::string> get_all_genomes_this_pangenome_snp_is_present() const {
        std::set<std::string> all_genomes_this_pangenome_snp_is_present;
        for (const GenomicCoordinate &genomic_coordinate : genomic_coordinates) {
            all_genomes_this_pangenome_snp_is_present.insert(genomic_coordinate.genome);
        }
        return all_genomes_this_pangenome_snp_is_present;
    }

    bool operator<(const PangenomeSNP &rhs) const {
        return std::tie(allele_1, allele_2, genomic_coordinates) <
               std::tie(rhs.allele_1, rhs.allele_2, rhs.genomic_coordinates);
    }

    bool operator>(const PangenomeSNP &rhs) const {
        return rhs < *this;
    }

    bool operator<=(const PangenomeSNP &rhs) const {
        return !(rhs < *this);
    }

    bool operator>=(const PangenomeSNP &rhs) const {
        return !(*this < rhs);
    }
};


class PangenomeSNPsIndex {
private:
    // Represents the key type of the index
    class PangenomeSNPsIndexKeyType {
    private:
        friend class boost::serialization::access;
        template<class Archive>
        void serialize(Archive & ar, const unsigned int version)
        {
            ar & genomic_coordinate;
            ar & allele_1;
            ar & allele_2;
        }

    public:
        const GenomicCoordinate genomic_coordinate;
        const char allele_1;
        const char allele_2;

        PangenomeSNPsIndexKeyType (const GenomicCoordinate &genomic_coordinate, const char allele_1, const char allele_2) :
                genomic_coordinate(genomic_coordinate), allele_1(allele_1), allele_2(allele_2) {}

        bool operator<(const PangenomeSNPsIndexKeyType &rhs) const {
            return std::tie(genomic_coordinate, allele_1, allele_2) <
                   std::tie(rhs.genomic_coordinate, rhs.allele_1, rhs.allele_2);
        }

        bool operator>(const PangenomeSNPsIndexKeyType &rhs) const {
            return rhs < *this;
        }

        bool operator<=(const PangenomeSNPsIndexKeyType &rhs) const {
            return !(rhs < *this);
        }

        bool operator>=(const PangenomeSNPsIndexKeyType &rhs) const {
            return !(*this < rhs);
        }
    };


    // represents the index
    using PangenomeSNP_Pointer = boost::shared_ptr<PangenomeSNP>;
    using PangenomeSNP_Raw_Pointer = const PangenomeSNP *;
    using Key_To_PangenomeSNP = std::map<PangenomeSNPsIndexKeyType, PangenomeSNP_Pointer>;
    using Key_To_PangenomeSNP_Pointer = boost::shared_ptr<Key_To_PangenomeSNP>;
    using Position_To_Key_To_PangenomeSNP =
            std::vector<Key_To_PangenomeSNP_Pointer>;
    Position_To_Key_To_PangenomeSNP position_To_Key_To_PangenomeSNP;


    PangenomeSNP_Pointer get_PangenomeSNP(const PangenomeSNPsIndexKeyType &key) const {
         /** @return: the PangenomeSNP on the coordinates given by pangenomeSNPsIndexKey, or nullptr if it does not exist */
        uint32_t position = key.genomic_coordinate.pos;
        Key_To_PangenomeSNP_Pointer key_To_PangenomeSNP_Pointer =
                position_To_Key_To_PangenomeSNP[position];

        if (key_To_PangenomeSNP_Pointer != nullptr)
            return key_To_PangenomeSNP_Pointer->at(key);
        else
            return nullptr;
    }


    void set_PangenomeSNP(const PangenomeSNPsIndexKeyType &key, const PangenomeSNP_Pointer &pangenomeSNP) {
        uint32_t position = key.genomic_coordinate.pos;
        Key_To_PangenomeSNP_Pointer key_To_PangenomeSNP_Pointer =
                position_To_Key_To_PangenomeSNP[position];

        if (key_To_PangenomeSNP_Pointer == nullptr) {
            key_To_PangenomeSNP_Pointer.reset(new Key_To_PangenomeSNP());
        }
        (*key_To_PangenomeSNP_Pointer)[key] = pangenomeSNP;
    }


    PangenomeSNP_Pointer merge_two_PangenomeSNPs(const PangenomeSNP_Pointer &previous_PangenomeSNP_in_key_of_snp_1,
                                 const PangenomeSNP_Pointer &previous_PangenomeSNP_in_key_of_snp_2,
                                 const char allele_1, const char allele_2) {
        PangenomeSNP_Pointer merged_SNP = boost::make_shared<PangenomeSNP>(previous_PangenomeSNP_in_key_of_snp_1->merge(*previous_PangenomeSNP_in_key_of_snp_2));

        for (const GenomicCoordinate &genomic_coordinate : merged_SNP->get_genomic_coordinates())
            set_PangenomeSNP(PangenomeSNPsIndexKeyType(genomic_coordinate, allele_1, allele_2), merged_SNP);

        return merged_SNP;
    }

    //TODO: there is just too much logic in here, maybe a candidate for refactoring
    void make_the_two_keys_point_to_the_same_pangenome_SNP (
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
                    pangenomeSNP_pointed_by_both_keys = boost::make_shared(merge_two_PangenomeSNPs(previous_PangenomeSNP_in_key_of_snp_1, previous_PangenomeSNP_in_key_of_snp_2));
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
        pangenomeSNP_pointed_by_both_keys->add_genomic_coordinate(key_of_snp_1.genomic_coordinate);
        pangenomeSNP_pointed_by_both_keys->add_genomic_coordinate(key_of_snp_2.genomic_coordinate);

    }


    std::set<PangenomeSNP_Raw_Pointer> get_all_unique_Pangenome_SNPs_as_raw_pointers() const {
        std::set<PangenomeSNP_Raw_Pointer> all_unique_Pangenome_SNPs_as_raw_pointers;
        for (const auto &key_To_PangenomeSNP_Pointer : position_To_Key_To_PangenomeSNP) {
            for (const auto& [key, pangenomeSNP_Pointer] : *key_To_PangenomeSNP_Pointer) {
                all_unique_Pangenome_SNPs_as_raw_pointers.insert(pangenomeSNP_Pointer.get());
            }
        }
        return all_unique_Pangenome_SNPs_as_raw_pointers;
    }


    friend class boost::serialization::access;
    template<class Archive>
    void serialize(Archive & ar, const unsigned int version)
    {
        ar & position_To_Key_To_PangenomeSNP;
    }

public:
    PangenomeSNPsIndex (uint32_t index_size) : position_To_Key_To_PangenomeSNP(index_size) {}

    void add_Positioned_SNP(
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
        add_Positioned_SNP(genomic_coordinate_1, genomic_coordinate_2, allele_1, allele_2);
    }

    void add_Positioned_SNP(const GenomicCoordinate &genomic_coordinate_1, const GenomicCoordinate &genomic_coordinate_2,
            const char allele_1, const char allele_2) {
        PangenomeSNPsIndexKeyType key_of_snp_1(genomic_coordinate_1, allele_1, allele_2);
        PangenomeSNPsIndexKeyType key_of_snp_2(genomic_coordinate_2, allele_1, allele_2);
        make_the_two_keys_point_to_the_same_pangenome_SNP(key_of_snp_1, key_of_snp_2, allele_1, allele_2);
    }


    std::map<std::string, uint32_t> get_nb_SNPs_that_can_be_found_with_the_given_genomes (const std::vector<std::string> &genomes) const {
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

    PangenomeSNPsIndex load (const std::string &filename) {
        std::ifstream ifs(filename);
        boost::archive::text_iarchive ia(ifs);
        ia >> *this;
    }

    PangenomeSNPsIndex save (const std::string &filename) {
        std::ofstream ofs(filename);
        boost::archive::text_oarchive oa(ofs);
        oa << *this;
    }
};
