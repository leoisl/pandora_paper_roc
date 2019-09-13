#ifndef PANDORA1_PAPER_PANGENOME_SNPS_H
#define PANDORA1_PAPER_PANGENOME_SNPS_H

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
private:
    std::string genome;
    std::string chrom;
    uint32_t pos{};

public:
    GenomicCoordinate() = default;

    GenomicCoordinate (const std::string &genome, const std::string &chrom, uint32_t pos) :
        genome(genome), chrom(chrom), pos(pos) {}


    bool operator<(const GenomicCoordinate &rhs) const {
        return std::tie(genome, chrom, pos) < std::tie(rhs.genome, rhs.chrom, rhs.pos);
    }

    bool operator>(const GenomicCoordinate &rhs) const {
        return rhs < *this;
    }

    bool operator<=(const GenomicCoordinate &rhs) const {
        return !(rhs < *this);
    }

    bool operator>=(const GenomicCoordinate &rhs) const {
        return !(*this < rhs);
    }

    bool operator==(const GenomicCoordinate &rhs) const {
        return std::tie(genome, chrom, pos) == std::tie(rhs.genome, rhs.chrom, rhs.pos);
    }

    bool operator!=(const GenomicCoordinate &rhs) const {
        return !(rhs == *this);
    }

    const std::string &get_genome() const {
        return genome;
    }

    const std::string &get_chrom() const {
        return chrom;
    }

    uint32_t get_pos() const {
        return pos;
    }

private:
    friend class boost::serialization::access;
    template<class Archive>
    inline void serialize(Archive &ar, const unsigned int version) {
        ar & genome;
        ar & chrom;
        ar & pos;
    }
};


// Represents a pangenome SNP: two alleles and several positions in the genomes where this SNP appears
class PangenomeSNP {
private:
    char allele_1{};
    char allele_2{};
    std::set<GenomicCoordinate> genomic_coordinates;

public:
    PangenomeSNP()= default;

    PangenomeSNP(const char allele_1, const char allele_2) : allele_1(allele_1), allele_2(allele_2), genomic_coordinates() {
        bool alleles_are_ordered = allele_1 < allele_2;
        assert(alleles_are_ordered);
    }

    void add_genomic_coordinate (const GenomicCoordinate &genomic_coordinate) {
        genomic_coordinates.insert(genomic_coordinate);
    }

    // Creates a new PangenomeSNP with the same alleles and positions merged
    PangenomeSNP merge(const PangenomeSNP &other) const;

    const std::set<GenomicCoordinate>& get_genomic_coordinates () const {
        return this->genomic_coordinates;
    }

    char get_allele_1() const {
        return allele_1;
    }

    char get_allele_2() const {
        return allele_2;
    }

    std::set<std::string> get_all_genomes_this_pangenome_snp_is_present() const;




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

    bool operator==(const PangenomeSNP &rhs) const {
        return std::tie(allele_1, allele_2, genomic_coordinates) ==
               std::tie(rhs.allele_1, rhs.allele_2, rhs.genomic_coordinates);
    }

    bool operator!=(const PangenomeSNP &rhs) const {
        return !(rhs == *this);
    }

private:
    friend class boost::serialization::access;
    template<class Archive>
    inline void serialize(Archive &ar, const unsigned int version) {
        ar & allele_1;
        ar & allele_2;
        ar & genomic_coordinates;
    }

};
using PangenomeSNP_Pointer = boost::shared_ptr<PangenomeSNP>;
using PangenomeSNP_Raw_Pointer = const PangenomeSNP *;


// Represents the key type of the index
class PangenomeSNPsIndexKeyType {
private:
    GenomicCoordinate genomic_coordinate;
    char allele_1;
    char allele_2;

public:
    PangenomeSNPsIndexKeyType(){}

    PangenomeSNPsIndexKeyType (const GenomicCoordinate &genomic_coordinate, const char allele_1, const char allele_2) :
            genomic_coordinate(genomic_coordinate), allele_1(allele_1), allele_2(allele_2) {}

    const GenomicCoordinate &get_genomic_coordinate() const {
        return genomic_coordinate;
    }

    char get_allele_1() const {
        return allele_1;
    }

    char get_allele_2() const {
        return allele_2;
    }

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

    bool operator==(const PangenomeSNPsIndexKeyType &rhs) const {
        return std::tie(genomic_coordinate, allele_1, allele_2) ==
               std::tie(rhs.genomic_coordinate, rhs.allele_1, rhs.allele_2);
    }

    bool operator!=(const PangenomeSNPsIndexKeyType &rhs) const {
        return !(rhs == *this);
    }

private:
    friend class boost::serialization::access;
    template<class Archive>
    inline void serialize(Archive &ar, const unsigned int version) {
        ar & genomic_coordinate;
        ar & allele_1;
        ar & allele_2;
    }

};
using Key_To_PangenomeSNP = std::map<PangenomeSNPsIndexKeyType, PangenomeSNP_Pointer>;
using Key_To_PangenomeSNP_Pointer = boost::shared_ptr<Key_To_PangenomeSNP>;


// TODO: extract methods to an abstract PangenomeSNPsMap class and make this an specialization of this class
class PangenomeSNPsMapWithFastPositionAccess {
private:
    using Position_To_Key_To_PangenomeSNP = std::vector<Key_To_PangenomeSNP_Pointer>;
    Position_To_Key_To_PangenomeSNP position_To_Key_To_PangenomeSNP;

    bool position_does_not_exist(uint32_t position) const {
        return position >= position_To_Key_To_PangenomeSNP.size();
    }

public:
    PangenomeSNPsMapWithFastPositionAccess() = default;
    explicit PangenomeSNPsMapWithFastPositionAccess (uint32_t number_of_indexed_positions) : position_To_Key_To_PangenomeSNP(number_of_indexed_positions) {}

    PangenomeSNP_Pointer get_PangenomeSNP(const PangenomeSNPsIndexKeyType &key) const;

    void set_PangenomeSNP(const PangenomeSNPsIndexKeyType &key, const PangenomeSNP_Pointer &pangenomeSNP);

    PangenomeSNP_Pointer merge_two_PangenomeSNPs(const PangenomeSNP_Pointer &previous_PangenomeSNP_in_key_of_snp_1,
                                                 const PangenomeSNP_Pointer &previous_PangenomeSNP_in_key_of_snp_2);

    std::set<PangenomeSNP_Raw_Pointer> get_all_unique_Pangenome_SNPs_as_raw_pointers() const;

    void make_the_two_keys_point_to_the_same_pangenome_SNP (
            const PangenomeSNPsIndexKeyType &key_of_snp_1,
            const PangenomeSNPsIndexKeyType &key_of_snp_2,
            const char allele_1, const char allele_2);

private:
    friend class boost::serialization::access;
    template<class Archive>
    inline void serialize(Archive &ar, const unsigned int version) {
        ar & position_To_Key_To_PangenomeSNP;
    }
};



class PangenomeSNPsIndex {
private:
    PangenomeSNPsMapWithFastPositionAccess pangenomeSNPsMap;
public:
    PangenomeSNPsIndex (){}
    PangenomeSNPsIndex (uint32_t index_size) : pangenomeSNPsMap(index_size) {}

    void add_Positioned_SNP_from_genomic_coordinates(const GenomicCoordinate &genomic_coordinate_1,
                                                     const GenomicCoordinate &genomic_coordinate_2,
                                                     const char allele_1, const char allele_2);

    void add_Positioned_SNP(
            const std::string &genome_1,
            const std::string &chrom_1,
            const uint32_t pos_1,
            const std::string &genome_2,
            const std::string &chrom_2,
            const uint32_t pos_2,
            const char allele_1,
            const char allele_2);

    std::map<std::string, uint32_t> get_nb_SNPs_that_can_be_found_with_the_given_genomes (const std::vector<std::string> &genomes) const;


    static PangenomeSNPsIndex load (const std::string &filename);

    PangenomeSNPsIndex save (const std::string &filename) const;
private:
    friend class boost::serialization::access;
    template<class Archive>
    inline void serialize(Archive &ar, const unsigned int version) {
        ar & pangenomeSNPsMap;
    }
};


#endif //PANDORA1_PAPER_PANGENOME_SNPS_H