#include <string>
#include <utility>
#include <set>
#include <memory>


class GenomicCoordinate {
public:
    const std::string genome;
    const std::string chrom;
    const uint32_t pos;
    GenomicCoordinate (const std::string &genome, const std::string &chrom, uint32_t pos) :
        genome(genome), chrom(chrom), pos(pos) {}
};

class PangenomeSNP {
    // Represents a pangenome SNP: two alleles and several positions in the genomes where this SNP appears
private:
    const char allele_1;
    const char allele_2;
    std::set<GenomicCoordinate> genomic_coordinates;
public:
    PangenomeSNP(const char allele_1, const char allele_2) {
        assert(allele_1 < allele_2);
        this->allele_1 = allele_1;
        this->allele_2 = allele_2;
    }

    void add_genomic_coordinate (const GenomicCoordinate &genomic_coordinate) {
        genomic_coordinates.insert(genomic_coordinate);
    }

    PangenomeSNP merge(const PangenomeSNP &other) const {
        /*
        Creates a new PangenomeSNP with the same allele and positions merged
        :param other: other PangenomeSNP
        :return: new PangenomeSNP with the same allele and positions merged
         */
        assert(this->allele_1 == other.allele_1 && this->allele_2 == other.allele_2); // we should not merge PangenomeSNPs that do not have the same alleles
        new_PangenomeSNP = PangenomeSNP(this->allele_1, this->allele_2);
        new_PangenomeSNP.genomic_coordinates.insert(this->genomic_coordinates.begin(), this->genomic_coordinates.end());
        new_PangenomeSNP.genomic_coordinates.insert(other.genomic_coordinates.begin(), other.genomic_coordinates.end());
        return new_PangenomeSNP;
    }

    bool operator< (const PangenomeSNP &other) const {
        if (this->allele_1 != other.allele_1) return this->allele_1 < other.allele_1;
        if (this->allele_2 != other.allele_2) return this->allele_2 < other.allele_2;
        if (this->genomic_coordinates != other.genomic_coordinates) return this->genomic_coordinates < other.genomic_coordinates;
        return false;
    }
};


/**
 * Class used to be the key in the PangenomeSNPsIndex
 */
class PangenomeSNPsIndexKeyType {
public:
    const GenomicCoordinate genomic_coordinate;
    const char allele_1;
    const char allele_2;
    PangenomeSNPsIndexKeyType (const GenomicCoordinate &genomic_coordinate, const char allele_1, const char allele_2) :
            genomic_coordinate(genomic_coordinate), allele_1(allele_1), allele_2(allele_2) {}
    bool operator< (const Positioned_SNP &other) const {
        if (this->allele_1 != other.allele_1) return this->allele_1 < other.allele_1;
        if (this->allele_2 != other.allele_2) return this->allele_2 < other.allele_2;
        if (this->genomic_coordinate != other.genomic_coordinate) return this->genomic_coordinate < other.genomic_coordinate;
        return false;
    }
};


class PangenomeSNPsIndex {
    // Represents an index of PositionedSNPs, where the key is PangenomeSNPsIndexKeyType and value is PangenomeSNP
private:
    using PangenomeSNPsIndexKey_To_PangenomeSNP = std::map<PangenomeSNPsIndexKeyType, PangenomeSNP>;
    using Position_To_PangenomeSNPsIndexKey_To_PangenomeSNP =
            std::vector<std::shared_pointer<PangenomeSNPsIndexKey_To_PangenomeSNP>>;

    PangenomeSNP* get_PangenomeSNP(const PangenomeSNPsIndexKeyType &pangenomeSNPsIndexKey) {
        /**
         * :return: the PangenomeSNP on the coordinates given by pangenomeSNPsIndexKey, or None if it does not exist
         */
        std::shared_pointer<PangenomeSNPsIndexKey_To_PangenomeSNP> pangenomeSNPsIndexKey_To_PangenomeSNP =
                position_To_PangenomeSNPsIndexKey_To_PangenomeSNP[pangenomeSNPsIndexKey.genomic_coordinate.pos];

        if (pangenomeSNPsIndexKey_To_PangenomeSNP != nullptr)
            return pangenomeSNPsIndexKey_To_PangenomeSNP->at(pangenomeSNPsIndexKey);
        else
            return nullptr;
    }


    void set_PositionedSNP(const PangenomeSNPsIndexKeyType &pangenomeSNPsIndexKey, const PangenomeSNP &pangenomeSNP) {
        /**
         * set the pangenomeSNP to the pangenomeSNPsIndexKey
         */
        std::shared_pointer<PangenomeSNPsIndexKey_To_PangenomeSNP> &pangenomeSNPsIndexKey_To_PangenomeSNP =
                position_To_PangenomeSNPsIndexKey_To_PangenomeSNP[pangenomeSNPsIndexKey.genomic_coordinate.pos];
        std::shared_pointer<PositionToPositionedSNPIndexKeyToPositionedSNPMap> &positionedSNPIndexKeyToPositionedSNPMap =
                positionToPositionedSNPIndexKeyToPositionedSNPMap[positionedSNPIndexKey.position.pos];
        if (positionedSNPIndexKeyToPositionedSNPMap == nullptr) {
            positionedSNPIndexKeyToPositionedSNPMap = make_shared<PositionedSNPIndexKeyToPositionedSNPMap>();
        }
        positionedSNPIndexKeyToPositionedSNPMap[positionedSNPIndexKey] = positioned_SNP;
    }



public:
    PositionedSNPsIndex (uint32_t length_of_longest_contig) {
        Position_To_PangenomeSNPsIndexKey_To_PangenomeSNP position_To_PangenomeSNPsIndexKey_To_PangenomeSNP(length_of_longest_contig);
    }

    void add_Positioned_SNP(
            const std::string &genome_1,
            const std::string &chrom_1,
            const uint32_t pos_1,
            const std::string &genome_2,
            const std::string &chrom_2,
            const uint32_t pos_2,
            const char allele_1,
            const char allele_2) {
        //Add the SNP represented by the positions and alleles to this index

        Position position_1(genome_1, chrom_1, pos_1);
        Position position_2(genome_2, chrom_2, pos_2);
        add_Positioned_SNP(position_1, position_2, allele_1, allele_2);
    }

    void add_Positioned_SNP(const Position &position1, const Position &position2, const char allele_1, const char allele_2) {
        /*
         *     Add the SNP represented by the positions and alleles to this index
                    :param position_1: position 1 of the SNP
                    :param position_2: position 2 of the SNP
                    :param allele_1: first base
                    :param allele_2: second base

         */
        // create the keys to be added
        PositionedSNPIndexKeyType PositionedSNPIndexKey_1(position_1, allele_1, allele_2);
        PositionedSNPIndexKeyType PositionedSNPIndexKey_2(position_2, allele_1, allele_2);

        // get the previous positioned SNPs, if any
        std::shared_pointer<Positioned_SNP> previous_PositionedSNP_1 = get_PositionedSNP(PositionedSNPIndexKey_1);
        std::shared_pointer<Positioned_SNP> previous_PositionedSNP_2 = get_PositionedSNP(PositionedSNPIndexKey_2);

        //denotes the Positioned_SNP to be updated
        std::shared_pointer<Positioned_SNP> current_PositionedSNP = nullptr;

        if (previous_PositionedSNP_1 == nullptr && previous_PositionedSNP_2 == nullptr) {
            // we need to create a new PositionedSNP
            current_PositionedSNP = std::make_shared<PositionedSNP>(allele_1, allele_2);

            // and associate it to these positions in the index
            set_PositionedSNPIndexKey_to_PositionedSNP(PositionedSNPIndexKey_1, current_PositionedSNP);
            set_PositionedSNPIndexKey_to_PositionedSNP(PositionedSNPIndexKey_2, current_PositionedSNP);
        }
        else if (previous_PositionedSNP_1 == nullptr) {
            // update PositionedSNPIndexKey_to_PositionedSNP[PositionedSNPIndexKey_1]
            current_PositionedSNP = previous_PositionedSNP_2;
            set_PositionedSNPIndexKey_to_PositionedSNP(PositionedSNPIndexKey_1, current_PositionedSNP);
        }
        else if (previous_PositionedSNP_2 == nullptr) {
            // update self.PositionedSNPIndexKey_to_PositionedSNP[PositionedSNPIndexKey_2]
            current_PositionedSNP = previous_PositionedSNP_1;
            set_PositionedSNPIndexKey_to_PositionedSNP(PositionedSNPIndexKey_2, current_PositionedSNP);
        }
        else if (previous_PositionedSNP_1 != previous_PositionedSNP_2) {
            // if both SNPs do not point to the same PositionedSNP, it means they have to be merged
            current_PositionedSNP = std::make_shared(previous_PositionedSNP_1->merge(previous_PositionedSNP_2));

            // we associate these positions in the index to the new merged PositionedSNP
            set_PositionedSNPIndexKey_to_PositionedSNP(PositionedSNPIndexKey_1, current_PositionedSNP);
            set_PositionedSNPIndexKey_to_PositionedSNP(PositionedSNPIndexKey_2, current_PositionedSNP);

            // and also all the previous positions now point to this PositionedSNP
            for (const Position &position : current_PositionedSNP->positions)
                set_PositionedSNPIndexKey_to_PositionedSNP(position, allele_1, allele_2,current_PositionedSNP);
        }
        else {
            assert(previous_PositionedSNP_1 == previous_PositionedSNP_2 && previous_PositionedSNP_1 != nullptr && previous_PositionedSNP_2 != nullptr)
            current_PositionedSNP = previous_PositionedSNP_1
        }

        // add the positions to current_PositionedSNP
        assert(current_PositionedSNP != nullptr)
        current_PositionedSNP->add_pos(position_1)
        current_PositionedSNP->add_pos(position_2)
    }


};
class PositionedSNPsIndex:
    def __get_PositionedSNP(self, positionedSNPIndexKey):
        '''
        :return: the Positioned SNP on the coordinates given by PositionedSNPIndexKey, or None if it does not exist
        '''
        dict_in_pos = self._PositionedSNPIndexKey_to_PositionedSNP[positionedSNPIndexKey.position.pos]
        if dict_in_pos is not None:
            return dict_in_pos.get((positionedSNPIndexKey.position.genome, positionedSNPIndexKey.position.chrom,
                                    positionedSNPIndexKey.allele_1, positionedSNPIndexKey.allele_2))
        return None

    def __set_PositionedSNPIndexKey_to_PositionedSNP(self, positionedSNPIndexKey, positionedSNP):
        '''
        set the positionedSNP to the positionedSNPIndexKey
        '''
        if self._PositionedSNPIndexKey_to_PositionedSNP[positionedSNPIndexKey.position.pos] is None:
            self._PositionedSNPIndexKey_to_PositionedSNP[positionedSNPIndexKey.position.pos] = {}
        dict_in_pos = self._PositionedSNPIndexKey_to_PositionedSNP[positionedSNPIndexKey.position.pos]
        dict_in_pos[(positionedSNPIndexKey.position.genome, positionedSNPIndexKey.position.chrom,
                         positionedSNPIndexKey.allele_1, positionedSNPIndexKey.allele_2)] = positionedSNP

    '''
    Represents an index of PositionedSNPs, where the key is PositionedSNPIndexKeyType and value is PositionedSNP
    '''
    def __init__(self, length_of_longest_contig):
        self._PositionedSNPIndexKey_to_PositionedSNP = [None] * length_of_longest_contig  # index the positioned SNPs


    def add_SNPs_from_csv(self, csv_file, genome_1, genome_2):
        '''
        :param csv_file: a csv file with the SNPs computed by get_SNPs_using_mummer rule
        :param genome_1: a string with genome_1 name
        :param genome_2: a string with genome_2 name
        '''

        print(f"[DEBUG_add_SNPs_from_csv]: add_SNPs_from_csv('{csv_file}', '{genome_1}', '{genome_2}')")

        '''
        # nice engineered solution, but too slow
        snps_dataframe = pd.read_csv(csv_file, sep = "\t")

        # populate self.PositionedSNPIndexKey_to_PositionedSNP
        for i, row in snps_dataframe.iterrows():
            # get the data
            position_1 = Position(genome=genome_1, chrom=row["ref_chrom"], pos=row["ref_pos"])
            position_2 = Position(genome=genome_2, chrom=row["query_chrom"], pos=row["query_pos"])
            allele_1 = min(row["ref_sub"], row["query_sub"])
            allele_2 = max(row["ref_sub"], row["query_sub"])
            self.add_PositionedSNP(position_1, position_2, allele_1, allele_2)
        '''

        # crude solution, but fast
        with open(csv_file) as fin:
            all_csv_lines = fin.readlines()

        # populate self.PositionedSNPIndexKey_to_PositionedSNP
        for i, line in enumerate(all_csv_lines):
            if i == 0: continue  # skip header

            # get the data
            line_split = line.split()
            position_1 = Position(genome=genome_1, chrom=line_split[13], pos=int(line_split[1]))
            position_2 = Position(genome=genome_2, chrom=line_split[14], pos=int(line_split[4]))
            ref_base, query_base = line_split[2], line_split[3]
            allele_1 = min(ref_base, query_base)
            allele_2 = max(ref_base, query_base)
            self.add_PositionedSNP(position_1, position_2, allele_1, allele_2)




    def get_all_unique_Positioned_SNPs(self):
        # as there are no Positioned_SNPs deep copies, their memory address is the unique identifier
        all_unique_positioned_SNPs = []
        all_unique_positioned_SNPs_ids = set()
        for dict_in_pos in self._PositionedSNPIndexKey_to_PositionedSNP:
            if dict_in_pos is not None:
                for positioned_snp in dict_in_pos.values():
                    if id(positioned_snp) not in all_unique_positioned_SNPs_ids:
                        all_unique_positioned_SNPs_ids.add(id(positioned_snp))
                        all_unique_positioned_SNPs.append(positioned_snp)
        return all_unique_positioned_SNPs

    def get_nb_SNPs_that_can_be_found_with_the_given_genomes(self, genomes):
        """
        :param genomes: the genomes to be found (list of string)
        :return: A dictionary with the genome names and the nb of SNPs that can be found in each. Additionally, an entry with "all" is added, with the nb of SNPs in the pangenome
        """
        genomes_to_nb_of_SNPs = {}
        all_unique_positioned_SNPs = self.get_all_unique_Positioned_SNPs()

        genomes_to_nb_of_SNPs["all"] = len(all_unique_positioned_SNPs)

        for genome in genomes:
            nb_of_SNPs = 0
            for positioned_SNP in all_unique_positioned_SNPs:
                for position in positioned_SNP.positions:
                    if position.genome == genome:
                        nb_of_SNPs += 1
                        break
            genomes_to_nb_of_SNPs[genome] = nb_of_SNPs

        return genomes_to_nb_of_SNPs

    @classmethod
    def load(cls, filename):
        with open(filename, "rb") as fin:
            positionedSNPs = pickle.load(fin)
        return positionedSNPs

    def save(self, filename):
        with open(filename, "wb") as fout:
            pickle.dump(self, fout)

    def to_dict(self):
        '''
        Transforms self._PositionedSNPIndexKey_to_PositionedSNP into a dictionary, mainly used for testing
        '''
        self_as_dict = {}
        for pos, dict_in_pos in enumerate(self._PositionedSNPIndexKey_to_PositionedSNP):
            if dict_in_pos is not None:
                self_as_dict[pos] = dict_in_pos
        return self_as_dict

    def __eq__(self, other):
        """Overrides the default implementation"""
        if isinstance(other, PositionedSNPsIndex):
            return self.__dict__ == other.__dict__
        return NotImplemented


    def __repr__(self):
        return str(self.__dict__)


'''
def test():
    positionedSNPsIndex = PositionedSNPsIndex(5572075)
    positionedSNPsIndex.add_SNPs_from_csv('assemblies_sample_out/genome.0/genome.0-SEP-genome.1.mummer.csv', 'genome.0', 'genome.1')
    positionedSNPsIndex.add_SNPs_from_csv('assemblies_sample_out/genome.0/genome.0-SEP-genome.2.mummer.csv', 'genome.0', 'genome.2')
    positionedSNPsIndex.add_SNPs_from_csv('assemblies_sample_out/genome.0/genome.0-SEP-genome.3.mummer.csv', 'genome.0', 'genome.3')
    positionedSNPsIndex.add_SNPs_from_csv('assemblies_sample_out/genome.1/genome.1-SEP-genome.2.mummer.csv', 'genome.1', 'genome.2')
    positionedSNPsIndex.add_SNPs_from_csv('assemblies_sample_out/genome.1/genome.1-SEP-genome.3.mummer.csv', 'genome.1', 'genome.3')
    positionedSNPsIndex.add_SNPs_from_csv('assemblies_sample_out/genome.2/genome.2-SEP-genome.3.mummer.csv', 'genome.2', 'genome.3')

test()
'''

'''
import cProfile
cProfile.run("test()", 'profiling_stats')
'''