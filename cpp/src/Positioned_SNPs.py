from collections import namedtuple
import pickle

Position = namedtuple('Position', ['genome', 'chrom', 'pos'])
class PositionedSNP:
    '''
    Represents a positioned SNP: two alleles and several positions in the genomes where this SNP appears
    '''
    def __init__(self, allele_1, allele_2):
        assert allele_1 < allele_2, "Alleles should be given in a canonical order"
        self.alleles = (allele_1, allele_2)
        self.positions = set()

    def add_pos(self, position):
        self.positions.add(position)

    def merge(self, other):
        '''
        Creates a new PositionedSNP with the same allele and positions merged
        :param other: other PositionedSNP
        :return: new PositionedSNP with the same allele and positions merged
        '''
        assert self.alleles == other.alleles  # we should not merge SNPs that does not have the same alleles
        new_PositionedSNP = PositionedSNP(*self.alleles)
        new_PositionedSNP.positions = self.positions.union(other.positions)
        return new_PositionedSNP

    @classmethod
    def from_dict(cls, dict):
        '''
        Factory method mainly used for testing
        :param dict: represents a dictionary to configure a new PositionedSNP object
        :return: new PositionedSNP object
        '''
        new_positioned_SNP = PositionedSNP("A", "C")
        new_positioned_SNP.__dict__ = dict
        return new_positioned_SNP

    def __eq__(self, other):
        """Overrides the default implementation"""
        if isinstance(other, PositionedSNP):
            return self.__dict__ == other.__dict__
        elif isinstance(other, dict):
            return self.__dict__ == other
        return NotImplemented

    def __repr__(self):
        return str(self.__dict__)


PositionedSNPIndexKeyType = namedtuple('PositionedSNPIndexKeyType', ['position', 'allele_1', 'allele_2'])
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

    def add_PositionedSNP(self, position_1, position_2, allele_1, allele_2):
        '''
        Add the SNP represented by the positions and alleles to this index
        :param position_1: position 1 of the SNP
        :param position_2: position 2 of the SNP
        :param allele_1: first base
        :param allele_2: second base
        '''
        # create the keys to be added
        PositionedSNPIndexKey_1 = PositionedSNPIndexKeyType(position=position_1, allele_1=allele_1, allele_2=allele_2)
        PositionedSNPIndexKey_2 = PositionedSNPIndexKeyType(position=position_2, allele_1=allele_1, allele_2=allele_2)

        # get the previous positioned SNPs, if any
        previous_PositionedSNP_1 = self.__get_PositionedSNP(PositionedSNPIndexKey_1)
        previous_PositionedSNP_2 = self.__get_PositionedSNP(PositionedSNPIndexKey_2)

        current_PositionedSNP = None
        if previous_PositionedSNP_1 is None and previous_PositionedSNP_2 is None:
            # we need to create a new PositionedSNP
            current_PositionedSNP = PositionedSNP(allele_1, allele_2)

            # and associate it to these positions in the index
            self.__set_PositionedSNPIndexKey_to_PositionedSNP(PositionedSNPIndexKey_1, current_PositionedSNP)
            self.__set_PositionedSNPIndexKey_to_PositionedSNP(PositionedSNPIndexKey_2, current_PositionedSNP)
        elif previous_PositionedSNP_1 is None:
            # update self.PositionedSNPIndexKey_to_PositionedSNP[PositionedSNPIndexKey_1]
            current_PositionedSNP = previous_PositionedSNP_2
            self.__set_PositionedSNPIndexKey_to_PositionedSNP(PositionedSNPIndexKey_1, current_PositionedSNP)
        elif previous_PositionedSNP_2 is None:
            # update self.PositionedSNPIndexKey_to_PositionedSNP[PositionedSNPIndexKey_2]
            current_PositionedSNP = previous_PositionedSNP_1
            self.__set_PositionedSNPIndexKey_to_PositionedSNP(PositionedSNPIndexKey_2, current_PositionedSNP)
        elif not (previous_PositionedSNP_1 is previous_PositionedSNP_2):
            # if both SNPs do not point to the same PositionedSNP, it means they have to be merged
            current_PositionedSNP = previous_PositionedSNP_1.merge(previous_PositionedSNP_2)

            # we associate these positions in the index to the new merged PositionedSNP
            self.__set_PositionedSNPIndexKey_to_PositionedSNP(PositionedSNPIndexKey_1, current_PositionedSNP)
            self.__set_PositionedSNPIndexKey_to_PositionedSNP(PositionedSNPIndexKey_2, current_PositionedSNP)

            # and also all the previous positions now point to this PositionedSNP
            for position in current_PositionedSNP.positions:
                self.__set_PositionedSNPIndexKey_to_PositionedSNP(
                    PositionedSNPIndexKeyType(position=position, allele_1=allele_1, allele_2=allele_2),
                    current_PositionedSNP)
        else:
            assert previous_PositionedSNP_1 is previous_PositionedSNP_2 and previous_PositionedSNP_1 is not None and previous_PositionedSNP_2 is not None
            current_PositionedSNP = previous_PositionedSNP_1

        # add the positions to current_PositionedSNP
        assert current_PositionedSNP is not None
        current_PositionedSNP.add_pos(position_1)
        current_PositionedSNP.add_pos(position_2)

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