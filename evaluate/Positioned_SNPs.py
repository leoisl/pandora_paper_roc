from collections import namedtuple
import pandas as pd
import pickle

Position = namedtuple('Position', ['genome', 'chrom', 'pos'])
class PositionedSNP:
    '''
    Represents a positioned SNP: two alleles and several positions in the genomes where this SNP appears
    '''
    def __init__(self, allele_1, allele_2):
        assert allele_1 < allele_2, "Alleles should be given in a canonical order"
        self.alleles = [allele_1, allele_2]
        self.positions = set()

    def add_pos (self, position):
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
        return NotImplemented

    def __repr__(self):
        return str(self.__dict__)


PositionedSNPIndexKeyType = namedtuple('PositionedSNPIndexKeyType', ['position', 'allele_1', 'allele_2'])
class PositionedSNPsIndex:
    '''
    Represents an index of PositionedSNPs, where the key is PositionedSNPIndexKeyType and value is PositionedSNP
    '''
    def __init__(self):
        self.PositionedSNPIndexKey_to_PositionedSNP = {}  # index the positioned SNPs

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
        previous_PositionedSNP_1 = self.PositionedSNPIndexKey_to_PositionedSNP.get(PositionedSNPIndexKey_1)
        previous_PositionedSNP_2 = self.PositionedSNPIndexKey_to_PositionedSNP.get(PositionedSNPIndexKey_2)

        current_PositionedSNP = None
        if previous_PositionedSNP_1 is None and previous_PositionedSNP_2 is None:
            # we need to create a new PositionedSNP
            current_PositionedSNP = PositionedSNP(allele_1, allele_2)

            # and associate it to these positions in the index
            self.PositionedSNPIndexKey_to_PositionedSNP[PositionedSNPIndexKey_1] = current_PositionedSNP
            self.PositionedSNPIndexKey_to_PositionedSNP[PositionedSNPIndexKey_2] = current_PositionedSNP
        elif previous_PositionedSNP_1 is None:
            # update self.PositionedSNPIndexKey_to_PositionedSNP[PositionedSNPIndexKey_1]
            current_PositionedSNP = previous_PositionedSNP_2
            self.PositionedSNPIndexKey_to_PositionedSNP[PositionedSNPIndexKey_1] = current_PositionedSNP
        elif previous_PositionedSNP_2 is None:
            # update self.PositionedSNPIndexKey_to_PositionedSNP[PositionedSNPIndexKey_2]
            current_PositionedSNP = previous_PositionedSNP_1
            self.PositionedSNPIndexKey_to_PositionedSNP[PositionedSNPIndexKey_2] = current_PositionedSNP
        elif not (previous_PositionedSNP_1 is previous_PositionedSNP_2):
            # if both SNPs do not point to the same PositionedSNP, it means they have to be merged
            current_PositionedSNP = previous_PositionedSNP_1.merge(previous_PositionedSNP_2)

            # we associate these positions in the index to the new merged PositionedSNP
            self.PositionedSNPIndexKey_to_PositionedSNP[PositionedSNPIndexKey_1] = current_PositionedSNP
            self.PositionedSNPIndexKey_to_PositionedSNP[PositionedSNPIndexKey_2] = current_PositionedSNP

            # and also all the previous positions now point to this PositionedSNP
            for position in current_PositionedSNP.positions:
                self.PositionedSNPIndexKey_to_PositionedSNP[PositionedSNPIndexKeyType(position=position, allele_1=allele_1, allele_2=allele_2)] = current_PositionedSNP
        else:
            assert previous_PositionedSNP_1 is previous_PositionedSNP_2 and previous_PositionedSNP_1 is not None and previous_PositionedSNP_2 is not None
            current_PositionedSNP = previous_PositionedSNP_1

        # add the positions to current_PositionedSNP
        assert current_PositionedSNP is not None
        current_PositionedSNP.add_pos(position_1)
        current_PositionedSNP.add_pos(position_2)
        pass

    def add_SNPs_from_csv(self, csv_file, genome_1, genome_2):
        '''
        :param csv_file: a csv file with the SNPs computed by get_SNPs_using_mummer rule
        :param genome_1: a string with genome_1 name
        :param genome_2: a string with genome_2 name
        '''
        snps_dataframe = pd.read_csv(csv_file, sep = "\t")

        # populate self.PositionedSNPIndexKey_to_PositionedSNP
        for i, row in snps_dataframe.iterrows():
            # get the data
            position_1 = Position(genome=genome_1, chrom=row["ref_chrom"], pos=row["ref_pos"])
            position_2 = Position(genome=genome_2, chrom=row["query_chrom"], pos=row["query_pos"])
            allele_1 = min(row["ref_sub"], row["query_sub"])
            allele_2 = max(row["ref_sub"], row["query_sub"])
            self.add_PositionedSNP(position_1, position_2, allele_1, allele_2)



    def _get_all_unique_Positioned_SNPs(self):
        # as there are no Positioned_SNPs deep copies, their memory address is the unique identifier
        all_unique_positioned_SNPs = []
        all_unique_positioned_SNPs_ids = set()
        for positioned_snp in self.PositionedSNPIndexKey_to_PositionedSNP.values():
            if id(positioned_snp) not in all_unique_positioned_SNPs_ids:
                all_unique_positioned_SNPs.append(positioned_snp)
        return all_unique_positioned_SNPs

    def get_nb_SNPs_in_pangenome(self):
        return len(self._get_all_unique_Positioned_SNPs())

    def get_nb_SNPs_that_can_be_found_with_a_given_genome(self, genome):
        all_positioned_SNPs = self._get_all_unique_Positioned_SNPs()
        nb = 0
        for positioned_SNP in all_positioned_SNPs:
            for position in positioned_SNP.positions:
                if position.genome == genome:
                    nb += 1
                    break
        return nb

    @classmethod
    def load(cls, filename):
        with open(filename, "rb") as fin:
            positionedSNPs = pickle.load(fin)
        return positionedSNPs

    def save(self, filename):
        with open(filename, "wb") as fout:
            pickle.dump(self, fout)

    @classmethod
    def from_dict(cls, dict):
        '''
        Factory method mainly used for testing
        :param dict: represents a dictionary to configure a new PositionedSNPsIndex object
        :return: new PositionedSNPsIndex object
        '''
        new_PositionedSNPsIndex = PositionedSNPsIndex()
        new_PositionedSNPsIndex.__dict__ = dict
        return new_PositionedSNPsIndex


    def __eq__(self, other):
        """Overrides the default implementation"""
        if isinstance(other, PositionedSNPsIndex):
            return self.__dict__ == other.__dict__
        return NotImplemented


    def __repr__(self):
        return str(self.__dict__)



# TODO: add unit test
# positionedSNPsIndex = PositionedSNPsIndex()
# positionedSNPsIndex.add_SNPs_from_csv('assemblies_sample_out/genome.1/genome.1-SEP-genome.2.mummer.csv', 'genome.1', 'genome.2')
# positionedSNPsIndex.add_SNPs_from_csv('assemblies_sample_out/genome.1/genome.1-SEP-genome.3.mummer.csv', 'genome.1', 'genome.3')
# positionedSNPsIndex.add_SNPs_from_csv('assemblies_sample_out/genome.2/genome.2-SEP-genome.3.mummer.csv', 'genome.2', 'genome.3')
# positionedSNPsIndex.add_SNPs_from_csv('assemblies_sample_out/genome.0/genome.0-SEP-genome.1.mummer.csv', 'genome.0', 'genome.1')
# positionedSNPsIndex.add_SNPs_from_csv('assemblies_sample_out/genome.0/genome.0-SEP-genome.2.mummer.csv', 'genome.0', 'genome.2')
# positionedSNPsIndex.add_SNPs_from_csv('assemblies_sample_out/genome.0/genome.0-SEP-genome.3.mummer.csv', 'genome.0', 'genome.3')
# positionedSNPsIndex.serialize("assemblies_sample_out/positionedSNPsIndex")
# positionedSNPs_loaded = positionedSNPsIndex.load("assemblies_sample_out/positionedSNPsIndex")
