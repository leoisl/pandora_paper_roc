import networkx as nx
from evaluate.mummer import ShowSNPsDataframe
from pathlib import Path
import re
import pickle
from typing import Tuple
import functools

class NotASNP(Exception):
    pass

@functools.total_ordering
class Allele:
    def __init__(self, genome: str, chrom: str, pos: int, sequence: str):
        self._genome = genome
        self._chrom = chrom
        self._pos = pos

        # sequence is just used to ensure this is a SNP, for now
        is_snp = len(sequence) == 1
        if not is_snp:
            raise NotASNP()

    @property
    def genome(self) -> str:
        return self._genome
    @property
    def chrom(self) -> str:
        return self._chrom
    @property
    def pos(self) -> int:
        return self._pos
    @property
    def data_tuple(self) -> Tuple[str,str,int]:
        """
        :return: a tuple with all the data in the allele
        """
        return self.genome, self.chrom, self.pos

    def __eq__(self, other: object) -> bool:
        if isinstance(other, Allele):
            return  self.data_tuple==other.data_tuple
        else:
            return False

    def __hash__(self) -> int:
        return hash(self.data_tuple)

    def __lt__(self, other: object) -> bool:
        if isinstance(other, Allele):
            return self.data_tuple < other.data_tuple
        else:
            raise TypeError()


class PairwiseVariation:
    def __init__(self, allele_1: Allele, allele_2: Allele):
        self._allele_1 = min(allele_1, allele_2)
        self._allele_2 = max(allele_1, allele_2)
    @property
    def allele_1(self) -> Allele:
        return self._allele_1
    @property
    def allele_2(self) -> Allele:
        return self._allele_2

    def __eq__(self, other: object) -> bool:
        if isinstance(other, PairwiseVariation):
            return self.allele_1 == other.allele_1 and self.allele_2 == other.allele_2
        else:
            return False

    def __hash__(self) -> int:
        return hash((self.allele_1, self.allele_2))

    def share_allele(self, other: "PairwiseVariation") -> bool:
        return self.allele_1 == other.allele_1 or self.allele_1 == other.allele_2 or \
               self.allele_2 == other.allele_1 or self.allele_2 == other.allele_2


class DeduplicationGraph:
    def __init__(self):
        self._graph = nx.Graph()
    @property
    def graph(self) -> nx.Graph:
        return self._graph

    @staticmethod
    def _get_ref_and_query_from_ShowSNPsDataframe_filepath(ShowSNPsDataframe_filepath: str) -> Tuple[str, str]:
        ShowSNPsDataframe_filepath = Path(ShowSNPsDataframe_filepath)
        ShowSNPsDataframe_filename = ShowSNPsDataframe_filepath.name
        matches = re.match(r"(.*)_and_(.*).snps_df.pickle", ShowSNPsDataframe_filename)
        ref = matches.group(1)
        query = matches.group(2)
        return ref, query

    # Note: not tested (trivial method)
    @staticmethod
    def _load_pickled_ShowSNPsDataframe(df_filepath: str) -> ShowSNPsDataframe:
        with open(df_filepath, "rb") as df_fh:
            return pickle.load(df_fh)

    # Note: not tested (trivial method)
    def _add_pairwise_variation(self, pairwise_variation: PairwiseVariation) -> None:
        self.graph.add_node(pairwise_variation)

    def _add_variants_from_ShowSNPsDataframe_core(self, ref: str, query: str, snps_df: ShowSNPsDataframe) -> None:
        for ref_chrom, ref_pos, ref_sub, query_chrom, query_pos, query_sub in \
            zip(snps_df["ref_chrom"], snps_df["ref_pos"], snps_df["ref_sub"],
                snps_df["query_chrom"], snps_df["query_pos"], snps_df["query_sub"]):
            is_snp = len(ref_sub)==1 and len(query_sub)==1
            if not is_snp:
                continue  # we just deal with SNPs as of now

            ref_allele = Allele(ref, ref_chrom, ref_pos, ref_sub)
            query_allele = Allele(query, query_chrom, query_pos, query_sub)
            pairwise_variation = PairwiseVariation(ref_allele, query_allele)
            self._add_pairwise_variation(pairwise_variation)


    # Note: not tested (trivial method)
    def add_variants_from_ShowSNPsDataframe(self, ShowSNPsDataframe_filepath: str) -> None:
        ref, query = DeduplicationGraph._get_ref_and_query_from_ShowSNPsDataframe_filepath(ShowSNPsDataframe_filepath)
        snps_df = DeduplicationGraph._load_pickled_ShowSNPsDataframe(ShowSNPsDataframe_filepath)
        self._add_variants_from_ShowSNPsDataframe_core(ref, query, snps_df)
