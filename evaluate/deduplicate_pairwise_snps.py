import networkx as nx
from evaluate.mummer import ShowSNPsDataframe
from pathlib import Path
import re
import pickle
from typing import Tuple, Iterable, Set, List, Generator, Dict, Optional, BinaryIO, TextIO
import functools
from collections import defaultdict
import pandas as pd
from .probe import Probe, ProbeHeader, ProbeInterval


class DeduplicatedVariationsDataframe(pd.DataFrame):
    def get_probes(self) -> Tuple[str, str]:
        ref_probes = []
        query_probes = []

        for index, row in self.iterrows():
            ref_probe, query_probe = self._get_ref_and_query_probe(row)
            ref_probes.append(str(ref_probe))
            query_probes.append(str(query_probe))

        return (
            "\n".join(probe for probe in ref_probes if probe),
            "\n".join(probe for probe in query_probes if probe),
        )


    @property
    def _constructor(self):
        return DeduplicatedVariationsDataframe


    @staticmethod
    def _get_ref_and_query_probe(row: pd.Series) -> Tuple[Probe, ...]:
        probes = []
        probe_prefixes = ["ref", "query"]

        # TODO: is this really " - 1" or len(ref/query_sub)??
        flank_width = int((len(row[f"{probe_prefixes[0]}_context"]) - 1) / 2)

        for prefix in probe_prefixes:
            core_sequence = row[f"{prefix}_sub"].replace(".", "")
            left_flank = row[f"{prefix}_context"][:flank_width].replace("-", "")
            right_flank = row[f"{prefix}_context"][flank_width + 1 :].replace("-", "")
            call_start_idx = len(left_flank)
            call_end_idx = call_start_idx + len(core_sequence)
            header = ProbeHeader(
                sample=row[f"{prefix}_genome"],
                chrom=row[f"{prefix}_chrom"],
                pos=row[f"{prefix}_pos"],
                ref_length=len(core_sequence),
                interval=ProbeInterval(call_start_idx, call_end_idx),
                pangenome_variation_id=row["pangenome_variation_id"],
                number_of_alleles=row["number_of_alleles"],
                allele_id=row[f"{prefix}_allele_id"],
                number_of_different_allele_sequences=row["number_of_different_allele_sequences"],
                allele_sequence_id=row[f"{prefix}_allele_sequence_id"],
            )
            full_sequence = left_flank + core_sequence + right_flank
            probes.append(Probe(header=header, full_sequence=full_sequence))

        return tuple(probes)



class DeduplicatePairwiseSNPsUtils:
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
            snps_df = pickle.load(df_fh)
        snps_df = snps_df.translate_to_FWD_strand()
        return snps_df



class MPHF:
    def __init__(self):
        self._object_to_id = {}
        self._id_to_object = []

    @property
    def object_to_id(self) -> Dict[object, int]:
        return self._object_to_id
    @property
    def id_to_object(self) -> List[object]:
        return self._id_to_object

    def add_object(self, object):
        if object not in self.object_to_id:
            new_id = self.get_number_of_objects()
            self.object_to_id[object] = new_id
            self.id_to_object.append(object)

    def get_id(self, object) -> int:
        return self.object_to_id[object]

    def get_object(self, object_id: int) -> object:
        return self.id_to_object[object_id]

    def get_number_of_objects(self) -> int:
        both_DS_have_the_same_length = len(self.object_to_id) == len(self.id_to_object)
        assert both_DS_have_the_same_length
        return len(self.object_to_id)


    # serialization
    def dump(self, file_with_nb_of_objects: TextIO, pickle_file: BinaryIO):
        file_with_nb_of_objects.write(str(self.get_number_of_objects()))
        pickle.dump(self, pickle_file)
    @staticmethod
    def load(file: BinaryIO) -> "AlleleMPHF":
        return pickle.load(file)



class NotASNP(Exception):
    pass

@functools.total_ordering
class Allele:
    def __init__(self, genome: str, chrom: str, pos: int, sequence: str):
        self._genome = genome
        self._chrom = chrom
        self._pos = pos  # TODO: do we need to care about strand when looking at pos?
        self._sequence = sequence

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
    def sequence(self) -> str:
        return self._sequence
    @property
    def data_tuple(self) -> Tuple[str,str,int, str]:
        """
        :return: a tuple with all the data in the allele
        """
        return self.genome, self.chrom, self.pos, self.sequence

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

    def __repr__(self):
        return str(vars(self))

    # Note: tested through DeduplicationGraph._add_variants_from_ShowSNPsDataframe_core()
    @staticmethod
    def get_alleles_from_ShowSNPsDataframe(ref: str, query: str, snps_df: ShowSNPsDataframe) -> Generator[Tuple["Allele", "Allele"], None, None]:
        for ref_chrom, ref_pos, ref_sub, query_chrom, query_pos, query_sub in \
            zip(snps_df["ref_chrom"], snps_df["ref_pos"], snps_df["ref_sub"],
                snps_df["query_chrom"], snps_df["query_pos"], snps_df["query_sub"]):
            is_snp = len(ref_sub)==1 and len(query_sub)==1
            if not is_snp:
                continue  # we just deal with SNPs as of now

            ref_allele = Allele(ref, ref_chrom, ref_pos, ref_sub)
            query_allele = Allele(query, query_chrom, query_pos, query_sub)
            yield ref_allele, query_allele



class AlleleMPHF(MPHF):
    def __init__(self):
        super().__init__()

    # helpers
    def _add_variants_from_ShowSNPsDataframe_core(self, ref: str, query: str, snps_df: ShowSNPsDataframe):
        for ref_allele, query_allele in Allele.get_alleles_from_ShowSNPsDataframe(ref, query, snps_df):
            self.add_object(ref_allele)
            self.add_object(query_allele)

    def _add_variants_from_ShowSNPsDataframe_filepath(self, ShowSNPsDataframe_filepath):
        ref, query = DeduplicatePairwiseSNPsUtils._get_ref_and_query_from_ShowSNPsDataframe_filepath(
            ShowSNPsDataframe_filepath)
        snps_df = DeduplicatePairwiseSNPsUtils._load_pickled_ShowSNPsDataframe(ShowSNPsDataframe_filepath)
        self._add_variants_from_ShowSNPsDataframe_core(ref, query, snps_df)

    @staticmethod
    def build_from_list_of_snps_dfs_filepaths(snps_dfs_filepaths: List[str]) -> "AlleleMPHF":
        allele_mphf = AlleleMPHF()
        for snps_df_filepath in snps_dfs_filepaths:
            allele_mphf._add_variants_from_ShowSNPsDataframe_filepath(snps_df_filepath)
        return allele_mphf



@functools.total_ordering
class PairwiseVariation:
    """
    Pairwise variation does not know alleles, only allele IDs
    """
    def __init__(self, ref_allele_id: int, query_allele_id: int,
                 allele_mphf: AlleleMPHF):
        # this ordering is done to facilitate the usage of this class, for the equal and hash functions
        self._allele_1_id = min(ref_allele_id, query_allele_id)
        self._allele_2_id = max(ref_allele_id, query_allele_id)
        self._original_ref_allele_id = ref_allele_id
        self._original_query_allele_id = query_allele_id
        self._allele_mphf = allele_mphf

    @property
    def allele_1_id(self) -> int:
        return self._allele_1_id
    @property
    def allele_2_id(self) -> int:
        return self._allele_2_id
    @property
    def original_ref_allele_id(self) -> int:
        return self._original_ref_allele_id
    @property
    def original_query_allele_id(self) -> int:
        return self._original_query_allele_id
    @property
    def allele_1(self) -> Allele:
        return self._allele_mphf.get_object(self.allele_1_id)
    @property
    def allele_2(self) -> Allele:
        return self._allele_mphf.get_object(self.allele_2_id)
    @property
    def original_ref_allele(self) -> Allele:
        return self._allele_mphf.get_object(self.original_ref_allele_id)
    @property
    def original_query_allele(self) -> Allele:
        return self._allele_mphf.get_object(self.original_query_allele_id)

    def __eq__(self, other: object) -> bool:
        if isinstance(other, PairwiseVariation):
            return (self.allele_1_id, self.allele_2_id) == (other.allele_1_id, other.allele_2_id)
        else:
            return False

    def __hash__(self) -> int:
        return hash((self.allele_1_id, self.allele_2_id))

    def __lt__(self, other: object) -> bool:
        if isinstance(other, PairwiseVariation):
            return (self.allele_1_id, self.allele_2_id) < (other.allele_1_id, other.allele_2_id)
        else:
            raise TypeError()

    def __repr__(self):
        return str(vars(self))

    def share_allele(self, other: "PairwiseVariation") -> bool:
        return self.allele_1_id == other.allele_1_id or self.allele_1_id == other.allele_2_id or \
               self.allele_2_id == other.allele_1_id or self.allele_2_id == other.allele_2_id

    # Note: tested through DeduplicationGraph._add_variants_from_ShowSNPsDataframe_core()
    @staticmethod
    def get_PairwiseVariation_from_ShowSNPsDataframe(ref: str, query: str, snps_df: ShowSNPsDataframe,
                                                     allele_mphf: AlleleMPHF) -> Generator["PairwiseVariation", None, None]:
        for ref_allele, query_allele in Allele.get_alleles_from_ShowSNPsDataframe(ref, query, snps_df):
            ref_allele_id = allele_mphf.get_id(ref_allele)
            query_allele_id = allele_mphf.get_id(query_allele)
            yield PairwiseVariation(ref_allele_id, query_allele_id, allele_mphf)



class PairwiseVariationMPHF(MPHF):
    def __init__(self):
        super().__init__()

    # helpers
    def _add_variants_from_ShowSNPsDataframe_core(self, ref: str, query: str, snps_df: ShowSNPsDataframe, allele_mphf: AlleleMPHF):
        for pairwise_variation in PairwiseVariation.get_PairwiseVariation_from_ShowSNPsDataframe(ref, query, snps_df,
                                                                                                 allele_mphf):
            self.add_object(pairwise_variation)

    def _add_variants_from_ShowSNPsDataframe_filepath(self, ShowSNPsDataframe_filepath: str, allele_mphf: AlleleMPHF):
        ref, query = DeduplicatePairwiseSNPsUtils._get_ref_and_query_from_ShowSNPsDataframe_filepath(
            ShowSNPsDataframe_filepath)
        snps_df = DeduplicatePairwiseSNPsUtils._load_pickled_ShowSNPsDataframe(ShowSNPsDataframe_filepath)
        self._add_variants_from_ShowSNPsDataframe_core(ref, query, snps_df, allele_mphf)

    @staticmethod
    def build_from_list_of_snps_dfs_filepaths(snps_dfs_filepaths: List[str], allele_mphf_filepath: str) -> "PairwiseVariationMPHF":
        with open(allele_mphf_filepath, "rb") as allele_mphf_file:
            allele_mphf = AlleleMPHF.load(allele_mphf_file)

        pairwise_variation_mphf = PairwiseVariationMPHF()
        for snps_df_filepath in snps_dfs_filepaths:
            pairwise_variation_mphf._add_variants_from_ShowSNPsDataframe_filepath(snps_df_filepath, allele_mphf)
        return pairwise_variation_mphf

    def get_pairwise_variation_id_to_alleles_id(self) -> List[Tuple[int, int]]:
        pairwise_variation_id_to_alleles_id = []
        for pairwise_variation in self.id_to_object:
            pairwise_variation_id_to_alleles_id.append((pairwise_variation.allele_1_id, pairwise_variation.allele_2_id))
        return pairwise_variation_id_to_alleles_id

    def dump(self, file_with_nb_of_objects: TextIO, pickle_file: BinaryIO, pairwise_variation_id_to_alleles_id_file: BinaryIO):
        super().dump(file_with_nb_of_objects, pickle_file)
        pairwise_variation_id_to_alleles_id = self.get_pairwise_variation_id_to_alleles_id()
        pickle.dump(pairwise_variation_id_to_alleles_id, pairwise_variation_id_to_alleles_id_file)



class DeduplicationGraph:
    def __init__(self, number_of_alleles: int, pairwise_variation_id_to_alleles_id: List[Tuple[int, int]]):
        self._number_of_alleles = number_of_alleles
        self._pairwise_variation_id_to_alleles_id = pairwise_variation_id_to_alleles_id
        self._build_graph_nodes()
        self._index()
        self._build_edges()
    @property
    def number_of_alleles(self) -> int:
        return self._number_of_alleles
    @property
    def pairwise_variation_id_to_alleles_id(self) -> List[Tuple[int, int]]:
        return self._pairwise_variation_id_to_alleles_id
    @property
    def number_of_pairwise_variations(self) -> int:
        return len(self.pairwise_variation_id_to_alleles_id)
    @property
    def graph(self) -> nx.Graph:
        return self._graph
    @property
    def nodes(self):
        return self.graph.nodes
    @property
    def edges(self):
        return self.graph.edges
    @property
    def allele_to_pairwise_variations(self) -> List[Set[int]]:
        return self._allele_to_pairwise_variations


    def _build_graph_nodes(self):
        self._graph = nx.Graph()
        self._graph.add_nodes_from(range(self.number_of_pairwise_variations))

    def _index(self):
        self._allele_to_pairwise_variations = [set() for _ in range(self.number_of_alleles)]
        for pairwise_variation_id, (allele_id_1, allele_id_2) in enumerate(self.pairwise_variation_id_to_alleles_id):
            self._allele_to_pairwise_variations[allele_id_1].add(pairwise_variation_id)
            self._allele_to_pairwise_variations[allele_id_2].add(pairwise_variation_id)

    def _add_edge(self, variant_1, variant_2) -> None:
        self._graph.add_edge(variant_1, variant_2)

    def _build_edges(self) -> None:
        """
        Note: Should be called after all nodes are added.
        """
        for pairwise_variations in self.allele_to_pairwise_variations:
            if len(pairwise_variations) > 1:
                # connect the variations with a path
                pairwise_variations_as_list = list(pairwise_variations)
                for pairwise_variation_1, pairwise_variation_2 in \
                    zip(pairwise_variations_as_list, pairwise_variations_as_list[1:]):
                    self._add_edge(pairwise_variation_1, pairwise_variation_2)


    # Note: not tested (trivial method)
    def _get_connected_components(self) -> Generator[Set[int], None, None]:
        return nx.connected_components(self.graph)

    def get_pangenome_variations_defined_by_allele_ids(self) -> List[Set[int]]:
        pangenome_variations_defined_by_allele_ids = []

        connected_components = self._get_connected_components()
        for connected_component_index, connected_component in enumerate(connected_components):
            allele_ids_in_connected_component = []
            for pairwise_variation_id in connected_component:
                allele_ids_in_connected_component.extend(self.pairwise_variation_id_to_alleles_id[pairwise_variation_id])
            pangenome_variation = set(allele_ids_in_connected_component)
            pangenome_variations_defined_by_allele_ids.append(pangenome_variation)

        return pangenome_variations_defined_by_allele_ids

    def __repr__(self):
        return str(vars(self))



class PangenomeVariation:
    # Note: trivial method, not tested
    def __init__(self, id: int, alleles: Iterable[Allele]):
        self._id = id
        self._alleles = sorted(list(set(alleles)))
        self._unique_allele_sequences = self._get_unique_allele_sequences()

    @property
    def id(self) -> int:
        return self._id
    @property
    def alleles(self) -> List[Allele]:
        return self._alleles
    @property
    def unique_allele_sequences(self) -> List[str]:
        return self._unique_allele_sequences

    # Note: trivial method, not tested
    def __eq__(self, other: object):
        if isinstance(other, PangenomeVariation):
            return self.id==other.id and self.alleles==other.alleles
        else:
            return False

    def _get_unique_allele_sequences(self) -> List[str]:
        """
        Get the set of different allele sequences in this Pangenome Variation.
        If this is a SNP A -> C, then we have 2 different allele sequences (A and C)
        If this Pangenome Variations has all possible SNPs, then we would have 4 different allele sequences (ACGT)
        """
        set_of_unique_allele_sequences = {allele.sequence for allele in self.alleles}
        unique_allele_sequences = sorted(list(set_of_unique_allele_sequences))
        return unique_allele_sequences

    def is_consistent(self) -> bool:
        genomes_to_chrom_and_pos = defaultdict(set)
        for allele in self.alleles:
            genomes_to_chrom_and_pos[allele.genome].add((allele.chrom, allele.pos))
        genomes_to_have_consistent_alleles = {
            genome: len(chrom_and_pos)==1 for genome, chrom_and_pos in genomes_to_chrom_and_pos.items()
        }
        all_genomes_have_consistent_alleles = all(genomes_to_have_consistent_alleles.values())
        return all_genomes_have_consistent_alleles

    # Note: trivial getters, not tested:
    def get_number_of_alleles(self) -> int:
        return len(self.alleles)
    def get_allele_index(self, allele: Allele) -> int:
        return self.alleles.index(allele)
    def get_number_of_different_allele_sequences(self) -> int:
        return len(self.unique_allele_sequences)
    def get_allele_sequence_index(self, allele: Allele) -> int:
        return self.unique_allele_sequences.index(allele.sequence)

    def __repr__(self):
        return str(vars(self))

# Note: trivial class, not tested
class PangenomeVariations:
    """
    Stores a list of PangenomeVariation
    """
    def __init__(self):
        self._pangenome_variations = []
    @property
    def pangenome_variations(self) -> List[PangenomeVariation]:
        return self._pangenome_variations
    def append(self, pangenome_variation: PangenomeVariation):
        self.pangenome_variations.append(pangenome_variation)
    def __eq__(self, other: object):
        if isinstance(other, PangenomeVariations):
            return self.pangenome_variations == other.pangenome_variations
        else:
            return False

    def __repr__(self):
        return str(vars(self))

    @staticmethod
    def build_from_pangenome_variations_defined_by_allele_ids(pangenome_variations_defined_by_allele_ids: List[Set[int]],
                                                              allele_mphf: AlleleMPHF) -> "PangenomeVariations":
        pangenome_variations = PangenomeVariations()
        for pangenome_variation_index, alleles_indexes in enumerate(pangenome_variations_defined_by_allele_ids):
            alleles = [allele_mphf.get_object(allele_index) for allele_index in alleles_indexes]
            pangenome_variation = PangenomeVariation(pangenome_variation_index, alleles)
            pangenome_variations.append(pangenome_variation)
        return pangenome_variations

class InconsistentPangenomeVariations(Exception):
    pass

class ConsistentPangenomeVariations:
    """
    Represents a list of ConsistentPangenomeVariations (built from PangenomeVariations)
    For the definition of consistent Pangenome Variations, see https://github.com/iqbal-lab/pandora1_paper/issues/144#issue-603283664
    """
    def __init__(self, pangenome_variations: PangenomeVariations):
        self._consistent_pangenome_variations = [
            pangenome_variation \
            for pangenome_variation in pangenome_variations.pangenome_variations \
            if pangenome_variation.is_consistent()
        ]

        # this allele to consistent_pangenome_variations indexing is to speedup some methods
        self._alleles_to_consistent_pangenome_variations = defaultdict(lambda: None)
        for consistent_pangenome_variation in self.consistent_pangenome_variations:
            for allele in consistent_pangenome_variation.alleles:
                self._alleles_to_consistent_pangenome_variations[allele] = consistent_pangenome_variation

    @property
    def consistent_pangenome_variations(self) -> List[PangenomeVariation]:
        return self._consistent_pangenome_variations
    @property
    def alleles_to_consistent_pangenome_variations(self) -> Dict[Allele, Optional[PangenomeVariation]]:
        return self._alleles_to_consistent_pangenome_variations


    def get_consistent_pangenome_variation(self, pairwise_variation: PairwiseVariation) -> Optional[PangenomeVariation]:
        # Note: pangenome_variation_of_allele_1/2 can be None
        pangenome_variation_of_allele_1 = self.alleles_to_consistent_pangenome_variations[pairwise_variation.allele_1]
        pangenome_variation_of_allele_2 = self.alleles_to_consistent_pangenome_variations[pairwise_variation.allele_2]

        both_alleles_have_the_same_pangenome_variation = pangenome_variation_of_allele_1 == pangenome_variation_of_allele_2
        if not both_alleles_have_the_same_pangenome_variation:
            raise InconsistentPangenomeVariations()

        return pangenome_variation_of_allele_1


    def _get_DeduplicatedVariationsDataframe(self, ref: str, query: str, snps_df: ShowSNPsDataframe,
                                             allele_mphf: AlleleMPHF) -> DeduplicatedVariationsDataframe:
        """
        Builds a DeduplicatedVariationsDataframe from a ShowSNPsDataframe with info computed from the ConsistentPangenomeVariations.
        ** WARNING: this also modifies snps_df parameter, we dont want to make a copy**

        Adds the following columns to snps_df:
        ref_genome: str
        query_genome: str
        present_in_a_consistent_pangenome_variation: bool
        pangenome_variation_id: int
        number_of_alleles: int
        ref_allele_id: int
        query_allele_id: int
        number_of_different_allele_sequences: int
        ref_allele_sequence_id: int
        query_allele_sequence_id: int

        :return the DeduplicatedVariationsDataframe
        """
        ref_genome = [ref]*len(snps_df)
        query_genome = [query]*len(snps_df)
        present_in_a_consistent_pangenome_variation = [False]*len(snps_df)
        pangenome_variation_id = [-1]*len(snps_df)
        number_of_alleles = [-1]*len(snps_df)
        ref_allele_id = [-1]*len(snps_df)
        query_allele_id = [-1] * len(snps_df)
        number_of_different_allele_sequences = [-1]*len(snps_df)
        ref_allele_sequence_id = [-1]*len(snps_df)
        query_allele_sequence_id = [-1] * len(snps_df)

        for index, pairwise_variation in enumerate(PairwiseVariation.get_PairwiseVariation_from_ShowSNPsDataframe(ref, query, snps_df, allele_mphf)):
            consistent_pangenome_variation = self.get_consistent_pangenome_variation(pairwise_variation)
            is_present = consistent_pangenome_variation is not None
            if is_present:
                present_in_a_consistent_pangenome_variation[index] = True
                pangenome_variation_id[index] = consistent_pangenome_variation.id
                number_of_alleles[index] = consistent_pangenome_variation.get_number_of_alleles()
                ref_allele_id[index] = consistent_pangenome_variation.get_allele_index(pairwise_variation.original_ref_allele)
                query_allele_id[index] = consistent_pangenome_variation.get_allele_index(pairwise_variation.original_query_allele)
                number_of_different_allele_sequences[index] = consistent_pangenome_variation.get_number_of_different_allele_sequences()
                ref_allele_sequence_id[index] = consistent_pangenome_variation.get_allele_sequence_index(pairwise_variation.original_ref_allele)
                query_allele_sequence_id[index] = consistent_pangenome_variation.get_allele_sequence_index(pairwise_variation.original_query_allele)

        deduplicated_snps_df = DeduplicatedVariationsDataframe(snps_df)
        deduplicated_snps_df["ref_genome"] = ref_genome
        deduplicated_snps_df["query_genome"] = query_genome
        deduplicated_snps_df["present_in_a_consistent_pangenome_variation"] = present_in_a_consistent_pangenome_variation
        deduplicated_snps_df["pangenome_variation_id"] = pangenome_variation_id
        deduplicated_snps_df["number_of_alleles"] = number_of_alleles
        deduplicated_snps_df["ref_allele_id"] = ref_allele_id
        deduplicated_snps_df["query_allele_id"] = query_allele_id
        deduplicated_snps_df["number_of_different_allele_sequences"] = number_of_different_allele_sequences
        deduplicated_snps_df["ref_allele_sequence_id"] = ref_allele_sequence_id
        deduplicated_snps_df["query_allele_sequence_id"] = query_allele_sequence_id
        return deduplicated_snps_df


    def build_DeduplicatedVariationsDataframe_from_ShowSNPsDataframe(self, ShowSNPsDataframe_filepath: str,
                                                                     allele_mphf: AlleleMPHF) -> DeduplicatedVariationsDataframe:
        """
        Loads a ShowSNPsDataframe, add all the relevant information about Consistent Pangenome Variations into it,
        builds the DeduplicatedVariationsDataframe, and filter out variations that are not in a Consistent Pangenome Variation
        """
        ref, query = DeduplicatePairwiseSNPsUtils._get_ref_and_query_from_ShowSNPsDataframe_filepath(ShowSNPsDataframe_filepath)
        snps_df = DeduplicatePairwiseSNPsUtils._load_pickled_ShowSNPsDataframe(ShowSNPsDataframe_filepath)
        deduplicated_snps_df = self._get_DeduplicatedVariationsDataframe(ref, query, snps_df, allele_mphf)
        filtered_snps_df = deduplicated_snps_df[deduplicated_snps_df.present_in_a_consistent_pangenome_variation == True]
        filtered_snps_df.reset_index(drop=True, inplace=True)
        return DeduplicatedVariationsDataframe(filtered_snps_df)

    def __repr__(self):
        return str(vars(self))
