import networkx as nx
from evaluate.mummer import ShowSNPsDataframe
from pathlib import Path
import re
import pickle
from typing import Tuple, Iterable, Set, List, Generator, Dict, Optional
import functools
import itertools
from collections import defaultdict


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
            return pickle.load(df_fh)


class NotASNP(Exception):
    pass

@functools.total_ordering
class Allele:
    def __init__(self, genome: str, chrom: str, pos: int, sequence: str):
        self._genome = genome
        self._chrom = chrom
        self._pos = pos
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


class PairwiseVariation:
    def __init__(self, ref_allele: Allele, query_allele: Allele):
        # this ordering is done to facilitate the usage of this class, for the equal and hash functions
        self._allele_1 = min(ref_allele, query_allele)
        self._allele_2 = max(ref_allele, query_allele)

        # keep the data given to this object anyway
        self._original_ref_allele = ref_allele
        self._original_query_allele = query_allele

    @property
    def allele_1(self) -> Allele:
        return self._allele_1
    @property
    def allele_2(self) -> Allele:
        return self._allele_2
    @property
    def original_ref_allele(self) -> Allele:
        return self._original_ref_allele
    @property
    def original_query_allele(self) -> Allele:
        return self._original_query_allele

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

    # Note: tested through DeduplicationGraph._add_variants_from_ShowSNPsDataframe_core()
    @staticmethod
    def get_PairwiseVariation_from_ShowSNPsDataframe(ref: str, query: str, snps_df: ShowSNPsDataframe) -> Generator[
        "PairwiseVariation", None, None]:
        for ref_allele, query_allele in Allele.get_alleles_from_ShowSNPsDataframe(ref, query, snps_df):
            yield PairwiseVariation(ref_allele, query_allele)

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


    def _enrich_ShowSNPsDataframe(self, ref: str, query: str, snps_df: ShowSNPsDataframe) -> ShowSNPsDataframe:
        """
        Enriches a ShowSNPsDataframe with info computed from the given ConsistentPangenomeVariations.
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

        :return the enriched snps_df
        """
        enriched_snps_df = snps_df
        ref_genome = [ref]*len(enriched_snps_df)
        query_genome = [query]*len(enriched_snps_df)
        present_in_a_consistent_pangenome_variation = [False]*len(enriched_snps_df)
        pangenome_variation_id = [-1]*len(enriched_snps_df)
        number_of_alleles = [-1]*len(enriched_snps_df)
        ref_allele_id = [-1]*len(enriched_snps_df)
        query_allele_id = [-1] * len(enriched_snps_df)
        number_of_different_allele_sequences = [-1]*len(enriched_snps_df)
        ref_allele_sequence_id = [-1]*len(enriched_snps_df)
        query_allele_sequence_id = [-1] * len(enriched_snps_df)

        for index, pairwise_variation in enumerate(PairwiseVariation.get_PairwiseVariation_from_ShowSNPsDataframe(ref, query, enriched_snps_df)):
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

        enriched_snps_df["ref_genome"] = ref_genome
        enriched_snps_df["query_genome"] = query_genome
        enriched_snps_df["present_in_a_consistent_pangenome_variation"] = present_in_a_consistent_pangenome_variation
        enriched_snps_df["pangenome_variation_id"] = pangenome_variation_id
        enriched_snps_df["number_of_alleles"] = number_of_alleles
        enriched_snps_df["ref_allele_id"] = ref_allele_id
        enriched_snps_df["query_allele_id"] = query_allele_id
        enriched_snps_df["number_of_different_allele_sequences"] = number_of_different_allele_sequences
        enriched_snps_df["ref_allele_sequence_id"] = ref_allele_sequence_id
        enriched_snps_df["query_allele_sequence_id"] = query_allele_sequence_id
        return enriched_snps_df


    def load_and_process_ShowSNPsDataframe(self, ShowSNPsDataframe_filepath: str) -> ShowSNPsDataframe:
        """
        Loads a ShowSNPsDataframe, add all the relevant information about Consistent Pangenome Variations into it,
        and filter out variations that are not in a Consistent Pangenome Variation
        """
        ref, query = DeduplicatePairwiseSNPsUtils._get_ref_and_query_from_ShowSNPsDataframe_filepath(ShowSNPsDataframe_filepath)
        snps_df = DeduplicatePairwiseSNPsUtils._load_pickled_ShowSNPsDataframe(ShowSNPsDataframe_filepath)
        enriched_snps_df = self._enrich_ShowSNPsDataframe(ref, query, snps_df)
        filtered_snps_df = enriched_snps_df[enriched_snps_df.present_in_a_consistent_pangenome_variation == True]
        filtered_snps_df.reset_index(drop=True, inplace=True)
        return ShowSNPsDataframe(filtered_snps_df)

    def __repr__(self):
        return str(vars(self))

class DeduplicationGraph:
    def __init__(self):
        self._graph = nx.Graph()
    @property
    def graph(self):
        return self._graph
    @property
    def nodes(self):
        return self.graph.nodes
    @property
    def edges(self):
        return self.graph.edges

    # Note: not tested (trivial method)
    def _add_pairwise_variation(self, pairwise_variation: PairwiseVariation) -> None:
        self._graph.add_node(pairwise_variation)


    def _add_variants_from_ShowSNPsDataframe_core(self, ref: str, query: str, snps_df: ShowSNPsDataframe) -> None:
        for pairwise_variation in PairwiseVariation.get_PairwiseVariation_from_ShowSNPsDataframe(ref, query, snps_df):
            self._add_pairwise_variation(pairwise_variation)


    # Note: not tested (trivial method)
    def add_variants_from_ShowSNPsDataframe_filepath(self, ShowSNPsDataframe_filepath: str) -> None:
        ref, query = DeduplicatePairwiseSNPsUtils._get_ref_and_query_from_ShowSNPsDataframe_filepath(ShowSNPsDataframe_filepath)
        snps_df = DeduplicatePairwiseSNPsUtils._load_pickled_ShowSNPsDataframe(ShowSNPsDataframe_filepath)
        self._add_variants_from_ShowSNPsDataframe_core(ref, query, snps_df)


    def _add_edge(self, variant_1, variant_2) -> None:
        self._graph.add_edge(variant_1, variant_2)


    def build_edges(self) -> None:
        """
        Note: Should be called after all nodes are added.
        """
        #TODO: the runtime is quadratic on the number of variants. In the 24-way will be too much. We will have 24*23/2=276
        #TODO: pairwise comparisons. These comparisons have around 100k variants, so the code inside this for will be executed
        #TODO: (276 * 100000) ** 2 = ~761 trillion times... probably the bottleneck will be here.
        #TODO: implement an indexing of variants per (genome, chrom, position/10000)
        #TODO: a variant only needs to check the buckets belonging to its alleles
        for variant_1, variant_2 in itertools.product(self.nodes, self.nodes):
            if variant_1 != variant_2:
                if variant_1.share_allele(variant_2):
                    self._add_edge(variant_1, variant_2)


    # Note: not tested (trivial method)
    def _get_connected_components(self) -> Generator[Set[PairwiseVariation], None, None]:
        return nx.connected_components(self.graph)

    def get_pangenome_variations(self) -> PangenomeVariations:
        pangenome_variations = PangenomeVariations()

        connected_components = self._get_connected_components()
        for connected_component_index, connected_component in enumerate(connected_components):
            alleles_in_connected_component = []
            for pairwise_variation in connected_component:
                alleles_in_connected_component.append(pairwise_variation.allele_1)
                alleles_in_connected_component.append(pairwise_variation.allele_2)
            pangenome_variation = PangenomeVariation(connected_component_index, alleles_in_connected_component)
            pangenome_variations.append(pangenome_variation)

        return pangenome_variations

    def __repr__(self):
        return str(vars(self))

class Deduplicator:
    # Note: not tested, this is API usage, to be invoked by the rule (TODO: maybe not add this to evaluate)?
    @staticmethod
    def deduplicate(ShowSNPsDataframe_filepaths: List[str]):
        # create the deduplication graph
        deduplication_graph = DeduplicationGraph()
        for ShowSNPsDataframe_filepath in ShowSNPsDataframe_filepaths:
            deduplication_graph.add_variants_from_ShowSNPsDataframe_filepath(ShowSNPsDataframe_filepath)
        deduplication_graph.build_edges()

        # create the consistent pangenome variations
        pangenome_variations = deduplication_graph.get_pangenome_variations()
        consistent_pangenome_variations = ConsistentPangenomeVariations(pangenome_variations)

        # write the enriched and filtered variations
        for ShowSNPsDataframe_filepath in ShowSNPsDataframe_filepaths:
            filtered_snps_df = consistent_pangenome_variations.load_and_process_ShowSNPsDataframe(ShowSNPsDataframe_filepath)
            with open(f"{ShowSNPsDataframe_filepath}.deduplicated.pickle", "wb") as filtered_snps_df_fh:
                pickle.dump(filtered_snps_df, file=filtered_snps_df_fh)
            filtered_snps_df.to_csv(f"{ShowSNPsDataframe_filepath}.deduplicated.csv", index=False)
