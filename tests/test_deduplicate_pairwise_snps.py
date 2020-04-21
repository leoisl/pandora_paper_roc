from evaluate.deduplicate_pairwise_snps import NotASNP, Allele, PairwiseVariation, DeduplicationGraph
import pandas as pd
import pytest
from io import StringIO

class TestAllele:
    def test_constructor_isSNP_constructorOK(self):
        Allele("genome_1", "chrom_1", 10, "A")
        assert True

    def test_constructor_isNotASNP_constructorRaisesNotASNP(self):
        with pytest.raises(NotASNP):
            Allele("genome_1", "chrom_1", 10, "AC")

    def test_constructor_isDeletion_constructorRaisesNotASNP(self):
        with pytest.raises(NotASNP):
            Allele("genome_1", "chrom_1", 10, "")

    def test_hash_equalAlleles_returnsEqualHashValues(self):
        allele_1 = Allele("genome_1", "chrom_1", 10, "A")
        allele_2 = Allele("genome_1", "chrom_1", 10, "A")
        assert hash(allele_1) == hash(allele_2)

    def test_hash_differentGenomes_returnsDifferentHashValues(self):
        allele_1 = Allele("genome_1", "chrom_1", 10, "A")
        allele_2 = Allele("genome_2", "chrom_1", 10, "A")
        assert hash(allele_1) != hash(allele_2)

    def test_hash_differentChroms_returnsDifferentHashValues(self):
        allele_1 = Allele("genome_1", "chrom_1", 10, "A")
        allele_2 = Allele("genome_1", "chrom_2", 10, "A")
        assert hash(allele_1) != hash(allele_2)

    def test_hash_differentPos_returnsDifferentHashValues(self):
        allele_1 = Allele("genome_1", "chrom_1", 10, "A")
        allele_2 = Allele("genome_1", "chrom_2", 11, "A")
        assert hash(allele_1) != hash(allele_2)

    def test___relational_operators___equalAlleles(self):
        allele_1 = Allele("genome_1", "chrom_1", 10, "A")
        allele_2 = Allele("genome_1", "chrom_1", 10, "A")
        assert (allele_1 == allele_2) == True
        assert (allele_1 != allele_2) == False
        assert (allele_1 <  allele_2) == False
        assert (allele_1 <= allele_2) == True
        assert (allele_1 >  allele_2) == False
        assert (allele_1 >= allele_2) == True

    def test___relational_operators___differentGenomes___first_is_smaller(self):
        allele_1 = Allele("genome_1", "chrom_1", 10, "A")
        allele_2 = Allele("genome_2", "chrom_1", 10, "A")
        assert (allele_1 == allele_2) == False
        assert (allele_1 != allele_2) == True
        assert (allele_1 < allele_2) == True
        assert (allele_1 <= allele_2) == True
        assert (allele_1 > allele_2) == False
        assert (allele_1 >= allele_2) == False

    def test___relational_operators___differentGenomes___first_is_larger(self):
        allele_1 = Allele("genome_2", "chrom_1", 10, "A")
        allele_2 = Allele("genome_1", "chrom_1", 10, "A")
        assert (allele_1 == allele_2) == False
        assert (allele_1 != allele_2) == True
        assert (allele_1 < allele_2) == False
        assert (allele_1 <= allele_2) == False
        assert (allele_1 > allele_2) == True
        assert (allele_1 >= allele_2) == True

    def test___relational_operators___differentChroms___first_is_smaller(self):
        allele_1 = Allele("genome_1", "chrom_1", 10, "A")
        allele_2 = Allele("genome_1", "chrom_2", 10, "A")
        assert (allele_1 == allele_2) == False
        assert (allele_1 != allele_2) == True
        assert (allele_1 < allele_2) == True
        assert (allele_1 <= allele_2) == True
        assert (allele_1 > allele_2) == False
        assert (allele_1 >= allele_2) == False

    def test___relational_operators___differentChroms___first_is_larger(self):
        allele_1 = Allele("genome_1", "chrom_2", 10, "A")
        allele_2 = Allele("genome_1", "chrom_1", 10, "A")
        assert (allele_1 == allele_2) == False
        assert (allele_1 != allele_2) == True
        assert (allele_1 < allele_2) == False
        assert (allele_1 <= allele_2) == False
        assert (allele_1 > allele_2) == True
        assert (allele_1 >= allele_2) == True

    def test___relational_operators___differentPos___first_is_smaller(self):
        allele_1 = Allele("genome_1", "chrom_1", 10, "A")
        allele_2 = Allele("genome_1", "chrom_1", 11, "A")
        assert (allele_1 == allele_2) == False
        assert (allele_1 != allele_2) == True
        assert (allele_1 < allele_2) == True
        assert (allele_1 <= allele_2) == True
        assert (allele_1 > allele_2) == False
        assert (allele_1 >= allele_2) == False

    def test___relational_operators___differentPos___first_is_larger(self):
        allele_1 = Allele("genome_1", "chrom_1", 11, "A")
        allele_2 = Allele("genome_1", "chrom_1", 10, "A")
        assert (allele_1 == allele_2) == False
        assert (allele_1 != allele_2) == True
        assert (allele_1 < allele_2) == False
        assert (allele_1 <= allele_2) == False
        assert (allele_1 > allele_2) == True
        assert (allele_1 >= allele_2) == True


class TestPairwiseVariation:
    def test___constructor___ordered_alleles(self):
        allele_1 = Allele("genome_1", "chrom_1", 10, "A")
        allele_2 = Allele("genome_2", "chrom_2", 20, "A")
        pairwise_variation = PairwiseVariation(allele_1, allele_2)
        assert pairwise_variation.allele_1 == allele_1
        assert pairwise_variation.allele_2 == allele_2

    def test___constructor___unordered_alleles(self):
        allele_1 = Allele("genome_2", "chrom_2", 20, "A")
        allele_2 = Allele("genome_1", "chrom_1", 10, "A")
        pairwise_variation = PairwiseVariation(allele_1, allele_2)
        assert pairwise_variation.allele_1 == allele_2
        assert pairwise_variation.allele_2 == allele_1

    def test___equality___equalPairwiseVariations(self):
        allele_1 = Allele("genome_2", "chrom_2", 20, "A")
        allele_2 = Allele("genome_1", "chrom_1", 10, "A")
        pairwise_variation_1 = PairwiseVariation(allele_1, allele_2)
        pairwise_variation_2 = PairwiseVariation(allele_2, allele_1)
        assert pairwise_variation_1 == pairwise_variation_2

    def test___equality___differentPairwiseVariations(self):
        allele_1 = Allele("genome_2", "chrom_2", 20, "A")
        allele_2 = Allele("genome_1", "chrom_1", 10, "A")
        allele_3 = Allele("genome_1", "chrom_1", 11, "A")
        pairwise_variation_1 = PairwiseVariation(allele_1, allele_2)
        pairwise_variation_2 = PairwiseVariation(allele_1, allele_3)
        assert pairwise_variation_1 != pairwise_variation_2

    def test_hash_equalPairwiseVariations_returnsEqualHashValues(self):
        allele_1 = Allele("genome_2", "chrom_2", 20, "A")
        allele_2 = Allele("genome_1", "chrom_1", 10, "A")
        pairwise_variation_1 = PairwiseVariation(allele_1, allele_2)
        pairwise_variation_2 = PairwiseVariation(allele_2, allele_1)
        assert hash(pairwise_variation_1) == hash(pairwise_variation_2)

    def test_hash_differentPairwiseVariations_returnsDifferentHashValues(self):
        allele_1 = Allele("genome_2", "chrom_2", 20, "A")
        allele_2 = Allele("genome_1", "chrom_1", 10, "A")
        allele_3 = Allele("genome_1", "chrom_1", 11, "A")
        pairwise_variation_1 = PairwiseVariation(allele_1, allele_2)
        pairwise_variation_2 = PairwiseVariation(allele_1, allele_3)
        assert hash(pairwise_variation_1) != hash(pairwise_variation_2)

    def test___share_allele___no_sharing(self):
        allele_1 = Allele("genome_1", "chrom_1", 1, "A")
        allele_2 = Allele("genome_1", "chrom_1", 2, "A")
        allele_3 = Allele("genome_1", "chrom_1", 3, "A")
        allele_4 = Allele("genome_1", "chrom_1", 4, "A")
        pairwise_variation_1 = PairwiseVariation(allele_1, allele_2)
        pairwise_variation_2 = PairwiseVariation(allele_3, allele_4)
        assert pairwise_variation_1.share_allele(pairwise_variation_2) == False

    def test___share_allele___one_allele_shared(self):
        allele_1 = Allele("genome_1", "chrom_1", 1, "A")
        allele_2 = Allele("genome_1", "chrom_1", 2, "A")
        allele_3 = Allele("genome_1", "chrom_1", 3, "A")
        pairwise_variation_1 = PairwiseVariation(allele_1, allele_2)
        pairwise_variation_2 = PairwiseVariation(allele_2, allele_3)
        assert pairwise_variation_1.share_allele(pairwise_variation_2) == True

    def test___share_allele___both_alleles_shared(self):
        allele_1 = Allele("genome_1", "chrom_1", 1, "A")
        allele_2 = Allele("genome_1", "chrom_1", 2, "A")
        pairwise_variation_1 = PairwiseVariation(allele_1, allele_2)
        pairwise_variation_2 = PairwiseVariation(allele_2, allele_1)
        assert pairwise_variation_1.share_allele(pairwise_variation_2) == True


class TestDeduplicationGraph:
    def test___get_ref_and_query_from_ShowSNPsDataframe_filepath___absolute_path(self):
        filepath = "/home/leandro/git/pandora1_paper/analysis_output/recall/snps_dfs/CFT073/CFT073_and_H131800734.snps_df.pickle"
        ref, query = DeduplicationGraph._get_ref_and_query_from_ShowSNPsDataframe_filepath(filepath)
        assert ref == "CFT073" and query=="H131800734"

    def test___get_ref_and_query_from_ShowSNPsDataframe_filepath___relative_path(self):
        filepath = "analysis_output/recall/snps_dfs/CFT073/CFT073_and_H131800734.snps_df.pickle"
        ref, query = DeduplicationGraph._get_ref_and_query_from_ShowSNPsDataframe_filepath(filepath)
        assert ref == "CFT073" and query=="H131800734"

    def test___get_ref_and_query_from_ShowSNPsDataframe_filepath___local_path(self):
        filepath = "CFT073_and_H131800734.snps_df.pickle"
        ref, query = DeduplicationGraph._get_ref_and_query_from_ShowSNPsDataframe_filepath(filepath)
        assert ref == "CFT073" and query=="H131800734"

    def test___get_ref_and_query_from_ShowSNPsDataframe_filepath___local_path___trivial_names(self):
        filepath = "A_and_B.snps_df.pickle"
        ref, query = DeduplicationGraph._get_ref_and_query_from_ShowSNPsDataframe_filepath(filepath)
        assert ref == "A" and query=="B"

    def test___add_variants_from_ShowSNPsDataframe_core___empty_df___no_variations_added(self):
        deduplication_graph = DeduplicationGraph()
        snps_df = pd.read_csv(StringIO(
"""ref_pos,ref_sub,query_sub,query_pos,nearest_mismatch,nearest_end,ref_len,query_len,ref_context,query_context,ref_strand,query_strand,ref_chrom,query_chrom
"""
        ))
        deduplication_graph._add_variants_from_ShowSNPsDataframe_core("genome_1", "genome_2", snps_df)
        graph_is_empty = len(deduplication_graph.graph)==0
        assert graph_is_empty

    def test___add_variants_from_ShowSNPsDataframe_core___one_variation_in_df___one_variation_added(self):
        deduplication_graph = DeduplicationGraph()
        snps_df = pd.read_csv(StringIO(
"""ref_pos,ref_sub,query_sub,query_pos,nearest_mismatch,nearest_end,ref_len,query_len,ref_context,query_context,ref_strand,query_strand,ref_chrom,query_chrom
1,A,C,2,0,0,0,0,ACGT,ACGT,1,1,chrom_1,chrom_2
"""
        ))
        deduplication_graph._add_variants_from_ShowSNPsDataframe_core("genome_1", "genome_2", snps_df)

        pairwise_variation = PairwiseVariation(Allele("genome_1", "chrom_1", 1, "A"),
                                               Allele("genome_2", "chrom_2", 2, "C"))
        assert len(deduplication_graph.graph)==1 and pairwise_variation in deduplication_graph.graph


    def test___add_variants_from_ShowSNPsDataframe_core___one_SNP_variation_four_non_SNP_variation_in_df___only_SNP_variation_added(self):
        deduplication_graph = DeduplicationGraph()
        snps_df = pd.read_csv(StringIO(
"""ref_pos,ref_sub,query_sub,query_pos,nearest_mismatch,nearest_end,ref_len,query_len,ref_context,query_context,ref_strand,query_strand,ref_chrom,query_chrom
1,,C,2,0,0,0,0,ACGT,ACGT,1,1,chrom_1,chrom_2
1,A,CC,2,0,0,0,0,ACGT,ACGT,1,1,chrom_1,chrom_2
1,A,C,2,0,0,0,0,ACGT,ACGT,1,1,chrom_1,chrom_2
1,AA,CC,2,0,0,0,0,ACGT,ACGT,1,1,chrom_1,chrom_2
1,AA,C,2,0,0,0,0,ACGT,ACGT,1,1,chrom_1,chrom_2
"""
        ), na_filter=False)
        deduplication_graph._add_variants_from_ShowSNPsDataframe_core("genome_1", "genome_2", snps_df)

        pairwise_variation = PairwiseVariation(Allele("genome_1", "chrom_1", 1, "A"),
                                               Allele("genome_2", "chrom_2", 2, "C"))
        assert len(deduplication_graph.graph)==1 and pairwise_variation in deduplication_graph.graph



    def test___add_variants_from_ShowSNPsDataframe_core___three_SNP_variation_four_non_SNP_variation_in_df___only_SNPs_variations_added(self):
        deduplication_graph = DeduplicationGraph()
        snps_df = pd.read_csv(StringIO(
"""ref_pos,ref_sub,query_sub,query_pos,nearest_mismatch,nearest_end,ref_len,query_len,ref_context,query_context,ref_strand,query_strand,ref_chrom,query_chrom
10,G,T,20,0,0,0,0,ACGT,ACGT,1,1,chrom_10,chrom_20
1,,C,2,0,0,0,0,ACGT,ACGT,1,1,chrom_1,chrom_2
1,A,CC,2,0,0,0,0,ACGT,ACGT,1,1,chrom_1,chrom_2
1,A,C,2,0,0,0,0,ACGT,ACGT,1,1,chrom_1,chrom_2
1,AA,CC,2,0,0,0,0,ACGT,ACGT,1,1,chrom_1,chrom_2
1,AA,C,2,0,0,0,0,ACGT,ACGT,1,1,chrom_1,chrom_2
1,A,C,2,0,0,0,0,ACGT,ACGT,1,1,chrom_1,chrom_2
100,A,T,200,0,0,0,0,ACGT,ACGT,1,1,chrom_100,chrom_200
"""
        ), na_filter=False)
        deduplication_graph._add_variants_from_ShowSNPsDataframe_core("genome_1", "genome_2", snps_df)


        pairwise_variation_1 = PairwiseVariation(Allele("genome_1", "chrom_10", 10, "G"),
                                                 Allele("genome_2", "chrom_20", 20, "T"))
        pairwise_variation_2 = PairwiseVariation(Allele("genome_1", "chrom_1", 1, "A"),
                                                 Allele("genome_2", "chrom_2", 2, "C"))
        pairwise_variation_3 = PairwiseVariation(Allele("genome_1", "chrom_100", 100, "A"),
                                                 Allele("genome_2", "chrom_200", 200, "T"))
        assert len(deduplication_graph.graph)==3 and \
               pairwise_variation_1 in deduplication_graph.graph and \
               pairwise_variation_2 in deduplication_graph.graph and \
               pairwise_variation_3 in deduplication_graph.graph
