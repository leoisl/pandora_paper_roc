from evaluate.deduplicate_pairwise_snps import DeduplicatePairwiseSNPsUtils, NotASNP, Allele, PairwiseVariation, \
    DeduplicationGraph, PangenomeVariation, PangenomeVariations, ConsistentPangenomeVariations, InconsistentPangenomeVariations
import pandas as pd
from io import StringIO
from unittest.mock import patch, PropertyMock, Mock
from unittest import TestCase
from collections import defaultdict
from evaluate.mummer import ShowSNPsDataframe


class TestDeduplicatePairwiseSNPsUtils(TestCase):
    def test___get_ref_and_query_from_ShowSNPsDataframe_filepath___absolute_path(self):
        filepath = "/home/leandro/git/pandora1_paper/analysis_output/recall/snps_dfs/CFT073/CFT073_and_H131800734.snps_df.pickle"
        ref, query = DeduplicatePairwiseSNPsUtils._get_ref_and_query_from_ShowSNPsDataframe_filepath(filepath)
        assert ref == "CFT073" and query=="H131800734"

    def test___get_ref_and_query_from_ShowSNPsDataframe_filepath___relative_path(self):
        filepath = "analysis_output/recall/snps_dfs/CFT073/CFT073_and_H131800734.snps_df.pickle"
        ref, query = DeduplicatePairwiseSNPsUtils._get_ref_and_query_from_ShowSNPsDataframe_filepath(filepath)
        assert ref == "CFT073" and query=="H131800734"

    def test___get_ref_and_query_from_ShowSNPsDataframe_filepath___local_path(self):
        filepath = "CFT073_and_H131800734.snps_df.pickle"
        ref, query = DeduplicatePairwiseSNPsUtils._get_ref_and_query_from_ShowSNPsDataframe_filepath(filepath)
        assert ref == "CFT073" and query=="H131800734"

    def test___get_ref_and_query_from_ShowSNPsDataframe_filepath___local_path___trivial_names(self):
        filepath = "A_and_B.snps_df.pickle"
        ref, query = DeduplicatePairwiseSNPsUtils._get_ref_and_query_from_ShowSNPsDataframe_filepath(filepath)
        assert ref == "A" and query=="B"



class TestAllele(TestCase):
    def test_constructor_isSNP_constructorOK(self):
        Allele("genome_1", "chrom_1", 10, "A")
        self.assertTrue(True)

    def test_constructor_isNotASNP_constructorRaisesNotASNP(self):
        with self.assertRaises(NotASNP):
            Allele("genome_1", "chrom_1", 10, "AC")

    def test_constructor_isDeletion_constructorRaisesNotASNP(self):
        with self.assertRaises(NotASNP):
            Allele("genome_1", "chrom_1", 10, "")

    def test_hash_equalAlleles_returnsEqualHashValues(self):
        allele_1 = Allele("genome_1", "chrom_1", 10, "A")
        allele_2 = Allele("genome_1", "chrom_1", 10, "A")
        self.assertEqual(hash(allele_1), hash(allele_2))

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
        allele_2 = Allele("genome_1", "chrom_1", 11, "A")
        assert hash(allele_1) != hash(allele_2)

    def test_hash_differentSeqs_returnsDifferentHashValues(self):
        allele_1 = Allele("genome_1", "chrom_1", 10, "A")
        allele_2 = Allele("genome_1", "chrom_1", 10, "C")
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

    def test___relational_operators___differentSeqs___first_is_smaller(self):
        allele_1 = Allele("genome_1", "chrom_1", 10, "A")
        allele_2 = Allele("genome_1", "chrom_1", 10, "C")
        assert (allele_1 == allele_2) == False
        assert (allele_1 != allele_2) == True
        assert (allele_1 < allele_2) == True
        assert (allele_1 <= allele_2) == True
        assert (allele_1 > allele_2) == False
        assert (allele_1 >= allele_2) == False

    def test___relational_operators___differentSeqs___first_is_larger(self):
        allele_1 = Allele("genome_1", "chrom_1", 10, "C")
        allele_2 = Allele("genome_1", "chrom_1", 10, "A")
        assert (allele_1 == allele_2) == False
        assert (allele_1 != allele_2) == True
        assert (allele_1 < allele_2) == False
        assert (allele_1 <= allele_2) == False
        assert (allele_1 > allele_2) == True
        assert (allele_1 >= allele_2) == True


class TestPairwiseVariation(TestCase):
    def test___constructor___ordered_alleles(self):
        allele_1 = Allele("genome_1", "chrom_1", 10, "A")
        allele_2 = Allele("genome_2", "chrom_2", 20, "A")
        pairwise_variation = PairwiseVariation(allele_1, allele_2)
        assert pairwise_variation.allele_1 == allele_1
        assert pairwise_variation.allele_2 == allele_2
        assert pairwise_variation.original_ref_allele == allele_1
        assert pairwise_variation.original_query_allele == allele_2

    def test___constructor___unordered_alleles(self):
        allele_1 = Allele("genome_2", "chrom_2", 20, "A")
        allele_2 = Allele("genome_1", "chrom_1", 10, "A")
        pairwise_variation = PairwiseVariation(allele_1, allele_2)
        assert pairwise_variation.allele_1 == allele_2
        assert pairwise_variation.allele_2 == allele_1
        assert pairwise_variation.original_ref_allele == allele_1
        assert pairwise_variation.original_query_allele == allele_2

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


class TestPangenomeVariation(TestCase):
    def test____get_unique_allele_sequences___one_allele(self, *mocks):
        pangenome_variation = PangenomeVariation(0, [Allele("genome_1", "chrom_1", 1, "A")])
        actual = pangenome_variation.unique_allele_sequences
        expected=["A"]
        self.assertListEqual(actual, expected)

    def test____get_unique_allele_sequences___two_unique_alleles_with_several_alleles(self, *mocks):
        pangenome_variation = PangenomeVariation(0, [Allele("genome_1", "chrom_1", 1, "G"),
                                                     Allele("genome_1", "chrom_1", 1, "G"),
                                                     Allele("genome_1", "chrom_1", 1, "G"),
                                                     Allele("genome_1", "chrom_1", 1, "A"),
                                                     Allele("genome_1", "chrom_1", 1, "G"),
                                                     Allele("genome_1", "chrom_1", 1, "A"),])
        actual = pangenome_variation.unique_allele_sequences
        expected=["A", "G"]
        self.assertListEqual(actual, expected)

    def test____get_unique_allele_sequences___three_unique_alleles_with_several_alleles(self, *mocks):
        pangenome_variation = PangenomeVariation(0, [Allele("genome_1", "chrom_1", 1, "G"),
                                                     Allele("genome_1", "chrom_1", 1, "G"),
                                                     Allele("genome_1", "chrom_1", 1, "G"),
                                                     Allele("genome_1", "chrom_1", 1, "A"),
                                                     Allele("genome_1", "chrom_1", 1, "G"),
                                                     Allele("genome_1", "chrom_1", 1, "A"),
                                                     Allele("genome_1", "chrom_1", 1, "C"),])
        actual = pangenome_variation.unique_allele_sequences
        expected=["A", "C", "G"]
        self.assertListEqual(actual, expected)

    def test___is_consistent___two_genomes___single_chrom_and_pos_for_each_genome(self):
        pangenome_variation = PangenomeVariation(0, [Allele("genome_1", "chrom_1", 1, "A"), Allele("genome_2", "chrom_2", 2, "A"),
                                                     Allele("genome_1", "chrom_1", 1, "A"), Allele("genome_2", "chrom_2", 2, "A"),
                                                     Allele("genome_1", "chrom_1", 1, "A"), Allele("genome_2", "chrom_2", 2, "A")])
        self.assertTrue(pangenome_variation.is_consistent())

    def test___is_consistent___two_genomes___different_chrom_for_genome_1(self):
        pangenome_variation = PangenomeVariation(0, [Allele("genome_1", "chrom_1", 1, "A"), Allele("genome_2", "chrom_2", 2, "A"),
                                                     Allele("genome_1", "chrom_2", 1, "A"), Allele("genome_2", "chrom_2", 2, "A"),
                                                     Allele("genome_1", "chrom_1", 1, "A"), Allele("genome_2", "chrom_2", 2, "A")])
        self.assertFalse(pangenome_variation.is_consistent())
    def test___is_consistent___two_genomes___different_chrom_for_genome_2(self):
        pangenome_variation = PangenomeVariation(0, [Allele("genome_1", "chrom_1", 1, "A"), Allele("genome_2", "chrom_2", 2, "A"),
                                                     Allele("genome_1", "chrom_1", 1, "A"), Allele("genome_2", "chrom_2", 2, "A"),
                                                     Allele("genome_1", "chrom_1", 1, "A"), Allele("genome_2", "chrom_1", 2, "A")])
        self.assertFalse(pangenome_variation.is_consistent())
    def test___is_consistent___two_genomes___different_pos_for_genome_1(self):
        pangenome_variation = PangenomeVariation(0, [Allele("genome_1", "chrom_1", 1, "A"), Allele("genome_2", "chrom_2", 2, "A"),
                                                     Allele("genome_1", "chrom_1", 1, "A"), Allele("genome_2", "chrom_2", 2, "A"),
                                                     Allele("genome_1", "chrom_1", 2, "A"), Allele("genome_2", "chrom_2", 2, "A")])
        self.assertFalse(pangenome_variation.is_consistent())
    def test___is_consistent___two_genomes___different_pos_for_genome_2(self):
        pangenome_variation = PangenomeVariation(0, [Allele("genome_1", "chrom_1", 1, "A"), Allele("genome_2", "chrom_2", 1, "A"),
                                                     Allele("genome_1", "chrom_1", 1, "A"), Allele("genome_2", "chrom_2", 2, "A"),
                                                     Allele("genome_1", "chrom_1", 1, "A"), Allele("genome_2", "chrom_2", 2, "A")])
        self.assertFalse(pangenome_variation.is_consistent())


class TestConsistentPangenomeVariations(TestCase):
    def setUp(self) -> None:
        self.dummy_consistent_pangenome_variations = ConsistentPangenomeVariations(PangenomeVariations())

    def test___constructor(self):
        # setup
        consistent_pangenome_variations = []
        alleles_to_consistent_pangenome_variations = {}
        for i in range(3):
            consistent_pangenome_variation = Mock()
            consistent_pangenome_variation.is_consistent.return_value = True
            consistent_pangenome_variation.alleles = [f"consistent_pangenome_variation_{i}.alleles"]
            alleles_to_consistent_pangenome_variations[f"consistent_pangenome_variation_{i}.alleles"] = consistent_pangenome_variation
            consistent_pangenome_variations.append(consistent_pangenome_variation)

        inconsistent_pangenome_variations = []
        for _ in range(3):
            inconsistent_pangenome_variation = Mock()
            inconsistent_pangenome_variation.is_consistent.return_value = False
            inconsistent_pangenome_variations.append(inconsistent_pangenome_variation)

        list_of_pangenome_variations = [
            consistent_pangenome_variations[0],
            consistent_pangenome_variations[1],
            inconsistent_pangenome_variations[0],
            inconsistent_pangenome_variations[1],
            consistent_pangenome_variations[2],
            inconsistent_pangenome_variations[2]
        ]
        pangenome_variations = PangenomeVariations()
        pangenome_variations._pangenome_variations = list_of_pangenome_variations
        actual_consistent_pangenome_variations = ConsistentPangenomeVariations(pangenome_variations)

        self.assertListEqual(actual_consistent_pangenome_variations.consistent_pangenome_variations, consistent_pangenome_variations)
        self.assertDictEqual(actual_consistent_pangenome_variations.alleles_to_consistent_pangenome_variations,
                             alleles_to_consistent_pangenome_variations)

    @patch.object(ConsistentPangenomeVariations, "alleles_to_consistent_pangenome_variations", new_callable=PropertyMock,
                  return_value=defaultdict(lambda: None, {"allele_1": "CPV1", "allele_2": "CPV1"}))
    def test___get_consistent_pangenome_variation___both_alleles_present_and_in_same_CPV(self, *mocks):
        pairwise_variation = Mock(allele_1="allele_1", allele_2="allele_2")
        actual = self.dummy_consistent_pangenome_variations.get_consistent_pangenome_variation(pairwise_variation)
        expected="CPV1"
        self.assertEqual(actual, expected)

    @patch.object(ConsistentPangenomeVariations, "alleles_to_consistent_pangenome_variations", new_callable=PropertyMock,
                  return_value=defaultdict(lambda: None, {"allele_1": "CPV1", "allele_2": "CPV2"}))
    def test___get_consistent_pangenome_variation___both_alleles_present_but_in_different_CPVs(self, *mocks):
        pairwise_variation = Mock(allele_1="allele_1", allele_2="allele_2")
        with self.assertRaises(InconsistentPangenomeVariations):
            self.dummy_consistent_pangenome_variations.get_consistent_pangenome_variation(pairwise_variation)

    @patch.object(ConsistentPangenomeVariations, "alleles_to_consistent_pangenome_variations", new_callable=PropertyMock,
                  return_value=defaultdict(lambda: None, {"allele_1": "CPV1"}))
    def test___get_consistent_pangenome_variation___only_first_allele_present(self, *mocks):
        pairwise_variation = Mock(allele_1="allele_1", allele_2="allele_2")
        with self.assertRaises(InconsistentPangenomeVariations):
            self.dummy_consistent_pangenome_variations.get_consistent_pangenome_variation(pairwise_variation)

    @patch.object(ConsistentPangenomeVariations, "alleles_to_consistent_pangenome_variations", new_callable=PropertyMock,
                  return_value=defaultdict(lambda: None, {"allele_2": "CPV1"}))
    def test___get_consistent_pangenome_variation___only_second_allele_present(self, *mocks):
        pairwise_variation = Mock(allele_1="allele_1", allele_2="allele_2")
        with self.assertRaises(InconsistentPangenomeVariations):
            self.dummy_consistent_pangenome_variations.get_consistent_pangenome_variation(pairwise_variation)

    @patch.object(ConsistentPangenomeVariations, "alleles_to_consistent_pangenome_variations", new_callable=PropertyMock,
                  return_value=defaultdict(lambda: None))
    def test___get_consistent_pangenome_variation___no_alleles_present(self, *mocks):
        pairwise_variation = Mock(allele_1="allele_1", allele_2="allele_2")
        actual = self.dummy_consistent_pangenome_variations.get_consistent_pangenome_variation(pairwise_variation)
        expected = None
        self.assertEqual(actual, expected)


    @patch.object(PangenomeVariation, "get_number_of_alleles", side_effect=[2,4])
    @patch.object(PangenomeVariation, "get_allele_index", side_effect=[1,0,2,3])
    @patch.object(PangenomeVariation, "get_number_of_different_allele_sequences", side_effect=[2,3])
    @patch.object(PangenomeVariation, "get_allele_sequence_index", side_effect=[0,0,2,1])
    @patch.object(ConsistentPangenomeVariations, "get_consistent_pangenome_variation",
                  side_effect=[None,PangenomeVariation(0, []), PangenomeVariation(1, []), None])
    @patch.object(PairwiseVariation, PairwiseVariation.get_PairwiseVariation_from_ShowSNPsDataframe.__name__,
                  return_value=[Mock(), Mock(), Mock(), Mock()])
    def test____enrich_ShowSNPsDataframe(self, *mocks):
        snps_df = pd.read_csv(StringIO(
"""dummy
0
1
2
3
"""
        ))
        actual = self.dummy_consistent_pangenome_variations._enrich_ShowSNPsDataframe("ref", "query", snps_df)

        expected=pd.read_csv(StringIO(
"""dummy,ref_genome,query_genome,present_in_a_consistent_pangenome_variation,pangenome_variation_id,number_of_alleles,ref_allele_id,query_allele_id,number_of_different_allele_sequences,ref_allele_sequence_id,query_allele_sequence_id
0,ref,query,False,-1,-1,-1,-1,-1,-1,-1
1,ref,query,True,0,2,1,0,2,0,0
2,ref,query,True,1,4,2,3,3,2,1
3,ref,query,False,-1,-1,-1,-1,-1,-1,-1
"""
        ))

        self.assertTrue(actual.equals(expected))


    @patch.object(DeduplicatePairwiseSNPsUtils, DeduplicatePairwiseSNPsUtils._get_ref_and_query_from_ShowSNPsDataframe_filepath.__name__,
                  return_value=(None, None))
    @patch.object(DeduplicatePairwiseSNPsUtils, DeduplicatePairwiseSNPsUtils._load_pickled_ShowSNPsDataframe.__name__)
    @patch.object(ConsistentPangenomeVariations, ConsistentPangenomeVariations._enrich_ShowSNPsDataframe.__name__)
    def test___load_and_process_ShowSNPsDataframe(self, _enrich_ShowSNPsDataframe_mock, *other_mocks):
        _enrich_ShowSNPsDataframe_mock.return_value = pd.read_csv(StringIO(
"""dummy,ref_genome,query_genome,present_in_a_consistent_pangenome_variation,pangenome_variation_id,number_of_alleles,ref_allele_id,query_allele_id,number_of_different_allele_sequences,ref_allele_sequence_id,query_allele_sequence_id
0,ref,query,False,-1,-1,-1,-1,-1,-1,-1
1,ref,query,True,0,2,1,0,2,0,0
2,ref,query,True,1,4,2,3,3,2,1
3,ref,query,False,-1,-1,-1,-1,-1,-1,-1
"""
        ))

        actual = self.dummy_consistent_pangenome_variations.load_and_process_ShowSNPsDataframe("dummy")
        expected=ShowSNPsDataframe(pd.read_csv(StringIO(
"""dummy,ref_genome,query_genome,present_in_a_consistent_pangenome_variation,pangenome_variation_id,number_of_alleles,ref_allele_id,query_allele_id,number_of_different_allele_sequences,ref_allele_sequence_id,query_allele_sequence_id
1,ref,query,True,0,2,1,0,2,0,0
2,ref,query,True,1,4,2,3,3,2,1
"""
        )))

        self.assertTrue(actual.equals(expected))



class TestDeduplicationGraph(TestCase):
    def setUp(self):
        self.alleles = [Allele(f"genome_{i}", "chrom_1", 1, "A") for i in range(10)]
        self.pairwise_mutations = [PairwiseVariation(self.alleles[i*2], self.alleles[i*2+1]) for i in range(5)]
        self.deduplication_graph = DeduplicationGraph()

    def test___add_variants_from_ShowSNPsDataframe_core___empty_df___no_variations_added(self):
        snps_df = pd.read_csv(StringIO(
"""ref_pos,ref_sub,query_sub,query_pos,nearest_mismatch,nearest_end,ref_len,query_len,ref_context,query_context,ref_strand,query_strand,ref_chrom,query_chrom
"""
        ))
        self.deduplication_graph._add_variants_from_ShowSNPsDataframe_core("genome_1", "genome_2", snps_df)
        graph_is_empty = self.deduplication_graph.graph.number_of_nodes()==0
        assert graph_is_empty

    def test___add_variants_from_ShowSNPsDataframe_core___one_variation_in_df___one_variation_added(self):
        snps_df = pd.read_csv(StringIO(
"""ref_pos,ref_sub,query_sub,query_pos,nearest_mismatch,nearest_end,ref_len,query_len,ref_context,query_context,ref_strand,query_strand,ref_chrom,query_chrom
1,A,C,2,0,0,0,0,ACGT,ACGT,1,1,chrom_1,chrom_2
"""
        ))
        self.deduplication_graph._add_variants_from_ShowSNPsDataframe_core("genome_1", "genome_2", snps_df)

        pairwise_variation = PairwiseVariation(Allele("genome_1", "chrom_1", 1, "A"),
                                               Allele("genome_2", "chrom_2", 2, "C"))
        assert self.deduplication_graph.graph.number_of_nodes()==1 and self.deduplication_graph.graph.has_node(pairwise_variation)


    def test___add_variants_from_ShowSNPsDataframe_core___one_SNP_variation_four_non_SNP_variation_in_df___only_SNP_variation_added(self):
        snps_df = pd.read_csv(StringIO(
"""ref_pos,ref_sub,query_sub,query_pos,nearest_mismatch,nearest_end,ref_len,query_len,ref_context,query_context,ref_strand,query_strand,ref_chrom,query_chrom
1,,C,2,0,0,0,0,ACGT,ACGT,1,1,chrom_1,chrom_2
1,A,CC,2,0,0,0,0,ACGT,ACGT,1,1,chrom_1,chrom_2
1,A,C,2,0,0,0,0,ACGT,ACGT,1,1,chrom_1,chrom_2
1,AA,CC,2,0,0,0,0,ACGT,ACGT,1,1,chrom_1,chrom_2
1,AA,C,2,0,0,0,0,ACGT,ACGT,1,1,chrom_1,chrom_2
"""
        ), na_filter=False)
        self.deduplication_graph._add_variants_from_ShowSNPsDataframe_core("genome_1", "genome_2", snps_df)

        pairwise_variation = PairwiseVariation(Allele("genome_1", "chrom_1", 1, "A"),
                                               Allele("genome_2", "chrom_2", 2, "C"))
        assert self.deduplication_graph.graph.number_of_nodes()==1 and self.deduplication_graph.graph.has_node(pairwise_variation)



    def test___add_variants_from_ShowSNPsDataframe_core___three_SNP_variation_four_non_SNP_variation_in_df___only_SNPs_variations_added(self):
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
        self.deduplication_graph._add_variants_from_ShowSNPsDataframe_core("genome_1", "genome_2", snps_df)


        pairwise_variation_1 = PairwiseVariation(Allele("genome_1", "chrom_10", 10, "G"),
                                                 Allele("genome_2", "chrom_20", 20, "T"))
        pairwise_variation_2 = PairwiseVariation(Allele("genome_1", "chrom_1", 1, "A"),
                                                 Allele("genome_2", "chrom_2", 2, "C"))
        pairwise_variation_3 = PairwiseVariation(Allele("genome_1", "chrom_100", 100, "A"),
                                                 Allele("genome_2", "chrom_200", 200, "T"))
        assert self.deduplication_graph.graph.number_of_nodes() == 3 and \
               self.deduplication_graph.graph.has_node(pairwise_variation_1) and \
               self.deduplication_graph.graph.has_node(pairwise_variation_2) and \
               self.deduplication_graph.graph.has_node(pairwise_variation_3)


    def test___add_variants_from_ShowSNPsDataframe_core___one_same_variation_in_three_dfs___one_variation_added(self):
        snps_df = pd.read_csv(StringIO(
"""ref_pos,ref_sub,query_sub,query_pos,nearest_mismatch,nearest_end,ref_len,query_len,ref_context,query_context,ref_strand,query_strand,ref_chrom,query_chrom
1,A,C,2,0,0,0,0,ACGT,ACGT,1,1,chrom_1,chrom_2
"""
        ))
        self.deduplication_graph._add_variants_from_ShowSNPsDataframe_core("genome_1", "genome_2", snps_df)
        self.deduplication_graph._add_variants_from_ShowSNPsDataframe_core("genome_1", "genome_2", snps_df)
        self.deduplication_graph._add_variants_from_ShowSNPsDataframe_core("genome_1", "genome_2", snps_df)

        pairwise_variation = PairwiseVariation(Allele("genome_1", "chrom_1", 1, "A"),
                                               Allele("genome_2", "chrom_2", 2, "C"))
        assert self.deduplication_graph.graph.number_of_nodes()==1 and self.deduplication_graph.graph.has_node(pairwise_variation)

    def test___add_variants_from_ShowSNPsDataframe_core___one_variation_in_three_dfs___three_variations_added(self):
        snps_df = pd.read_csv(StringIO(
"""ref_pos,ref_sub,query_sub,query_pos,nearest_mismatch,nearest_end,ref_len,query_len,ref_context,query_context,ref_strand,query_strand,ref_chrom,query_chrom
1,A,C,2,0,0,0,0,ACGT,ACGT,1,1,chrom_1,chrom_2
"""
        ))
        self.deduplication_graph._add_variants_from_ShowSNPsDataframe_core("genome_1", "genome_2", snps_df)
        self.deduplication_graph._add_variants_from_ShowSNPsDataframe_core("genome_3", "genome_2", snps_df)
        self.deduplication_graph._add_variants_from_ShowSNPsDataframe_core("genome_4", "genome_2", snps_df)

        pairwise_variation_1 = PairwiseVariation(Allele("genome_1", "chrom_1", 1, "A"),
                                                 Allele("genome_2", "chrom_2", 2, "C"))
        pairwise_variation_2 = PairwiseVariation(Allele("genome_3", "chrom_1", 1, "A"),
                                                 Allele("genome_2", "chrom_2", 2, "C"))
        pairwise_variation_3 = PairwiseVariation(Allele("genome_4", "chrom_1", 1, "A"),
                                                 Allele("genome_2", "chrom_2", 2, "C"))
        assert self.deduplication_graph.graph.number_of_nodes()==3 and \
               self.deduplication_graph.graph.has_node(pairwise_variation_1) and \
               self.deduplication_graph.graph.has_node(pairwise_variation_2) and \
               self.deduplication_graph.graph.has_node(pairwise_variation_3)

    @patch.object(PairwiseVariation, PairwiseVariation.share_allele.__name__, return_value=False)
    @patch.object(DeduplicationGraph, "nodes", new_callable=PropertyMock)
    def test___build_edges___no_shared_alleles___no_edge_built(self, nodes_mock, share_allele_mock):
        nodes_mock.return_value = self.pairwise_mutations
        self.deduplication_graph.build_edges()
        self.assertEqual(self.deduplication_graph.graph.number_of_edges(), 0)

    @patch.object(PairwiseVariation, PairwiseVariation.share_allele.__name__, autospec=True)
    @patch.object(DeduplicationGraph, "nodes", new_callable=PropertyMock)
    def test___build_edges___four_shared_alleles___four_edges_built(self, nodes_mock, share_allele_mock):
        # set up
        nodes_mock.return_value = self.pairwise_mutations
        # create a random matrix with 3 variations with shared alleles
        allele_sharing_matrix = defaultdict(lambda: defaultdict(int))
        allele_sharing_matrix[self.pairwise_mutations[0]][self.pairwise_mutations[2]] = 1  # add both ways
        allele_sharing_matrix[self.pairwise_mutations[2]][self.pairwise_mutations[0]] = 1  # add both ways
        allele_sharing_matrix[self.pairwise_mutations[1]][self.pairwise_mutations[0]] = 1  # add one way
        allele_sharing_matrix[self.pairwise_mutations[3]][self.pairwise_mutations[3]] = 1  # self loop, not added
        allele_sharing_matrix[self.pairwise_mutations[4]][self.pairwise_mutations[2]] = 1  # add one way
        def share_allele_return_mock(variant_1, variant_2):
            return allele_sharing_matrix[variant_1][variant_2]
        share_allele_mock.side_effect = share_allele_return_mock

        self.deduplication_graph.build_edges()
        self.assertTrue(self.deduplication_graph.graph.has_edge(self.pairwise_mutations[0], self.pairwise_mutations[2]))
        self.assertTrue(self.deduplication_graph.graph.has_edge(self.pairwise_mutations[2], self.pairwise_mutations[0]))
        self.assertTrue(self.deduplication_graph.graph.has_edge(self.pairwise_mutations[1], self.pairwise_mutations[0]))
        self.assertTrue(self.deduplication_graph.graph.has_edge(self.pairwise_mutations[0], self.pairwise_mutations[1]))
        self.assertTrue(self.deduplication_graph.graph.has_edge(self.pairwise_mutations[4], self.pairwise_mutations[2]))
        self.assertTrue(self.deduplication_graph.graph.has_edge(self.pairwise_mutations[2], self.pairwise_mutations[4]))
        self.assertEqual(self.deduplication_graph.graph.number_of_edges(), 3)


    @patch.object(DeduplicationGraph, DeduplicationGraph._get_connected_components.__name__)
    def test___get_pangenome_variations___totally_unconnected_graph(self, get_connected_components_mock):
        get_connected_components_mock.return_value = [
            {self.pairwise_mutations[0]},
            {self.pairwise_mutations[1]},
            {self.pairwise_mutations[2]},
            {self.pairwise_mutations[3]},
            {self.pairwise_mutations[4]}
        ]
        pangenome_variations_actual = self.deduplication_graph.get_pangenome_variations()

        pangenome_variations_expected = PangenomeVariations()
        for i in range(5):
            pangenome_variations_expected.append(PangenomeVariation(i, [
                self.pairwise_mutations[i].allele_1, self.pairwise_mutations[i].allele_2]))

        self.assertEqual(pangenome_variations_actual, pangenome_variations_expected)


    @patch.object(DeduplicationGraph, DeduplicationGraph._get_connected_components.__name__)
    def test___get_pangenome_variations___totally_connected_graph(self, get_connected_components_mock):
        get_connected_components_mock.return_value = [
            {self.pairwise_mutations[0], self.pairwise_mutations[1], self.pairwise_mutations[2],
             self.pairwise_mutations[3], self.pairwise_mutations[4]}
        ]
        pangenome_variations_actual = self.deduplication_graph.get_pangenome_variations()

        pangenome_variations_expected = PangenomeVariations()
        pangenome_variations_expected.append(PangenomeVariation(0, self.alleles))

        self.assertEqual(pangenome_variations_actual, pangenome_variations_expected)

    @patch.object(DeduplicationGraph, DeduplicationGraph._get_connected_components.__name__)
    def test___get_pangenome_variations___isolated_node_and_a_path(self, get_connected_components_mock):
        get_connected_components_mock.return_value = [
            {self.pairwise_mutations[0], self.pairwise_mutations[1], self.pairwise_mutations[3], self.pairwise_mutations[4]},
            {self.pairwise_mutations[2]}
        ]
        pangenome_variations_actual = self.deduplication_graph.get_pangenome_variations()

        pangenome_variations_expected = PangenomeVariations()
        pangenome_variations_expected.append(PangenomeVariation(0, [
            self.pairwise_mutations[0].allele_1, self.pairwise_mutations[0].allele_2,
            self.pairwise_mutations[1].allele_1, self.pairwise_mutations[1].allele_2,
            self.pairwise_mutations[3].allele_1, self.pairwise_mutations[3].allele_2,
            self.pairwise_mutations[4].allele_1, self.pairwise_mutations[4].allele_2,
        ]))
        pangenome_variations_expected.append(PangenomeVariation(1, [
            self.pairwise_mutations[2].allele_1, self.pairwise_mutations[2].allele_2,
        ]))

        self.assertEqual(pangenome_variations_actual, pangenome_variations_expected)

    @patch.object(DeduplicationGraph, DeduplicationGraph._get_connected_components.__name__)
    def test___get_pangenome_variations___two_and_three_sized_components(self, get_connected_components_mock):
        get_connected_components_mock.return_value = [
            {self.pairwise_mutations[0], self.pairwise_mutations[1], self.pairwise_mutations[3]},
            {self.pairwise_mutations[2], self.pairwise_mutations[4]}
        ]
        pangenome_variations_actual = self.deduplication_graph.get_pangenome_variations()

        pangenome_variations_expected = PangenomeVariations()
        pangenome_variations_expected.append(PangenomeVariation(0, [
            self.pairwise_mutations[0].allele_1, self.pairwise_mutations[0].allele_2,
            self.pairwise_mutations[1].allele_1, self.pairwise_mutations[1].allele_2,
            self.pairwise_mutations[3].allele_1, self.pairwise_mutations[3].allele_2,
        ]))
        pangenome_variations_expected.append(PangenomeVariation(1, [
            self.pairwise_mutations[2].allele_1, self.pairwise_mutations[2].allele_2,
            self.pairwise_mutations[4].allele_1, self.pairwise_mutations[4].allele_2,
        ]))

        self.assertEqual(pangenome_variations_actual, pangenome_variations_expected)