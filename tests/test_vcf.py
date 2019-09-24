from evaluate.vcf import VCF
from .common import retrieve_entry_from_test_vcf
from unittest.mock import patch, PropertyMock, MagicMock


class Test_VCF:
    def test_init_retrieveVCFFromFile(self):
        entry = retrieve_entry_from_test_vcf(1)
        sample = "sample"

        actual = VCF(entry, sample)
        expected = MagicMock(
            genotype=1,
            genotype_confidence=262.757,
            variant_sequence="TTGGGGGAAGGCTCTGCACTGCCCGTTGGC",
            svtype="COMPLEX",
            mean_coverage_forward=6,
            mean_coverage_reverse=7,
            pos=1,
            chrom="GC00000001_155",
        )

        assert actual == expected

    @patch.object(VCF, "genotype", new_callable=PropertyMock, return_value=None)
    def test_isInvalidVcfEntry_withNoneGenotype_returnTrue(self, *mocks):
        vcf = VCF()

        actual = vcf.is_invalid_vcf_entry
        expected = True

        assert actual == expected

    @patch.object(
        VCF,
        "highest_likelihood_indexes",
        new_callable=PropertyMock,
        return_value=[0, 2],
    )
    @patch.object(VCF, "genotype", new_callable=PropertyMock, return_value=1)
    def test_isInvalidVcfEntry_withWronglyCalledGenotypeGenotype_returnTrue(
        self, *mocks
    ):
        vcf = VCF()

        actual = vcf.is_invalid_vcf_entry
        expected = True

        assert actual == expected

    @patch.object(
        VCF, "highest_likelihood_indexes", new_callable=PropertyMock, return_value=[1]
    )
    @patch.object(VCF, "genotype", new_callable=PropertyMock, return_value=1)
    def test_isInvalidVcfEntry_withGenotype1_returnFalse(self, *mocks):
        vcf = VCF()

        actual = vcf.is_invalid_vcf_entry
        expected = False

        assert actual == expected

    def test_genotypeConfidence(self):
        sample_name = "sample"
        sample_mock = MagicMock(get=MagicMock(return_value=262.757))
        variant_mock = MagicMock(samples={sample_name: sample_mock})
        vcf = VCF(variant=variant_mock, sample=sample_name)

        actual = vcf.genotype_confidence
        expected = 262.757

        sample_mock.get.assert_called_once_with("GT_CONF", 0)
        assert actual == expected

    def test_svtype(self):
        entry = retrieve_entry_from_test_vcf(0)
        sample = "sample"
        vcf = VCF(entry, sample)

        actual = vcf.svtype
        expected = "COMPLEX"

        assert actual == expected

    def test_meanCoverageForward(self):
        entry = retrieve_entry_from_test_vcf(2)
        sample = "sample"
        vcf = VCF(entry, sample)

        actual = vcf.mean_coverage_forward
        expected = 24

        assert actual == expected

    def test_meanCoverageReverse(self):
        entry = retrieve_entry_from_test_vcf(1)
        sample = "sample"
        vcf = VCF(entry, sample)

        actual = vcf.mean_coverage_reverse
        expected = 7

        assert actual == expected

    def test_meanCoverage(self):
        entry = retrieve_entry_from_test_vcf(1)
        sample = "sample"
        vcf = VCF(entry, sample)

        actual = vcf.mean_coverage
        expected = 13

        assert actual == expected

    def test_genotype_genotypeNone_returnNone(self):
        entry = retrieve_entry_from_test_vcf(0)
        sample = "sample"
        vcf = VCF(entry, sample)

        actual = vcf.genotype
        expected = None

        assert actual == expected

    def test_genotype_genotype1_return1(self):
        entry = retrieve_entry_from_test_vcf(1)
        sample = "sample"
        vcf = VCF(entry, sample)

        actual = vcf.genotype
        expected = 1

        assert actual == expected

    def test_variantSequence_genotypeNone_returnRef(self):
        entry = retrieve_entry_from_test_vcf(0)
        sample = "sample"
        vcf = VCF(entry, sample)

        actual = vcf.variant_sequence
        expected = "CTGCCCGTTGGC"

        assert actual == expected

    def test_variantSequence_genotypeOne_returnFirstAlt(self):
        entry = retrieve_entry_from_test_vcf(1)
        sample = "sample"
        vcf = VCF(entry, sample)

        actual = vcf.variant_sequence
        expected = "TTGGGGGAAGGCTCTGCACTGCCCGTTGGC"

        assert actual == expected

    def test_variantSequence_genotypeZero_returnRef(self):
        entry = retrieve_entry_from_test_vcf(2)
        sample = "sample"
        vcf = VCF(entry, sample)

        actual = vcf.variant_sequence
        expected = "CTGCCCGTTGGC"

        assert actual == expected

    def test_variantLength_genotypeNone_returnRef(self):
        entry = retrieve_entry_from_test_vcf(0)
        sample = "sample"
        vcf = VCF(entry, sample)

        actual = vcf.variant_length
        expected = 12

        assert actual == expected

    def test_variantLength_genotypeOne_returnFirstAlt(self):
        entry = retrieve_entry_from_test_vcf(1)
        sample = "sample"
        vcf = VCF(entry, sample)

        actual = vcf.variant_length
        expected = 30

        assert actual == expected

    def test_variantLength_genotypeZero_returnRef(self):
        entry = retrieve_entry_from_test_vcf(2)
        sample = "sample"
        vcf = VCF(entry, sample)

        actual = vcf.variant_length
        expected = 12

        assert actual == expected

    def test_likelihoods_fromVCFFile(self):
        entry = retrieve_entry_from_test_vcf(1)
        sample = "sample"
        vcf = VCF(entry, sample)

        actual = vcf.likelihoods
        expected = [-63.3221, -326.079, -432.546]

        assert actual == expected

    @patch.object(
        VCF, "likelihoods", new_callable=PropertyMock, return_value=[-100, -200, -2]
    )
    def test_highestLikelihoodIndexes_oneHighestIndex(self, *mocks):
        vcf = VCF()

        actual = vcf.highest_likelihood_indexes
        expected = [2]

        assert actual == expected

    def test_gaps_fromVCFFile(self):
        entry = retrieve_entry_from_test_vcf(1)
        sample = "sample"
        vcf = VCF(entry, sample)

        actual = vcf.gaps
        expected = 0.75

        assert actual == expected
