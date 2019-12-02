from evaluate.vcf import PandoraVCF, BuggedVCFError, NullVCFError, VCFFactory
from .common import retrieve_entry_from_test_vcf
from unittest.mock import patch, PropertyMock, Mock, MagicMock
import pytest

def build_PandoraVCF_bypassing_check(variant=None, sample=None):
    vcf = PandoraVCF()
    vcf.variant = variant
    vcf.sample = sample
    return vcf


class Test_PandoraVCF:
    @patch.object(PandoraVCF, "is_null_call", new_callable=PropertyMock, return_value=True)
    def test___from_VariantRecord_and_Sample___null_PandoraVCF_raises_NullVCFError(self, *mocks):
        with pytest.raises(NullVCFError):
            VCFFactory.create_Pandora_VCF_from_VariantRecord_and_Sample()

    @patch.object(PandoraVCF, "is_null_call", new_callable=PropertyMock, return_value=False)
    @patch.object(PandoraVCF, "has_genotype_bug", new_callable=PropertyMock, return_value=True)
    def test___from_VariantRecord_and_Sample___non_null_PandoraVCF_but_has_genotype_bug_raises_BuggedVCFError(self, *mocks):
        with pytest.raises(BuggedVCFError):
            VCFFactory.create_Pandora_VCF_from_VariantRecord_and_Sample()

    @patch.object(PandoraVCF, "is_null_call", new_callable=PropertyMock, return_value=False)
    @patch.object(PandoraVCF, "has_genotype_bug", new_callable=PropertyMock, return_value=False)
    def test___from_VariantRecord_and_Sample___valid_PandoraVCF(self, *mocks):
        variant_mock = Mock()
        sample_mock = Mock()
        vcf = VCFFactory.create_Pandora_VCF_from_VariantRecord_and_Sample(variant_mock, sample_mock)

        assert vcf.variant == variant_mock
        assert vcf.sample == sample_mock


    def test___genotype___one_gt(self):
        sample_name = "sample"
        sample_mock = MagicMock(get=MagicMock(return_value=[2]))
        variant_mock = MagicMock(samples={sample_name: sample_mock})
        vcf = build_PandoraVCF_bypassing_check(variant=variant_mock, sample=sample_name)

        actual = vcf.genotype
        expected = 2

        sample_mock.get.assert_called_once_with("GT")
        assert actual == expected

    @pytest.mark.xfail(strict=True)
    def test___genotype___no_gts___expects_death(self):
        sample_name = "sample"
        sample_mock = MagicMock(get=MagicMock(return_value=[]))
        variant_mock = MagicMock(samples={sample_name: sample_mock})
        vcf = build_PandoraVCF_bypassing_check(variant=variant_mock, sample=sample_name)

        actual = vcf.genotype

    @pytest.mark.xfail(strict=True)
    def test___genotype___two_gts___expects_death(self):
        sample_name = "sample"
        sample_mock = MagicMock(get=MagicMock(return_value=[1,2]))
        variant_mock = MagicMock(samples={sample_name: sample_mock})
        vcf = build_PandoraVCF_bypassing_check(variant=variant_mock, sample=sample_name)

        actual = vcf.genotype

    def test___genotype___genotype_is_one___returns_one___real_files_test(self):
        entry = retrieve_entry_from_test_vcf(1)
        sample = "sample"
        vcf = build_PandoraVCF_bypassing_check(entry, sample)

        actual = vcf.genotype
        expected = 1

        assert actual == expected

    @patch.object(PandoraVCF, "genotype", new_callable=PropertyMock, return_value=None)
    def test___is_null_call___with_None_genotype___returns_true(self, *mocks):
        vcf = PandoraVCF()

        actual = vcf.is_null_call
        expected = True

        assert actual == expected

    @patch.object(PandoraVCF, "genotype", new_callable=PropertyMock, return_value=0)
    def test___is_null_call___with_genotype_0___returns_false(self, *mocks):
        vcf = PandoraVCF()

        actual = vcf.is_null_call
        expected = False

        assert actual == expected

    @patch.object(
        PandoraVCF,
        "_highest_likelihood_indexes",
        new_callable=PropertyMock,
        return_value=[0, 2],
    )
    @patch.object(PandoraVCF, "genotype", new_callable=PropertyMock, return_value=1)
    def test___has_genotype_bug___with_wrongly_called_genotype___returns_true(
        self, *mocks
    ):
        vcf = PandoraVCF()

        actual = vcf.has_genotype_bug
        expected = True

        assert actual == expected

    @patch.object(
        PandoraVCF,
        "_highest_likelihood_indexes",
        new_callable=PropertyMock,
        return_value=[1],
    )
    @patch.object(PandoraVCF, "genotype", new_callable=PropertyMock, return_value=1)
    def test___has_genotype_bug___correctly_called_genotype___returns_false(
        self, *mocks
    ):
        vcf = PandoraVCF()

        actual = vcf.has_genotype_bug
        expected = False

        assert actual == expected

    def test___genotype_confidence(self):
        sample_name = "sample"
        sample_mock = MagicMock(get=MagicMock(return_value=262.757))
        variant_mock = MagicMock(samples={sample_name: sample_mock})
        vcf = build_PandoraVCF_bypassing_check(variant=variant_mock, sample=sample_name)

        actual = vcf.genotype_confidence
        expected = 262.757

        sample_mock.get.assert_called_once_with("GT_CONF")
        assert actual == expected


    def test___called_variant_sequence___genotype_is_zero___returns_ref(self):
        entry = retrieve_entry_from_test_vcf(2)
        sample = "sample"
        vcf = build_PandoraVCF_bypassing_check(entry, sample)

        actual = vcf.called_variant_sequence
        expected = "CTGCCCGTTGGC"

        assert actual == expected


    def test___called_variant_sequence___genotype_is_one___returns_first_alt(self):
        entry = retrieve_entry_from_test_vcf(1)
        sample = "sample"
        vcf = build_PandoraVCF_bypassing_check(entry, sample)

        actual = vcf.called_variant_sequence
        expected = "TTGGGGGAAGGCTCTGCACTGCCCGTTGGC"

        assert actual == expected


    def test___called_variant_length___genotype_is_zero___returns_length_ref(self):
        entry = retrieve_entry_from_test_vcf(2)
        sample = "sample"
        vcf = build_PandoraVCF_bypassing_check(entry, sample)

        actual = vcf.called_variant_length
        expected = 12

        assert actual == expected


    def test___called_variant_length___genotype_is_one___returns_length_first_alt(self):
        entry = retrieve_entry_from_test_vcf(1)
        sample = "sample"
        vcf = build_PandoraVCF_bypassing_check(entry, sample)

        actual = vcf.called_variant_length
        expected = 30

        assert actual == expected


    def test___svtype(self):
        entry = retrieve_entry_from_test_vcf(0)
        sample = "sample"
        vcf = build_PandoraVCF_bypassing_check(entry, sample)

        actual = vcf.svtype
        expected = "COMPLEX"

        assert actual == expected


    def test_meanCoverageForward(self):
        entry = retrieve_entry_from_test_vcf(2)
        sample = "sample"
        vcf = build_PandoraVCF_bypassing_check(entry, sample)

        actual = vcf._mean_coverage_forward
        expected = 24

        assert actual == expected

    def test_meanCoverageReverse(self):
        entry = retrieve_entry_from_test_vcf(1)
        sample = "sample"
        vcf = build_PandoraVCF_bypassing_check(entry, sample)

        actual = vcf._mean_coverage_reverse
        expected = 7

        assert actual == expected

    def test_meanCoverage(self):
        entry = retrieve_entry_from_test_vcf(1)
        sample = "sample"
        vcf = build_PandoraVCF_bypassing_check(entry, sample)

        actual = vcf.coverage
        expected = 13

        assert actual == expected


    def test___likelihoods___from_PandoraVCF_file(self):
        entry = retrieve_entry_from_test_vcf(1)
        sample = "sample"
        vcf = build_PandoraVCF_bypassing_check(entry, sample)

        actual = vcf._likelihoods
        expected = [-63.3221, -326.079, -432.546]

        assert actual == expected

    @patch.object(
        PandoraVCF, "_likelihoods", new_callable=PropertyMock, return_value=[-100, -200, -2]
    )
    def test___highest_likelihood_indexes___one_highest_index(self, *mocks):
        vcf = build_PandoraVCF_bypassing_check()

        actual = vcf._highest_likelihood_indexes
        expected = [2]

        assert actual == expected


    @patch.object(
        PandoraVCF, "_likelihoods", new_callable=PropertyMock, return_value=[-2, -200, -2]
    )
    def test___highest_likelihood_indexes___two_highest_indexes(self, *mocks):
        vcf = build_PandoraVCF_bypassing_check()

        actual = vcf._highest_likelihood_indexes
        expected = [0, 2]

        assert actual == expected


    def test___gaps___from_PandoraVCF_file(self):
        entry = retrieve_entry_from_test_vcf(1)
        sample = "sample"
        vcf = build_PandoraVCF_bypassing_check(entry, sample)

        actual = vcf._gaps
        expected = 0.75

        assert actual == expected
