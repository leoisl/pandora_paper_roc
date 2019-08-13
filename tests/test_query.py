from evaluate.query import *
import pytest
from pathlib import Path

TEST_CASES = Path("tests/test_cases")
TEST_VCF = TEST_CASES / "test.vcf"
TEST_PANEL = TEST_CASES / "test_panel.fa"
TEST_REF_SEQ = TEST_CASES / "test_reference.fa"
TEST_TMP_PANEL = "/tmp/deleteme.fa"
TEST_MAKE_PROBE_VCF = TEST_CASES / "test_make_probe.vcf"
TEST_QUERY_VCF = TEST_CASES / "test_query.vcf"
TEST_QUERY_REF = TEST_CASES / "test_query.fa"


def retrieve_entry_from_test_vcf(idx: int) -> pysam.VariantRecord:
    with pysam.VariantFile(TEST_VCF) as vcf:
        for i, record in enumerate(vcf):
            if i == idx:
                return record
    raise IndexError("You asked for an index that is beyond the number in the test VCF")


def retrieve_entry_from_test_query_vcf(idx: int) -> pysam.VariantRecord:
    with pysam.VariantFile(TEST_QUERY_VCF) as vcf:
        for i, record in enumerate(vcf):
            if i == idx:
                return record
    raise IndexError("You asked for an index that is beyond the number in the test VCF")


def create_sam_header(name: str, length: int) -> pysam.AlignmentHeader:
    return pysam.AlignmentHeader.from_text(
        f"@SQ	SN:{name}	LN:{length}\n@PG	ID:bwa	PN:bwa	VN:0.7.17-r1188	CL:bwa mem -t 1 panel.fa -"
    )


class TestQuery:
    def test_calculateProbeBoundariesForEntry_variantShorterThanMinLen_returnProbeOfMinLen(
        self
    ):
        flank_width = 10
        sample = "sample"
        query = Query(TEST_QUERY_VCF, TEST_PANEL, [sample], flank_width=flank_width)
        variant = retrieve_entry_from_test_query_vcf(1)

        expected = (9, 32)
        actual = query.calculate_probe_boundaries_for_entry(variant)

        assert actual == expected

    def test_calculateProbeBoundariesForEntry_variantAtStartOfGene_returnZeroLenLeftProbe(
        self
    ):
        flank_width = 10
        sample = "sample"
        query = Query(TEST_QUERY_VCF, TEST_PANEL, [sample], flank_width=flank_width)
        variant = retrieve_entry_from_test_query_vcf(0)

        expected = (0, 13)
        actual = query.calculate_probe_boundaries_for_entry(variant)

        assert actual == expected

    def test_calculateProbeBoundariesForEntry_variantLongerThanMinLen_returnZeroLenProbes(
        self
    ):
        flank_width = 2
        sample = "sample"
        query = Query(TEST_QUERY_VCF, TEST_PANEL, [sample], flank_width=flank_width)
        variant = retrieve_entry_from_test_query_vcf(1)

        expected = (17, 24)
        actual = query.calculate_probe_boundaries_for_entry(variant)

        assert actual == expected

    def test_createProbesForGeneVariants_emptyVariants_returnEmptyProbes(self):
        samples = ["sample"]
        query = Query(TEST_QUERY_VCF, TEST_QUERY_REF, samples)

        expected = {sample: "" for sample in samples}
        actual = query._create_probes_for_gene_variants(pysam.FastxRecord(), [])

        assert actual == expected

    def test_makeProbes_emptyVariantsReturnsEmptyProbes(self):
        vcf = TEST_CASES / "empty.vcf"
        genes = TEST_QUERY_REF
        samples = ["sample"]
        query = Query(vcf, genes, samples)

        actual = query.make_probes()
        expected = {s: "" for s in samples}

        assert actual == expected

    def test_makeProbes_emptyGenesReturnsEmptyProbes(self):
        vcf = TEST_CASES / "empty.vcf"
        genes = TEST_CASES / "empty.fa"
        samples = ["sample"]
        query = Query(vcf, genes, samples)

        actual = query.make_probes()
        expected = {s: "" for s in samples}

        assert actual == expected

    def test_makeProbes_oneGeneOneVcfRecordNotInGeneReturnsEmptyProbes(self):
        vcf = TEST_CASES / "empty.vcf"
        genes = TEST_QUERY_REF
        samples = ["sample"]
        query = Query(vcf, genes, samples)

        actual = query.make_probes()
        expected = {s: "" for s in samples}

        assert actual == expected

    def test_makeProbes_oneGeneOneVcfRecordInGeneReturnsOneProbe(self):
        vcf = TEST_CASES / "make_probes_1.vcf"
        genes = TEST_CASES / "make_probes_1.fa"
        flank_width = 3
        samples = ["sample"]
        query = Query(vcf, genes, samples, flank_width)

        actual = query.make_probes()
        expected = {
            "sample": (
                ">gene1_SAMPLE=sample_POS=4_INTERVAL=(0,7)_SVTYPE=COMPLEX_MEAN_"
                "FWD_COVG=6_MEAN_REV_COVG=7_GT_CONF=262.757\nxxxFxxx\n"
            )
        }

        assert actual == expected

    @pytest.mark.xfail(reason="We are ignoring overlapping reads at the moment")
    def test_makeProbes_oneGeneTwoOverlappingVcfRecordsInGeneRaisesException(self):
        vcf = TEST_CASES / "make_probes_2.vcf"
        genes = TEST_CASES / "make_probes_1.fa"
        min_probe_length = 3
        query = Query(vcf, genes, min_probe_length)

        with pytest.raises(OverlappingRecordsError):
            query.make_probes()

    @pytest.mark.xfail(reason="We are ignoring overlapping reads at the moment")
    def test_makeProbes_oneGeneTwoCloseVcfRecordsInGeneReturnsOneProbe(self):
        vcf = TEST_CASES / "make_probes_3.vcf"
        genes = TEST_CASES / "make_probes_2.fa"
        min_probe_length = 7
        query = Query(vcf, genes, min_probe_length)

        actual = query.make_probes()
        expected = ">gene1_interval=(1, 11)\nxxFOOxxFOOxx\n"

        assert actual == expected

    def test_makeProbes_oneGeneTwoNonCloseVcfRecordsInGeneReturnsTwoProbes(self):
        vcf = TEST_CASES / "make_probes_3.vcf"
        genes = TEST_CASES / "make_probes_2.fa"
        flank_width = 5
        samples = ["sample"]
        query = Query(vcf, genes, samples, flank_width)

        actual = query.make_probes()
        expected = {
            "sample": (
                ">gene1_SAMPLE=sample_POS=4_INTERVAL=(0,10)_SVTYPE=COMPLEX_MEAN_"
                "FWD_COVG=6_MEAN_REV_COVG=7_GT_CONF=262.757\nxxxFOOxxxxx\n"
                ">gene1_SAMPLE=sample_POS=8_INTERVAL=(2,14)_SVTYPE=COMPLEX_MEAN_"
                "FWD_COVG=6_MEAN_REV_COVG=7_GT_CONF=262.757\nxxxxxFOOxxxxx\n"
            )
        }

        assert actual == expected

    def test_makeProbes_twoGenesTwoNonCloseVcfRecordsInOneGeneReturnsTwoProbes(self):
        vcf = TEST_CASES / "make_probes_3.vcf"
        genes = TEST_CASES / "make_probes_3.fa"
        flank_width = 5
        samples = ["sample"]
        query = Query(vcf, genes, samples, flank_width)

        actual = query.make_probes()
        expected = {
            "sample": (
                ">gene1_SAMPLE=sample_POS=4_INTERVAL=(0,10)_SVTYPE=COMPLEX_MEAN_"
                "FWD_COVG=6_MEAN_REV_COVG=7_GT_CONF=262.757\nxxxFOOxxxxx\n"
                ">gene1_SAMPLE=sample_POS=8_INTERVAL=(2,14)_SVTYPE=COMPLEX_MEAN_"
                "FWD_COVG=6_MEAN_REV_COVG=7_GT_CONF=262.757\nxxxxxFOOxxx\n"
            )
        }

        assert actual == expected

    def test_makeProbes_twoGenesTwoVcfRecordsOneInEachGeneReturnsTwoProbes(self):
        vcf = TEST_CASES / "make_probes_4.vcf"
        genes = TEST_CASES / "make_probes_3.fa"
        flank_width = 5
        samples = ["sample"]
        query = Query(vcf, genes, samples, flank_width)

        actual = query.make_probes()
        expected = {
            "sample": (
                ">gene1_SAMPLE=sample_POS=4_INTERVAL=(0,10)_SVTYPE=COMPLEX_MEAN_"
                "FWD_COVG=6_MEAN_REV_COVG=7_GT_CONF=262.757\nxxxFOOxxxxx\n"
                ">gene2_SAMPLE=sample_POS=2_INTERVAL=(0,8)_SVTYPE=COMPLEX_MEAN_"
                "FWD_COVG=6_MEAN_REV_COVG=7_GT_CONF=262.757\nxFOOx\n"
            )
        }

        assert actual == expected

    def test_createProbeHeader(self):
        sample = "sample"
        variant = retrieve_entry_from_test_vcf(2)
        interval = (1, 2)

        actual = Query.create_probe_header(sample, variant, interval)
        expected = (
            "GC00000001_155_SAMPLE=sample_POS=1_INTERVAL=(1,2)_SVTYPE=COMPLEX_"
            "MEAN_FWD_COVG=24_MEAN_REV_COVG=30_GT_CONF=262.757"
        )

        assert actual == expected


def test_isInvalidVcfEntry_withNoneGenotype_returnTrue():
    entry = retrieve_entry_from_test_vcf(0)
    sample = "sample"
    assert is_invalid_vcf_entry(entry, sample)


def test_isInvalidVcfEntry_withGenotype1_returnFalse():
    entry = retrieve_entry_from_test_vcf(1)
    sample = "sample"
    assert not is_invalid_vcf_entry(entry, sample)

def test_getGenotypeConfidence():
    entry = retrieve_entry_from_test_vcf(0)
    sample = "sample"
    assert get_genotype_confidence(entry, sample) == 262.757


def test_getSvtype():
    entry = retrieve_entry_from_test_vcf(0)

    actual = get_svtype(entry)
    expected = "COMPLEX"

    assert actual == expected


def test_getMeanCoverageForward():
    entry = retrieve_entry_from_test_vcf(2)
    sample = "sample"

    actual = get_mean_coverage_forward(entry, sample)
    expected = 24

    assert actual == expected


def test_getMeanCoverageReverse():
    entry = retrieve_entry_from_test_vcf(1)
    sample = "sample"

    actual = get_mean_coverage_reverse(entry, sample)
    expected = 7

    assert actual == expected


def test_getGenotype_genotypeNone_returnNone():
    entry = retrieve_entry_from_test_vcf(0)
    sample = "sample"
    assert get_genotype(entry, sample) is None


def test_getGenotype_genotype1_return1():
    entry = retrieve_entry_from_test_vcf(1)
    sample = "sample"
    assert get_genotype(entry, sample) == 1


def test_getVariantSequence_genotypeNone_returnRef():
    entry = retrieve_entry_from_test_vcf(0)
    sample = "sample"

    actual = get_variant_sequence(entry, sample)
    expected = "CTGCCCGTTGGC"

    assert actual == expected


def test_getVariantSequence_genotypeOne_returnFirstAlt():
    entry = retrieve_entry_from_test_vcf(1)
    sample = "sample"

    actual = get_variant_sequence(entry, sample)
    expected = "TTGGGGGAAGGCTCTGCACTGCCCGTTGGC"

    assert actual == expected


def test_getVariantSequence_genotypeZero_returnRef():
    entry = retrieve_entry_from_test_vcf(2)
    sample = "sample"

    actual = get_variant_sequence(entry, sample)
    expected = "CTGCCCGTTGGC"

    assert actual == expected


def test_getVariantLength_genotypeNone_returnRef():
    entry = retrieve_entry_from_test_vcf(0)
    sample = "sample"

    actual = get_variant_length(entry, sample)
    expected = 12

    assert actual == expected


def test_getVariantLength_genotypeOne_returnFirstAlt():
    entry = retrieve_entry_from_test_vcf(1)
    sample = "sample"

    actual = get_variant_length(entry, sample)
    expected = 30

    assert actual == expected


def test_getVariantLength_genotypeZero_returnRef():
    entry = retrieve_entry_from_test_vcf(2)
    sample = "sample"

    actual = get_variant_length(entry, sample)
    expected = 12

    assert actual == expected


def test_NoOverlappingIntervals_NoChange():
    intervals = [[2, 4], [6, 9], [11, 12]]

    result = merge_overlap_intervals(intervals)
    expected = [(2, 4), (6, 9), (11, 12)]
    assert result == expected


def test_TwoIntervalsEqualEndStart_NoChange():
    intervals = [[6, 9], [9, 12]]

    result = merge_overlap_intervals(intervals)
    expected = [(6, 9), (9, 12)]
    assert result == expected


def test_TwoIntervalsOverlap_Merge():
    intervals = [[6, 9], [8, 12]]

    result = merge_overlap_intervals(intervals)
    expected = [(6, 12)]
    assert result == expected


def test_ThreeIntervalsOverlap_Merge():
    intervals = [[6, 9], [8, 12], [11, 14]]

    result = merge_overlap_intervals(intervals)
    expected = [(6, 14)]
    assert result == expected


def test_ThreeIntervalsOverlapTwoEqualsEndStart_MergeOverlapDontMergeEquals():
    intervals = [[6, 9], [8, 12], [11, 14], [14, 16]]

    result = merge_overlap_intervals(intervals)
    expected = [(6, 14), (14, 16)]

    assert result == expected


def test_FindIndexInIntervals_emptyIntervalsReturnsNegativeOne():
    intervals = []
    query = 2

    actual = find_index_in_intervals(intervals, query)
    expected = -1

    assert actual == expected


def test_FindIndexInIntervals_queryNotInIntervalsReturnsNegativeOne():
    intervals = [(3, 7), (9, 20)]
    query = 2

    actual = find_index_in_intervals(intervals, query)
    expected = -1

    assert actual == expected


def test_FindIndexInIntervals_queryInFirstIntervalsReturnsZero():
    intervals = [(3, 7), (9, 20)]
    query = 5

    actual = find_index_in_intervals(intervals, query)
    expected = 0

    assert actual == expected


def test_FindIndexInIntervals_queryInSecondIntervalsReturnsOne():
    intervals = [(3, 7), (9, 20)]
    query = 14

    actual = find_index_in_intervals(intervals, query)
    expected = 1

    assert actual == expected


def test_FindIndexInIntervals_queryEqualsStartOfFirstIntervalReturnsZero():
    intervals = [(3, 7), (9, 20)]
    query = 3

    actual = find_index_in_intervals(intervals, query)
    expected = 0

    assert actual == expected


def test_FindIndexInIntervals_queryEqualsEndOfFirstIntervalReturnsZero():
    intervals = [(3, 7), (9, 20)]
    query = 7

    actual = find_index_in_intervals(intervals, query)
    expected = 0

    assert actual == expected


def test_FindIndexInIntervals_queryGreaterthanLastIntervalReturnsNegativeOne():
    intervals = [(3, 7), (9, 20)]
    query = 70

    actual = find_index_in_intervals(intervals, query)
    expected = -1

    assert actual == expected
