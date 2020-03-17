from evaluate.query import *
from tests.common import (
    TEST_CASES,
    TEST_QUERY_VCF,
    TEST_QUERY_REF,
    retrieve_entry_from_test_vcf,
    retrieve_entry_from_test_query_vcf,
)
from evaluate.vcf_file import VCFFile, VCFFactory
import pysam


class TestQuery:
    def test_calculateProbeBoundariesForEntry_variantShorterThanMinLen_returnProbeOfMinLen(
        self
    ):
        flank_width = 10
        query = Query(None, None, None, flank_width=flank_width)
        variant = retrieve_entry_from_test_query_vcf(1)

        expected = (9, 32)
        actual = query.calculate_probe_boundaries_for_entry(variant)

        assert actual == expected

    def test_calculateProbeBoundariesForEntry_variantAtStartOfGene_returnZeroLenLeftProbe(
        self
    ):
        flank_width = 10
        query = Query(None, None, None, flank_width=flank_width)
        variant = retrieve_entry_from_test_query_vcf(0)

        expected = (0, 13)
        actual = query.calculate_probe_boundaries_for_entry(variant)

        assert actual == expected

    def test_calculateProbeBoundariesForEntry_variantLongerThanMinLen_returnZeroLenProbes(
        self
    ):
        flank_width = 2
        query = Query(None, None, None, flank_width=flank_width)
        variant = retrieve_entry_from_test_query_vcf(1)

        expected = (17, 24)
        actual = query.calculate_probe_boundaries_for_entry(variant)

        assert actual == expected

    def test_createProbesForGeneVariants_emptyVariants_returnEmptyProbes(self):
        samples = ["sample"]
        pysam_variant_file = pysam.VariantFile(TEST_QUERY_VCF, "r")
        vcf_file = VCFFile(pysam_variant_file, VCFFactory.create_Pandora_VCF_from_VariantRecord_and_Sample)
        pysam_variant_file.close()

        query = Query(vcf_file, TEST_QUERY_REF, samples)

        expected = {sample: "" for sample in samples}
        actual = query._create_probes_for_gene_variants(pysam.FastxRecord(), [])

        assert actual == expected

    def test_makeProbes_emptyVariantsReturnsEmptyProbes(self):
        samples = ["sample"]
        pysam_variant_file = pysam.VariantFile(TEST_CASES / "empty.vcf", "r")
        vcf_file = VCFFile(pysam_variant_file, VCFFactory.create_Pandora_VCF_from_VariantRecord_and_Sample)
        pysam_variant_file.close()
        genes = TEST_QUERY_REF

        query = Query(vcf_file, genes, samples)

        actual = query.make_probes()
        expected = {s: "" for s in samples}

        assert actual == expected

    def test_makeProbes_emptyGenesReturnsEmptyProbes(self):
        samples = ["sample"]
        pysam_variant_file = pysam.VariantFile(TEST_CASES / "empty.vcf", "r")
        vcf_file = VCFFile(pysam_variant_file, VCFFactory.create_Pandora_VCF_from_VariantRecord_and_Sample)
        pysam_variant_file.close()
        genes = TEST_CASES / "empty.fa"

        query = Query(vcf_file, genes, samples)

        actual = query.make_probes()
        expected = {s: "" for s in samples}

        assert actual == expected

    def test_makeProbes_oneGeneOneVcfRecordNotInGeneReturnsEmptyProbes(self):
        samples = ["sample"]
        pysam_variant_file = pysam.VariantFile(TEST_CASES / "empty.vcf", "r")
        vcf_file = VCFFile(pysam_variant_file, VCFFactory.create_Pandora_VCF_from_VariantRecord_and_Sample)
        pysam_variant_file.close()
        genes = TEST_QUERY_REF

        query = Query(vcf_file, genes, samples)

        actual = query.make_probes()
        expected = {s: "" for s in samples}

        assert actual == expected

    def test_makeProbes_oneGeneOneVcfRecordInGeneReturnsOneProbe(self):
        samples = ["sample"]
        pysam_variant_file = pysam.VariantFile(TEST_CASES / "make_probes_1.vcf", "r")
        vcf_file = VCFFile(pysam_variant_file, VCFFactory.create_Pandora_VCF_from_VariantRecord_and_Sample)
        pysam_variant_file.close()
        genes = TEST_CASES / "make_probes_1.fa"
        flank_width = 3

        query = Query(vcf_file, genes, samples, flank_width)

        actual = query.make_probes()
        expected = {
            "sample": str(
                Probe(
                    ProbeHeader(
                        chrom="gene1",
                        sample="sample",
                        pos=4,
                        interval=ProbeInterval(3, 4),
                        svtype="COMPLEX",
                        gt_conf=262.757,
                    ),
                    full_sequence="xxxFxxx",
                )
            )
            + "\n"
        }

        assert actual == expected

    def test_makeProbes_oneGeneTwoNonCloseVcfRecordsInGeneReturnsTwoProbes(self):
        samples = ["sample"]
        pysam_variant_file = pysam.VariantFile(TEST_CASES / "make_probes_3.vcf", "r")
        vcf_file = VCFFile(pysam_variant_file, VCFFactory.create_Pandora_VCF_from_VariantRecord_and_Sample)
        pysam_variant_file.close()
        genes = TEST_CASES / "make_probes_2.fa"
        flank_width = 5

        query = Query(vcf_file, genes, samples, flank_width)

        actual = query.make_probes()
        expected = {
            "sample": str(
                Probe(
                    ProbeHeader(
                        chrom="gene1",
                        sample="sample",
                        pos=4,
                        interval=ProbeInterval(3, 6),
                        svtype="COMPLEX",
                        gt_conf=262.757,
                    ),
                    full_sequence="xxxFOOxxxxx",
                )
            )
            + "\n"
            + str(
                Probe(
                    ProbeHeader(
                        chrom="gene1",
                        sample="sample",
                        pos=8,
                        interval=ProbeInterval(5, 8),
                        svtype="COMPLEX",
                        gt_conf=262.757,
                    ),
                    full_sequence="xxxxxFOOxxxxx",
                )
            )
            + "\n"
        }

        assert actual == expected

    def test_makeProbes_twoGenesTwoNonCloseVcfRecordsInOneGeneReturnsTwoProbes(self):
        samples = ["sample"]
        pysam_variant_file = pysam.VariantFile(TEST_CASES / "make_probes_3.vcf", "r")
        vcf_file = VCFFile(pysam_variant_file, VCFFactory.create_Pandora_VCF_from_VariantRecord_and_Sample)
        pysam_variant_file.close()
        genes = TEST_CASES / "make_probes_3.fa"
        flank_width = 5

        query = Query(vcf_file, genes, samples, flank_width)

        actual = query.make_probes()
        expected = {
            "sample": str(
                Probe(
                    ProbeHeader(
                        chrom="gene1",
                        sample="sample",
                        pos=4,
                        interval=ProbeInterval(3, 6),
                        svtype="COMPLEX",
                        gt_conf=262.757,
                    ),
                    full_sequence="xxxFOOxxxxx",
                )
            )
            + "\n"
            + str(
                Probe(
                    ProbeHeader(
                        chrom="gene1",
                        sample="sample",
                        pos=8,
                        interval=ProbeInterval(5, 8),
                        svtype="COMPLEX",
                        gt_conf=262.757,
                    ),
                    full_sequence="xxxxxFOOxxx",
                )
            )
            + "\n"
        }

        assert actual == expected

    def test_makeProbes_twoGenesTwoVcfRecordsOneInEachGeneReturnsTwoProbes(self):
        samples = ["sample"]
        pysam_variant_file = pysam.VariantFile(TEST_CASES / "make_probes_4.vcf", "r")
        vcf_file = VCFFile(pysam_variant_file, VCFFactory.create_Pandora_VCF_from_VariantRecord_and_Sample)
        pysam_variant_file.close()
        genes = TEST_CASES / "make_probes_3.fa"
        flank_width = 5

        query = Query(vcf_file, genes, samples, flank_width)

        actual = query.make_probes()
        expected = {
            "sample": str(
                Probe(
                    ProbeHeader(
                        chrom="gene1",
                        sample="sample",
                        pos=4,
                        interval=ProbeInterval(3, 6),
                        svtype="COMPLEX",
                        gt_conf=262.757,
                    ),
                    full_sequence="xxxFOOxxxxx",
                )
            )
            + "\n"
            + str(
                Probe(
                    ProbeHeader(
                        chrom="gene2",
                        sample="sample",
                        pos=2,
                        interval=ProbeInterval(1, 4),
                        svtype="COMPLEX",
                        gt_conf=262.757,
                    ),
                    full_sequence="xFOOx",
                )
            )
            + "\n"
        }

        assert actual == expected

    def test_makeProbes_oneGeneTwoVcfRecordsInTheSameIntervalWithDifferentGTConfReturnsProbeWithHighestGTConf(
        self
    ):
        samples = ["sample"]
        pysam_variant_file = pysam.VariantFile(TEST_CASES / "make_probes_6.vcf", "r")
        vcf_file = VCFFile(pysam_variant_file, VCFFactory.create_Pandora_VCF_from_VariantRecord_and_Sample)
        pysam_variant_file.close()
        genes = TEST_CASES / "make_probes_6.fa"
        flank_width = 5

        query = Query(vcf_file, genes, samples, flank_width)

        actual = query.make_probes()
        expected = {
            "sample": str(
                Probe(
                    ProbeHeader(
                        chrom="gene1",
                        sample="sample",
                        pos=4,
                        interval=ProbeInterval(3, 4),
                        svtype="COMPLEX",
                        gt_conf=20.0,
                    ),
                    full_sequence="xxxFxxxxx",
                )
            )
            + "\n"
        }

        assert actual == expected

    def test_createProbeHeader(self):
        sample = "sample"
        flank_width = 3
        query = Query(None, None, None, flank_width=flank_width)
        variant = retrieve_entry_from_test_vcf(1)
        vcf = VCFFactory.create_Pandora_VCF_from_VariantRecord_and_Sample(variant, sample)
        interval = query.calculate_probe_boundaries_for_entry(vcf)

        actual = Query._create_probe_header(sample, vcf, interval)
        expected = ProbeHeader(
            chrom="GC00000001_155",
            sample="sample",
            pos=1,
            interval=ProbeInterval(0, 30),
            svtype="COMPLEX",
            gt_conf=262.757,
        )

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
