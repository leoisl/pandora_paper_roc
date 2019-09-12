from evaluate.masker import Masker
from intervaltree import IntervalTree, Interval
from io import StringIO


class TestMasker:
    def test_equality_twoMaskersTheSameReturnsTrue(self):
        ivs = [(1, 4), (6, 9)]
        m1 = Masker(tree=IntervalTree.from_tuples(ivs))
        m2 = Masker(tree=IntervalTree.from_tuples(ivs))

        assert m1 == m2

    def test_equality_twoMaskersNotTheSameReturnsFalse(self):
        ivs = [(1, 4), (6, 9)]
        m1 = Masker(tree=IntervalTree.from_tuples(ivs))
        m2 = Masker()

        assert m1 != m2

    def test_fromBed_emptyBedReturnsEmpty(self):
        bed = StringIO()

        actual = Masker.from_bed(bed)
        expected = Masker()

        assert actual == expected

    def test_fromBed_oneLineBedReturnsMaskerWithOneInterval(self):
        bed = StringIO("chrom\t3\t7")

        actual = Masker.from_bed(bed)
        expected = Masker(IntervalTree([Interval(3, 7, "chrom")]))

        assert actual == expected

    def test_fromBed_twoLinesBedReturnsMaskerWithTwoIntervals(self):
        bed = StringIO("chrom\t3\t7\nchrom\t8\t10")

        actual = Masker.from_bed(bed)
        expected = Masker(
            IntervalTree([Interval(3, 7, "chrom"), Interval(8, 10, "chrom")])
        )

        assert actual == expected

    def test_fromBed_twoLinesBedSameIntervalDiffChromosomeReturnsMaskerWithTwoIntervals(
        self
    ):
        bed = StringIO("chrom1\t3\t7\nchrom2\t3\t7")

        actual = Masker.from_bed(bed)
        expected = Masker(
            IntervalTree([Interval(3, 7, "chrom1"), Interval(3, 7, "chrom2")])
        )

        assert actual == expected

    def test_fromBed_twoLinesBedSameIntervalAndChromosomeReturnsMaskerWithOneInterval(
        self
    ):
        bed = StringIO("chrom\t3\t7\nchrom\t3\t7")

        actual = Masker.from_bed(bed)
        expected = Masker(IntervalTree([Interval(3, 7, "chrom")]))

        assert actual == expected
