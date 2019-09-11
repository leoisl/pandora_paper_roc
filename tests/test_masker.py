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

# todo: add tests for different chromsomes
