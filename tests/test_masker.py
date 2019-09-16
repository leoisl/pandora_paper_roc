from evaluate.masker import Masker, PrecisionMasker
from unittest.mock import patch, PropertyMock
from intervaltree import IntervalTree, Interval
from evaluate.probe import ProbeHeader, Probe
from evaluate.classification import (
    Classification,
    RecallClassification,
    PrecisionClassification,
)
from evaluate.aligned_pairs import AlignedPairs
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

    def test_filterRecords_noRecordsNoMaskReturnsEmpty(self):
        records = []
        masker = Masker()

        actual = masker.filter_records(records)
        expected = []

        assert actual == expected

    @patch.object(Masker, "record_overlaps_mask", return_value=False)
    def test_filterRecords_oneRecordDoesNotOverlapMaskReturnsRecord(self, *mock):
        records = [Classification()]
        masker = Masker()

        actual = masker.filter_records(records)
        expected = [Classification()]

        assert actual == expected

    @patch.object(Masker, "record_overlaps_mask", return_value=True)
    def test_filterRecords_oneRecordDoesOverlapMaskReturnsEmpty(self, *mock):
        records = [Classification()]
        masker = Masker()

        actual = masker.filter_records(records)
        expected = []

        assert actual == expected

    @patch.object(Masker, "record_overlaps_mask", side_effect=[True, False])
    def test_filterRecords_twoRecordsOneDoesOverlapMaskOneDoesntReturnsOneRecord(
        self, *mock
    ):
        record = Classification()
        record.query_probe = Probe(header=ProbeHeader(pos=4))
        records = [Classification(), record]
        masker = Masker()

        actual = masker.filter_records(records)
        expected = [record]

        assert actual == expected

    @patch.object(Masker, "record_overlaps_mask", return_value=False)
    def test_filterRecords_twoRecordsNoneOverlapMaskReturnsTwoRecords(self, *mock):
        record = Classification()
        record.query_probe = Probe(header=ProbeHeader(pos=4))
        records = [Classification(), record]
        masker = Masker()

        actual = masker.filter_records(records)
        expected = records

        assert actual == expected

    @patch.object(Masker, "record_overlaps_mask", return_value=True)
    def test_filterRecords_twoRecordsAllOverlapMaskReturnsNoRecords(self, *mock):
        record = Classification()
        record.query_probe = Probe(header=ProbeHeader(pos=4))
        records = [Classification(), record]
        masker = Masker()

        actual = masker.filter_records(records)
        expected = []

        assert actual == expected

    @patch.object(
        Masker,
        "get_interval_where_probe_aligns_to_truth",
        return_value=Interval(30, 40, "chrom1"),
    )
    def test_recordOverlapsMask_recordDoesNotOverlapReturnsFalse(self, *mock):
        masker = Masker(tree=IntervalTree([Interval(10, 20, "chrom1")]))
        record = Classification()

        actual = masker.record_overlaps_mask(record)
        expected = False

        assert actual == expected

    @patch.object(
        Masker,
        "get_interval_where_probe_aligns_to_truth",
        return_value=Interval(10, 20, "chrom2"),
    )
    def test_recordOverlapsMask_recordOverlapsButDifferentChromReturnsFalse(
        self, *mock
    ):
        masker = Masker(tree=IntervalTree([Interval(10, 20, "chrom1")]))
        record = Classification()

        actual = masker.record_overlaps_mask(record)
        expected = False

        assert actual == expected

    @patch.object(
        Masker,
        "get_interval_where_probe_aligns_to_truth",
        return_value=Interval(10, 20, "chrom1"),
    )
    def test_recordOverlapsMask_maskExactlyRecordReturnsTrue(self, *mock):
        masker = Masker(tree=IntervalTree([Interval(10, 20, "chrom1")]))
        record = Classification()

        actual = masker.record_overlaps_mask(record)
        expected = True

        assert actual == expected

    @patch.object(
        Masker,
        "get_interval_where_probe_aligns_to_truth",
        return_value=Interval(12, 18, "chrom1"),
    )
    def test_recordOverlapsMask_maskEnvelopsRecordReturnsTrue(self, *mock):
        masker = Masker(tree=IntervalTree([Interval(10, 20, "chrom1")]))
        record = Classification()

        actual = masker.record_overlaps_mask(record)
        expected = True

        assert actual == expected

    @patch.object(
        Masker,
        "get_interval_where_probe_aligns_to_truth",
        return_value=Interval(2, 23, "chrom1"),
    )
    def test_recordOverlapsMask_maskSpannedByRecordReturnsTrue(self, *mock):
        masker = Masker(tree=IntervalTree([Interval(10, 20, "chrom1")]))
        record = Classification()

        actual = masker.record_overlaps_mask(record)
        expected = True

        assert actual == expected

    @patch.object(
        Masker,
        "get_interval_where_probe_aligns_to_truth",
        return_value=Interval(2, 11, "chrom1"),
    )
    def test_recordOverlapsMask_recordOverlapsLeftEdgeOfMaskReturnsTrue(self, *mock):
        masker = Masker(tree=IntervalTree([Interval(10, 20, "chrom1")]))
        record = Classification()

        actual = masker.record_overlaps_mask(record)
        expected = True

        assert actual == expected

    @patch.object(
        Masker,
        "get_interval_where_probe_aligns_to_truth",
        return_value=Interval(2, 10, "chrom1"),
    )
    def test_recordOverlapsMask_recordMissesLeftEdgeOfMaskByOneReturnsFalse(
        self, *mock
    ):
        masker = Masker(tree=IntervalTree([Interval(10, 20, "chrom1")]))
        record = Classification()

        actual = masker.record_overlaps_mask(record)
        expected = False

        assert actual == expected

    @patch.object(
        Masker,
        "get_interval_where_probe_aligns_to_truth",
        return_value=Interval(19, 31, "chrom1"),
    )
    def test_recordOverlapsMask_recordOverlapsRightEdgeOfMaskReturnsTrue(self, *mock):
        masker = Masker(tree=IntervalTree([Interval(10, 20, "chrom1")]))
        record = Classification()

        actual = masker.record_overlaps_mask(record)
        expected = True

        assert actual == expected

    @patch.object(
        Masker,
        "get_interval_where_probe_aligns_to_truth",
        return_value=Interval(20, 30, "chrom1"),
    )
    def test_recordOverlapsMask_recordMissesRightEdgeOfMaskByOneReturnsFalse(
        self, *mock
    ):
        masker = Masker(tree=IntervalTree([Interval(10, 20, "chrom1")]))
        record = Classification()

        actual = masker.record_overlaps_mask(record)
        expected = False

        assert actual == expected

    @patch.object(Masker, "get_interval_where_probe_aligns_to_truth", return_value=None)
    def test_recordOverlapsMask_recordIsUnmappedReturnsFalse(self, *mock):
        masker = Masker(tree=IntervalTree([Interval(10, 20, "chrom1")]))
        record = Classification()

        actual = masker.record_overlaps_mask(record)
        expected = False

        assert actual == expected


class TestPrecisionMasker:
    @patch.object(
        Classification, "is_unmapped", return_value=True, new_callable=PropertyMock
    )
    def test_getIntervalWhereProbeAlignsToTruth_ProbeIsUnmappedReturnsNone(self, *mock):
        classification = Classification()

        actual = PrecisionMasker.get_interval_where_probe_aligns_to_truth(
            classification
        )
        expected = None

        assert actual == expected

    @patch.object(
        Classification, "is_unmapped", return_value=False, new_callable=PropertyMock
    )
    @patch.object(Probe, "chrom", return_value="chrom1", new_callable=PropertyMock)
    @patch.object(
        Classification,
        "get_probe_aligned_pairs",
        return_value=AlignedPairs(
            [(0, 34, "A"), (1, None, None), (2, 35, "A"), (None, 36, "A")]
        ),
    )
    def test_getIntervalWhereProbeAlignsToTruth_probeMapsReturnsInterval(self, *mock):
        classification = Classification()

        actual = PrecisionMasker.get_interval_where_probe_aligns_to_truth(
            classification
        )
        expected = Interval(34, 37, "chrom1")

        assert actual == expected

    @patch.object(
        Classification, "is_unmapped", return_value=False, new_callable=PropertyMock
    )
    @patch.object(Probe, "chrom", return_value="chrom1", new_callable=PropertyMock)
    @patch.object(
        Classification,
        "get_probe_aligned_pairs",
        return_value=AlignedPairs([(10, None, None), (11, None, None)]),
    )
    @patch.object(
        Classification,
        "get_aligned_pairs",
        return_value=AlignedPairs(
            [
                (8, 50, "A"),
                (9, None, None),
                (10, None, None),
                (11, None, None),
                (12, 51, "C"),
            ]
        ),
    )
    def test_getIntervalWhereProbeAlignsToTruth_probeIsInsertionReturnsIntervalAroundInsertion(
        self, *mock
    ):
        classification = Classification()

        actual = PrecisionMasker.get_interval_where_probe_aligns_to_truth(
            classification
        )
        expected = Interval(50, 52, "chrom1")

        assert actual == expected
