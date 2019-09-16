import pytest
from more_itertools import side_effect
from intervaltree import Interval

from evaluate.aligned_pairs import AlignedPair, AlignmentType, AlignedPairs
from unittest.mock import patch, PropertyMock


class TestAlignedPair:
    def test_getAlignmentType_pairIsInsertion(self):
        aligned_pair = AlignedPair(query_pos=None, ref_pos=10, ref_base="A")

        actual = aligned_pair.get_alignment_type()
        expected = AlignmentType.INSERTION

        assert actual == expected

    def test_getAlignmentType_pairIsDeletion(self):
        aligned_pair = AlignedPair(query_pos=10, ref_pos=None, ref_base=None)

        actual = aligned_pair.get_alignment_type()
        expected = AlignmentType.DELETION

        assert actual == expected

    def test_getAlignmentType_pairIsMatch(self):
        aligned_pair = AlignedPair(query_pos=10, ref_pos=15, ref_base="A")

        actual = aligned_pair.get_alignment_type()
        expected = AlignmentType.MATCH

        assert actual == expected

    def test_getAlignmentType_pairIsMismatch(self):
        aligned_pair = AlignedPair(query_pos=10, ref_pos=15, ref_base="a")

        actual = aligned_pair.get_alignment_type()
        expected = AlignmentType.MISMATCH

        assert actual == expected


class TestAlignedPairs:
    def test_initAlignedPairsNoArguments(self):
        actual = AlignedPairs()
        expected = []

        assert actual == expected

    def test_initAlignedPairsEmptyList(self):
        actual = AlignedPairs([])
        expected = []

        assert actual == expected

    def test_initAlignedPairsListWith4AlignedPairs(self):
        aligned_pairs_as_tuples = [
            (0, 34, "A"),
            (1, None, None),
            (2, 35, "A"),
            (None, 36, "A"),
        ]

        actual = AlignedPairs(aligned_pairs_as_tuples)
        expected = [
            AlignedPair(*aligned_pair_as_tuple)
            for aligned_pair_as_tuple in aligned_pairs_as_tuples
        ]

        assert actual == expected

    @patch.object(
        AlignedPair,
        "get_alignment_type",
        side_effect=[
            AlignmentType.MATCH,
            AlignmentType.MISMATCH,
            AlignmentType.DELETION,
            AlignmentType.INSERTION,
        ],
    )
    def test_getAlignmentTypes(self, alignedPairMocked):
        aligned_pairs = AlignedPairs([AlignedPair()] * 4)

        actual = aligned_pairs.get_alignment_types()
        expected = [
            AlignmentType.MATCH,
            AlignmentType.MISMATCH,
            AlignmentType.DELETION,
            AlignmentType.INSERTION,
        ]

        assert actual == expected

    def test_getQueryPositions_emptyAlignedPairsReturnsEmptyList(self, *mocks):
        aligned_pairs = AlignedPairs()

        actual = aligned_pairs.get_query_positions()
        expected = []

        assert actual == expected

    def test_getQueryPositions_allAlignedPairsAreNoneNoTransformReturnsNone(self):
        aligned_pairs = AlignedPairs([AlignedPair(None, 10, "A")])

        actual = aligned_pairs.get_query_positions()
        expected = [None]

        assert actual == expected

    @patch.object(
        AlignedPairs, "transform_Nones_to_halfway_positions", side_effect=ValueError
    )
    def test_getQueryPositions_allAlignedPairsAreNoneAndTransformRaisesException(
        self, *mocks
    ):
        aligned_pairs = AlignedPairs([AlignedPair(None, 10, "A")])
        with pytest.raises(ValueError):
            aligned_pairs.get_query_positions(
                transform_Nones_into_halfway_positions=True
            )

    def test_getQueryPositions_twoPairsOneNoneNoTransformReturnsTwoPositionsOneNone(
        self
    ):
        aligned_pairs = AlignedPairs(
            [AlignedPair(None, 10, "A"), AlignedPair(5, 11, "C")]
        )

        actual = aligned_pairs.get_query_positions()
        expected = [None, 5]

        assert actual == expected

    @patch.object(
        AlignedPairs, "transform_Nones_to_halfway_positions", return_value=[4.5, 5.0]
    )
    def test_getQueryPositions_twoPairsOneNoneTransformReturnsTwoPositionsOneHalfway(
        self, *mocks
    ):
        aligned_pairs = AlignedPairs(
            [AlignedPair(None, 10, "A"), AlignedPair(5, 11, "C")]
        )

        actual = aligned_pairs.get_query_positions(
            transform_Nones_into_halfway_positions=True
        )
        expected = [4.5, 5.0]

        assert actual == expected

    def test_getRefPositions_emptyAlignedPairsReturnsEmptyList(self, *mocks):
        aligned_pairs = AlignedPairs()

        actual = aligned_pairs.get_ref_positions()
        expected = []

        assert actual == expected

    def test_getRefPositions_allAlignedPairsAreNoneNoTransformReturnsNone(self):
        aligned_pairs = AlignedPairs([AlignedPair(1, None, None)])

        actual = aligned_pairs.get_ref_positions()
        expected = [None]

        assert actual == expected

    @patch.object(
        AlignedPairs, "transform_Nones_to_halfway_positions", side_effect=ValueError
    )
    def test_getRefPositions_allAlignedPairsAreNoneAndTransformRaisesException(
        self, *mocks
    ):
        aligned_pairs = AlignedPairs([AlignedPair(10, None, None)])
        with pytest.raises(ValueError):
            aligned_pairs.get_ref_positions(transform_Nones_into_halfway_positions=True)

    def test_getRefPositions_twoPairsOneNoneNoTransformReturnsTwoPositionsOneNone(self):
        aligned_pairs = AlignedPairs(
            [AlignedPair(4, None, None), AlignedPair(5, 11, "C")]
        )

        actual = aligned_pairs.get_ref_positions()
        expected = [None, 11]

        assert actual == expected

    @patch.object(
        AlignedPairs, "transform_Nones_to_halfway_positions", return_value=[4.5, 5.0]
    )
    def test_getRefPositions_twoPairsOneNoneTransformReturnsTwoPositionsOneHalfway(
        self, *mocks
    ):
        aligned_pairs = AlignedPairs(
            [AlignedPair(1, None, None), AlignedPair(2, 5, "C")]
        )

        actual = aligned_pairs.get_ref_positions(
            transform_Nones_into_halfway_positions=True
        )
        expected = [4.5, 5.0]

        assert actual == expected

    def test_transformNonesToHalfwayPositions_emptyReturnsEmpty(self):
        actual = AlignedPairs.transform_Nones_to_halfway_positions([])
        expected = []

        assert actual == expected

    def test_transformNonesToHalfwayPositions_noNonesReturnsInput(self):
        array = [1, 2]
        actual = AlignedPairs.transform_Nones_to_halfway_positions(array)
        expected = array

        assert actual == expected

    def test_transformNonesToHalfwayPositions_severalNonesReturnsTransformedArray(self):
        array = [None, 5, None, None, 6, None, None]
        actual = AlignedPairs.transform_Nones_to_halfway_positions(array)
        expected = [4.5, 5, 5.5, 5.5, 6, 6.5, 6.5]

        assert actual == expected

    def test_transformNonesToHalfwayPositions_allNonesButOneReturnsTransformedArray(
        self
    ):
        array = [None, None, 0, None, None]
        actual = AlignedPairs.transform_Nones_to_halfway_positions(array)
        expected = [-0.5, -0.5, 0, 0.5, 0.5]

        assert actual == expected

    def test_transformNonesToHalfwayPositions_allNonesRaisesException(self):
        with pytest.raises(ValueError):
            AlignedPairs.transform_Nones_to_halfway_positions([None, None])

    def test_getPairsInQueryInterval_nullIntervalReturnsEmpty(self):
        interval = Interval(2, 2)

        aligned_pairs = AlignedPairs(
            [AlignedPair(1, 30, "A"), AlignedPair(2, 31, "C"), AlignedPair(3, 32, "T")]
        )
        actual = aligned_pairs.get_pairs_in_query_interval(interval)
        expected = AlignedPairs()

        assert actual == expected

    def test_getPairsInQueryInterval_intervalNotInPairsReturnsEmpty(self):
        interval = Interval(5, 10)

        aligned_pairs = AlignedPairs(
            [AlignedPair(1, 30, "A"), AlignedPair(2, 31, "C"), AlignedPair(3, 32, "T")]
        )
        actual = aligned_pairs.get_pairs_in_query_interval(interval)
        expected = AlignedPairs()

        assert actual == expected

    def test_getPairsInQueryInterval_intervalOverlapsLeftOfPairs(self):
        interval = Interval(5, 12)

        aligned_pairs = AlignedPairs(
            [
                AlignedPair(None, 30, "A"),
                AlignedPair(11, 31, "C"),
                AlignedPair(12, 32, "T"),
            ]
        )
        actual = aligned_pairs.get_pairs_in_query_interval(interval)
        expected = aligned_pairs[:2]

        assert actual == expected

    def test_getPairsInQueryInterval_intervalOverlapsRightOfPairs(self):
        interval = Interval(4, 12)

        aligned_pairs = AlignedPairs(
            [
                AlignedPair(3, 30, "A"),
                AlignedPair(4, 31, "C"),
                AlignedPair(None, 32, "T"),
            ]
        )
        actual = aligned_pairs.get_pairs_in_query_interval(interval)
        expected = aligned_pairs[1:]

        assert actual == expected

    def test_getPairsInQueryInterval_intervalSpansPairs(self):
        interval = Interval(1, 10)

        aligned_pairs = AlignedPairs(
            [
                AlignedPair(4, 30, "A"),
                AlignedPair(None, 31, "A"),
                AlignedPair(5, 32, "A"),
                AlignedPair(None, 33, "A"),
                AlignedPair(6, 34, "A"),
                AlignedPair(None, 35, "A"),
            ]
        )
        actual = aligned_pairs.get_pairs_in_query_interval(interval)
        expected = aligned_pairs

        assert actual == expected

    def test_getPairsInQueryInterval_intervalEnvelopedInPairs(self):
        interval = Interval(5, 7)

        aligned_pairs = AlignedPairs(
            [
                AlignedPair(4, 30, "A"),
                AlignedPair(None, 31, "A"),
                AlignedPair(5, 32, "A"),
                AlignedPair(None, 33, "A"),
                AlignedPair(6, 34, "A"),
                AlignedPair(None, 35, "A"),
            ]
        )
        actual = aligned_pairs.get_pairs_in_query_interval(interval)
        expected = aligned_pairs[2:5]

        assert actual == expected

    def test_getIndexOfQueryInterval_nullIntervalReturnsEmpty(self):
        interval = Interval(2, 2)

        aligned_pairs = AlignedPairs(
            [AlignedPair(1, 30, "A"), AlignedPair(2, 31, "C"), AlignedPair(3, 32, "T")]
        )
        actual = aligned_pairs.get_index_of_query_interval(interval)
        expected = (1, 1)

        assert actual == expected

    def test_getIndexOfQueryInterval_intervalNotInPairsReturnsEmpty(self):
        interval = Interval(5, 10)

        aligned_pairs = AlignedPairs(
            [AlignedPair(1, 30, "A"), AlignedPair(2, 31, "C"), AlignedPair(3, 32, "T")]
        )
        actual = aligned_pairs.get_index_of_query_interval(interval)
        expected = (3, 3)

        assert actual == expected

    def test_getIndexOfQueryInterval_intervalOverlapsLeftOfPairs(self):
        interval = Interval(5, 12)

        aligned_pairs = AlignedPairs(
            [
                AlignedPair(None, 30, "A"),
                AlignedPair(11, 31, "C"),
                AlignedPair(12, 32, "T"),
            ]
        )
        actual = aligned_pairs.get_index_of_query_interval(interval)
        expected = (0, 2)

        assert actual == expected

    def test_getIndexOfQueryInterval_intervalOverlapsRightOfPairs(self):
        interval = Interval(4, 12)

        aligned_pairs = AlignedPairs(
            [
                AlignedPair(3, 30, "A"),
                AlignedPair(4, 31, "C"),
                AlignedPair(None, 32, "T"),
            ]
        )
        actual = aligned_pairs.get_index_of_query_interval(interval)
        expected = (1, 3)

        assert actual == expected

    def test_getIndexOfQueryInterval_intervalSpansPairs(self):
        interval = Interval(1, 10)

        aligned_pairs = AlignedPairs(
            [
                AlignedPair(4, 30, "A"),
                AlignedPair(None, 31, "A"),
                AlignedPair(5, 32, "A"),
                AlignedPair(None, 33, "A"),
                AlignedPair(6, 34, "A"),
                AlignedPair(None, 35, "A"),
            ]
        )
        actual = aligned_pairs.get_index_of_query_interval(interval)
        expected = (0, 6)

        assert actual == expected

    def test_getIndexOfQueryInterval_intervalEnvelopedInPairs(self):
        interval = Interval(5, 7)

        aligned_pairs = AlignedPairs(
            [
                AlignedPair(4, 30, "A"),
                AlignedPair(None, 31, "A"),
                AlignedPair(5, 32, "A"),
                AlignedPair(None, 33, "A"),
                AlignedPair(6, 34, "A"),
                AlignedPair(None, 35, "A"),
            ]
        )
        actual = aligned_pairs.get_index_of_query_interval(interval)
        expected = (2, 5)

        assert actual == expected
