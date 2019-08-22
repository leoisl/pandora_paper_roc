from pathlib import Path
from evaluate.utils import *


def test_stripExtensions_noExtensionsReturnsSameAsInput():
    path = Path("path/to/foo")

    actual = strip_extensions(path)
    expected = path

    assert actual == expected


def test_stripExtensions_oneExtensionReturnsInputWithoutExtension():
    path = Path("path/to/foo.fa")

    actual = strip_extensions(path)
    expected = Path("path/to/foo")

    assert actual == expected


def test_stripExtensions_twoExtensionsReturnsInputWithoutAnyExtensions():
    path = Path("path/to/foo.fa.gz")

    actual = strip_extensions(path)
    expected = Path("path/to/foo")

    assert actual == expected


def test_consecutive_notConsecutiveReturnsFalse():
    x = 2
    y = 6

    assert consecutive(x, y) is False


def test_consecutive_consecutiveIncreasingReturnsTrue():
    x = 2
    y = 3

    assert consecutive(x, y)


def test_consecutive_consecutiveDecreasingReturnsTrue():
    x = 2
    y = 1

    assert consecutive(x, y)


def test_consecutive_consecutiveSameReturnsTrue():
    x = 2
    y = 2

    assert consecutive(x, y)


def test_consecutive_offByOneReturnsFalse():
    x = 2
    y = 4

    assert consecutive(x, y) is False


def test_collapseRanges_emptyInEmptyOut():
    ranges = []

    actual = collapse_ranges(ranges)
    expected = []

    assert actual == expected


def test_collapseRanges_singleIndexRangeReturnsSingleCollapsedRange():
    ranges = [[0]]

    actual = collapse_ranges(ranges)
    expected = [(0, 1)]

    assert actual == expected


def test_collapseRanges_rangeWithManyEntriesReturnsSingleCollapsedRange():
    ranges = [[0, 1, 2, 3, 4]]

    actual = collapse_ranges(ranges)
    expected = [(0, 5)]

    assert actual == expected


def test_collapseRanges_rangeWithManyEntriesAndSingleRangeReturnsCollapsedRanges():
    ranges = [[0, 1, 2, 3, 4], [10]]

    actual = collapse_ranges(ranges)
    expected = [(0, 5), (10, 11)]

    assert actual == expected


def test_collapseRanges_singleRangeIsEmptyReturnsEmpty():
    ranges = [[]]

    actual = collapse_ranges(ranges)
    expected = []

    assert actual == expected


def test_collapseRanges_singleRangeIsEmptyAndNormalRangesReturnsJustNormalRanges():
    ranges = [[0, 1, 2, 3, 4], [10], []]

    actual = collapse_ranges(ranges)
    expected = [(0, 5), (10, 11)]

    assert actual == expected


def test_argRanges_emptyInEmptyOut():
    array = []

    actual = arg_ranges(array)
    expected = []

    assert actual == expected


def test_argRanges_singleElementReturnsSingleRange():
    array = [1]

    actual = arg_ranges(array)
    expected = [(0, 1)]

    assert actual == expected


def test_argRanges_singleRepeatedElementReturnsSingleRange():
    array = [1, 1, 1]

    actual = arg_ranges(array)
    expected = [(0, 3)]

    assert actual == expected


def test_argRanges_singleElementAndConsecutiveReturnsTwiRanges():
    array = [1, 5, 6, 6]

    actual = arg_ranges(array)
    expected = [(0, 1), (1, 4)]

    assert actual == expected


def test_argRanges_multipleConsecutiveRanges():
    array = [0, 2, 2, 4, 7, 8, 9, 15, 4, 4, 7, 8, 11, 16, 4, 3, 2]

    actual = arg_ranges(array)
    expected = [
        (0, 1),
        (1, 3),
        (3, 4),
        (4, 7),
        (7, 8),
        (8, 10),
        (10, 12),
        (12, 13),
        (13, 14),
        (14, 17),
    ]

    assert actual == expected
