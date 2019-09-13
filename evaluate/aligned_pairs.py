from enum import Enum
from typing import Optional, NamedTuple, Iterable
from collections import UserList
from intervaltree import Interval
from bisect import bisect_left, bisect_right


class AlignmentType(Enum):
    INSERTION = 0
    DELETION = 1
    MISMATCH = 2
    MATCH = 3


class AlignedPair(NamedTuple):
    query_pos: Optional[int] = None
    ref_pos: Optional[int] = None
    ref_base: Optional[str] = None

    def get_alignment_type(self) -> AlignmentType:
        if self.ref_base is None:
            return AlignmentType.DELETION
        elif self.query_pos is None:
            return AlignmentType.INSERTION
        elif self.ref_base.islower():
            return AlignmentType.MISMATCH
        else:
            return AlignmentType.MATCH


class AlignedPairs(UserList):
    def __init__(self, aligned_pairs: Optional[Iterable[tuple or AlignedPair]] = None):
        if aligned_pairs is None:
            aligned_pairs = []
        self.data = [AlignedPair(*aligned_pair) for aligned_pair in aligned_pairs]

    @staticmethod
    def transform_Nones_to_halfway_positions(
        array: Iterable[Optional[int]]
    ) -> Iterable[float]:
        if not array:
            return []
        try:
            first_int = next(item for item in array if item is not None)
        except StopIteration:
            raise ValueError("All values are None")

        previous_position = first_int - 1

        transformed_array = []
        for position in array:
            if position is None:
                transformed_array.append(previous_position + 0.5)
            else:
                transformed_array.append(position)
                previous_position = position

        return transformed_array

    def get_query_positions(
        self, transform_Nones_into_halfway_positions: bool = False
    ) -> Iterable[float]:
        query_positions = [pair.query_pos for pair in self.data]

        if transform_Nones_into_halfway_positions:
            return self.transform_Nones_to_halfway_positions(query_positions)
        else:
            return query_positions

    def get_alignment_types(self):
        return [aligned_pair.get_alignment_type() for aligned_pair in self]

    def get_pairs_in_query_interval(self, interval: Interval) -> "AlignedPairs":
        query_positions = self.get_query_positions(
            transform_Nones_into_halfway_positions=True
        )
        query_start = bisect_left(query_positions, interval.begin)
        query_stop = bisect_right(query_positions, interval.end - 1)
        return self[query_start:query_stop]
