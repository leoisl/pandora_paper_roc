from evaluate.filter import Filter
from typing import Iterable
from collections import Counter

class UniqueSamRecordsFilter(Filter):
    def __init__(self, records: Iterable):
        all_query_names = [record.query_name for record in records]
        query_names_counter = Counter(all_query_names)

        self._all_unique_query_names = []
        for query_name, occurence in query_names_counter.items():
            if occurence == 1:
                self._all_unique_query_names.append(query_name)
        self._all_unique_query_names = sorted(self._all_unique_query_names)

    @property
    def all_unique_query_names(self):
        return self._all_unique_query_names

    def record_should_be_filtered_out(self, record) -> bool:
        query_name = record.query_name
        return query_name not in self.all_unique_query_names