from typing import Iterable, List


class Filter:
    def filter_records(self, records: Iterable) -> List:
        return [
            record
            for record in records
            if not self.record_should_be_filtered_out(record)
        ]

    def record_should_be_filtered_out(self, record) -> bool:
        raise NotImplementedError()
