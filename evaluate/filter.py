from typing import Iterable
import pysam


class Filter:
    def filter_records(
        self, records: Iterable[pysam.AlignedSegment]
    ) -> Iterable[pysam.AlignedSegment]:
        return [
            record
            for record in records
            if not self.record_should_be_filtered_out(record)
        ]

    def record_should_be_filtered_out(self, record: pysam.AlignedSegment) -> bool:
        raise NotImplementedError()
