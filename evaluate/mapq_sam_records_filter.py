from evaluate.filter import Filter
from typing import List, Optional, Dict
from collections import defaultdict
from pysam import AlignedSegment

class MAPQSamRecordsFilter(Filter):
    """
    Filters records based on MapQ, which tells us how confident we can be that the read comes from the reported position.

    If there are multiple alignments, keeps only the one with best MAPQ, based on the following constraints:
        -If the difference between the best and the second best MAPQ is too small (i.e. < 10), then discard all alignments to that probe (i.e. probe is ignored in the evaluation)
            -The reason is that we just have ambiguous data to say where it really came from, so the best is to avoid the evaluation of this probe
        -Otherwise, keeps only the one with best MAPQ

    Reasoning from https://www.biostars.org/p/179457/#299885:
    Alignment score is a metric that tells you how similar the read is to the reference. AS increases with the number of
    matches and decreases with the number of mismatches and gaps (rewards and penalties for matches and mismatches depend
    on the scoring matrix you use). MAPQ is a metric that tells you how confident you can be that the read comes from the
    reported position.
    You can have high AS and low MAPQ if the read aligns perfectly at multiple positions, and you can have low AS and
    high MAPQ if the read aligns with mismatches but still the reported position is still much more probable than any other.
    """
    def __init__(self, records: List[AlignedSegment], mapping_quality_threshold: int = 10):
        self._mapping_quality_threshold = mapping_quality_threshold

        query_name_to_best_record = self.get_query_name_to_best_record(records)
        self._records_to_keep = set(query_name_to_best_record.values())



    def get_record_with_highest_mapping_quality (self, records: List[AlignedSegment]) -> Optional[AlignedSegment]:
        if len(records) == 0:
            return None
        all_mapqs = [record.mapping_quality for record in records]
        highest_mapping_quality = max(all_mapqs)
        for record in records:
            if record.mapping_quality == highest_mapping_quality:
                return record
        assert False, "Should never reach here"

    def get_query_name_to_best_record(self, records: List) -> Dict[str, AlignedSegment]:
        query_name_to_records = defaultdict(list)
        for record in records:
            query_name_to_records[record.query_name].append(record)

        query_name_to_best_record = {}
        for query_name, records in query_name_to_records.items():
            records_copy = list(records)
            record_with_highest_mapping_quality = self.get_record_with_highest_mapping_quality(records_copy)
            records_copy.remove(record_with_highest_mapping_quality)
            record_with_second_highest_mapping_quality = self.get_record_with_highest_mapping_quality(
                records_copy)

            mapping_quality_gap_is_too_small = \
                record_with_highest_mapping_quality.mapping_quality - record_with_second_highest_mapping_quality.mapping_quality <= \
                self._mapping_quality_threshold

            if not mapping_quality_gap_is_too_small:
                query_name_to_best_record[query_name] = record_with_highest_mapping_quality

        return query_name_to_best_record

    @property
    def records_to_keep(self):
        return self._records_to_keep

    def record_should_be_filtered_out(self, record) -> bool:
        return record not in self.records_to_keep