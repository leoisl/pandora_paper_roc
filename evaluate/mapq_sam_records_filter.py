from evaluate.filter import Filter
from evaluate.classification import Classification
from typing import List, Optional, Dict, Tuple
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

        query_name_to_best_record = self._get_query_name_to_best_record(records)
        self._records_to_keep = set(query_name_to_best_record.values())


    def _get_only_records_that_cover_the_allele(self, records: List[AlignedSegment]) -> List[AlignedSegment]:
        classifications = [Classification(record) for record in records]
        return self._get_only_records_that_cover_the_allele_core(classifications)

    def _get_only_records_that_cover_the_allele_core(self, classifications: List[Classification]) -> List[AlignedSegment]:
        records_that_covers_the_allele = [classification.record for classification in classifications if
                                          classification._whole_query_probe_maps()]
        return records_that_covers_the_allele


    def _get_record_with_highest_mapping_quality (self, records: List[AlignedSegment]) -> Optional[AlignedSegment]:
        if len(records) == 0:
            return None
        all_mapqs = [record.mapping_quality for record in records]
        highest_mapping_quality = max(all_mapqs)
        for record in records:
            if record.mapping_quality == highest_mapping_quality:
                return record
        assert False, "Should never reach here"


    def _get_first_and_second_records_with_highest_mapping_qualities_that_covers_the_allele(self, records: List[AlignedSegment]) -> Tuple[Optional[AlignedSegment], Optional[AlignedSegment]]:
        records_that_covers_the_allele = self._get_only_records_that_cover_the_allele(records)

        if len(records_that_covers_the_allele) == 0:
            return None, None

        record_with_highest_mapping_quality_that_covers_the_allele = self._get_record_with_highest_mapping_quality(
            records_that_covers_the_allele)
        records_that_covers_the_allele.remove(record_with_highest_mapping_quality_that_covers_the_allele)
        record_with_second_highest_mapping_quality_that_covers_the_allele = self._get_record_with_highest_mapping_quality(
            records_that_covers_the_allele)
        return record_with_highest_mapping_quality_that_covers_the_allele, record_with_second_highest_mapping_quality_that_covers_the_allele

    def _get_query_name_to_best_record(self, records: List) -> Dict[str, AlignedSegment]:
        query_name_to_records = defaultdict(list)
        for record in records:
            query_name_to_records[record.query_name].append(record)

        query_name_to_best_record = {}
        for query_name, records in query_name_to_records.items():
            there_is_only_one_record = len(records) == 1

            if there_is_only_one_record:
                # Mapping found only one record, we don't know if it is a FP or TP, but should be evaluated
                query_name_to_best_record[query_name] = records[0]
            else:
                # More than 1 record for this query, choose the best that covers the allele
                record_with_highest_mapping_quality_that_covers_the_allele, record_with_second_highest_mapping_quality_that_covers_the_allele = \
                self._get_first_and_second_records_with_highest_mapping_qualities_that_covers_the_allele(records)

                there_is_at_least_one_record_that_covers_the_allele = record_with_highest_mapping_quality_that_covers_the_allele is not None

                if there_is_at_least_one_record_that_covers_the_allele:
                    there_is_no_second_record_that_covers_the_allele = record_with_second_highest_mapping_quality_that_covers_the_allele is None
                    if there_is_no_second_record_that_covers_the_allele:
                        # There is only one record covering the allele, thus let's select it
                        query_name_to_best_record[query_name] = record_with_highest_mapping_quality_that_covers_the_allele
                    else:
                        # Several records covering the allele, check mapping quality gap to choose a single one
                        mapping_quality_gap_is_large_enough = \
                            record_with_highest_mapping_quality_that_covers_the_allele.mapping_quality -\
                            record_with_second_highest_mapping_quality_that_covers_the_allele.mapping_quality > \
                            self._mapping_quality_threshold

                        if mapping_quality_gap_is_large_enough:
                            # We have a clear better mapping that covers the allele, select it
                            query_name_to_best_record[query_name] = record_with_highest_mapping_quality_that_covers_the_allele
                        else:
                            # Here we have 2+ good mappings that cover entirely the allele
                            # It is blurry to know what to do here.
                            # If we choose the one with best precision (i.e. the one where the allele bases matches the ref the most),
                            # we could be favouring the caller if that was not the call.
                            # If we choose the one with worst precision, we could be disfavouring the caller.
                            # This is too ambiguous, we don't have enough information to make a decision, so we just
                            # remove this case from the evaluation
                            pass
                else:
                    # No record covers the allele, so it is a FP.
                    # Chooses a random record to be evaluated as FP, as this is clearly a mistake from the variant caller.
                    query_name_to_best_record[query_name] = records[0]

        return query_name_to_best_record

    @property
    def records_to_keep(self):
        return self._records_to_keep

    def record_should_be_filtered_out(self, record) -> bool:
        return record not in self.records_to_keep