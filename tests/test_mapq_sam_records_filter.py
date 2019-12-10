from unittest.mock import patch, Mock, PropertyMock
from evaluate.mapq_sam_records_filter import MAPQSamRecordsFilter


mapping_quality_60_mock = Mock(query_name="query_name_1", mapping_quality = 60)
mapping_quality_51_mock = Mock(query_name="query_name_1", mapping_quality = 51)
mapping_quality_50_mock = Mock(query_name="query_name_1", mapping_quality = 50)
mapping_quality_49_mock = Mock(query_name="query_name_1", mapping_quality = 49)
mapping_quality_40_mock = Mock(query_name="query_name_1", mapping_quality = 49)

class TestMAPQSamRecordsFilter:
    @patch.object(MAPQSamRecordsFilter, MAPQSamRecordsFilter.get_query_name_to_best_record.__name__, return_value={})
    def test___constructor___no_sam_records___records_to_keep_is_empty(self, get_query_name_to_best_record_mock):
        records_mock = Mock()
        mapq_sam_records_filter = MAPQSamRecordsFilter(records_mock)
        assert mapq_sam_records_filter.records_to_keep == set()
        assert get_query_name_to_best_record_mock.called_once_with(records_mock)

    @patch.object(MAPQSamRecordsFilter, MAPQSamRecordsFilter.get_query_name_to_best_record.__name__, return_value={"record_1": "best_record_1"})
    def test___constructor___one_sam_record___records_to_keep_has_sam_record(self, get_query_name_to_best_record_mock):
        records_mock = Mock()
        mapq_sam_records_filter = MAPQSamRecordsFilter(records_mock)
        assert mapq_sam_records_filter.records_to_keep == {"best_record_1"}
        assert get_query_name_to_best_record_mock.called_once_with(records_mock)

    @patch.object(MAPQSamRecordsFilter, MAPQSamRecordsFilter.get_query_name_to_best_record.__name__,
                  return_value={"record_1": "best_record_1",
                                "record_2": "best_record_2",
                                "record_3": "best_record_3"})
    def test___constructor___three_sam_records___records_to_keep_has_the_three_sam_records(self, get_query_name_to_best_record_mock):
        records_mock = Mock()
        mapq_sam_records_filter = MAPQSamRecordsFilter(records_mock)
        assert mapq_sam_records_filter.records_to_keep == {"best_record_1", "best_record_2", "best_record_3"}
        assert get_query_name_to_best_record_mock.called_once_with(records_mock)


    def test___get_record_with_highest_mapping_quality___no_records___returns_None(self):
        mapq_sam_records_filter = MAPQSamRecordsFilter([])
        record_with_highest_mapping_quality = mapq_sam_records_filter.get_record_with_highest_mapping_quality([])
        assert record_with_highest_mapping_quality is None

    def test___get_record_with_highest_mapping_quality___one_record___returns_given_record(self):
        mapq_sam_records_filter = MAPQSamRecordsFilter([])
        record_with_highest_mapping_quality = mapq_sam_records_filter.get_record_with_highest_mapping_quality([mapping_quality_49_mock])
        assert record_with_highest_mapping_quality == mapping_quality_49_mock

    def test___get_record_with_highest_mapping_quality___several_records(self):
        mapq_sam_records_filter = MAPQSamRecordsFilter([])
        record_with_highest_mapping_quality = mapq_sam_records_filter.get_record_with_highest_mapping_quality(
            [mapping_quality_49_mock, mapping_quality_50_mock, mapping_quality_51_mock, mapping_quality_49_mock,
             mapping_quality_49_mock, mapping_quality_51_mock, mapping_quality_60_mock])
        assert record_with_highest_mapping_quality == mapping_quality_60_mock

    def test___get_query_name_to_best_record___no_records___return_empty_dict(self):
        mapq_sam_records_filter = MAPQSamRecordsFilter([])
        actual = mapq_sam_records_filter.get_query_name_to_best_record([])
        assert actual == {}

    def test___get_query_name_to_best_record___one_query_name_with_several_records___difference_exceeds_threshold___returns_best_mapping(self):
        mapq_sam_records_filter = MAPQSamRecordsFilter([], mapping_quality_threshold=10)
        actual = mapq_sam_records_filter.get_query_name_to_best_record([mapping_quality_49_mock, mapping_quality_40_mock,
           mapping_quality_60_mock, mapping_quality_40_mock, mapping_quality_40_mock])
        assert actual == {mapping_quality_60_mock.query_name: mapping_quality_60_mock}

    def test___get_query_name_to_best_record___one_query_name_with_several_records___difference_does_not_exceed_threshold___returns_empty_dict(self):
        mapq_sam_records_filter = MAPQSamRecordsFilter([], mapping_quality_threshold=10)
        actual = mapq_sam_records_filter.get_query_name_to_best_record([mapping_quality_50_mock, mapping_quality_60_mock])
        assert actual == {}

    def test___get_query_name_to_best_record___several_query_names_with_several_records(self):
        query_name_2_mapping_quality_60_mock = Mock(query_name="query_name_2", mapping_quality=60)
        query_name_2_mapping_quality_51_mock = Mock(query_name="query_name_2", mapping_quality=51)
        query_name_3_mapping_quality_60_mock = Mock(query_name="query_name_3", mapping_quality=60)
        query_name_3_mapping_quality_0_mock = Mock(query_name="query_name_3", mapping_quality=0)
        query_name_3_mapping_quality_30_mock = Mock(query_name="query_name_3", mapping_quality=30)
        query_name_4_mapping_quality_60_mock = Mock(query_name="query_name_4", mapping_quality=60)

        mapq_sam_records_filter = MAPQSamRecordsFilter([], mapping_quality_threshold=10)
        actual = mapq_sam_records_filter.get_query_name_to_best_record(
            [   # query_name_1 - exceeds threshold
                mapping_quality_49_mock, mapping_quality_40_mock,
                mapping_quality_60_mock, mapping_quality_40_mock, mapping_quality_40_mock,
                # query_name_2 - does not exceed threshold
                query_name_2_mapping_quality_51_mock, query_name_2_mapping_quality_60_mock,
                # query_name_3 - exceeds threshold
                query_name_3_mapping_quality_60_mock, query_name_3_mapping_quality_0_mock, query_name_3_mapping_quality_30_mock,
                # query_name_4 - does not exceed threshold
                query_name_4_mapping_quality_60_mock, query_name_4_mapping_quality_60_mock, query_name_4_mapping_quality_60_mock
             ])
        assert actual == {"query_name_1": mapping_quality_60_mock,
                          "query_name_3": query_name_3_mapping_quality_60_mock}

    @patch.object(MAPQSamRecordsFilter, "records_to_keep", new_callable=PropertyMock,
                  return_value=[mapping_quality_60_mock, mapping_quality_51_mock])
    def test___record_should_be_filtered_out___record_is_in_records_to_keep___should_NOT_be_filtered_out(self, *mocks):
        mapq_sam_records_filter = MAPQSamRecordsFilter([])
        assert not mapq_sam_records_filter.record_should_be_filtered_out(mapping_quality_51_mock)

    @patch.object(MAPQSamRecordsFilter, "records_to_keep", new_callable=PropertyMock,
                  return_value=[mapping_quality_60_mock, mapping_quality_51_mock])
    def test___record_should_be_filtered_out___record_is_NOT_in_records_to_keep___should_be_filtered_out(self, *mocks):
        mapq_sam_records_filter = MAPQSamRecordsFilter([])
        assert mapq_sam_records_filter.record_should_be_filtered_out(mapping_quality_50_mock)
