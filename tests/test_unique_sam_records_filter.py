from unittest.mock import patch, PropertyMock, Mock
from evaluate.unique_sam_records_filter import UniqueSamRecordsFilter


query_name_1_mock = Mock(query_name = "query_name_1")
query_name_2_mock = Mock(query_name = "query_name_2")
query_name_3_mock = Mock(query_name = "query_name_3")
query_name_4_mock = Mock(query_name = "query_name_4")

class TestUniqueSamRecordsFilter:
    def test___constructor___no_sam_records___all_unique_query_names_is_empty(self):
        records = []
        unique_sam_records_filter = UniqueSamRecordsFilter(records)
        assert unique_sam_records_filter.all_unique_query_names == []

    def test___constructor___one_unique_sam_record___all_unique_query_names_contains_one_record(self):
        records = [query_name_1_mock]
        unique_sam_records_filter = UniqueSamRecordsFilter(records)
        assert unique_sam_records_filter.all_unique_query_names == [query_name_1_mock.query_name]

    def test___constructor___one_duplicated_sam_record___all_unique_query_names_is_empty(self):
        records = [query_name_1_mock, query_name_1_mock]
        unique_sam_records_filter = UniqueSamRecordsFilter(records)
        assert unique_sam_records_filter.all_unique_query_names == []

    def test___constructor___one_quintuplicated_sam_record___all_unique_query_names_is_empty(self):
        records = [query_name_1_mock, query_name_1_mock, query_name_1_mock, query_name_1_mock, query_name_1_mock]
        unique_sam_records_filter = UniqueSamRecordsFilter(records)
        assert unique_sam_records_filter.all_unique_query_names == []

    def test___constructor___one_unique_and_one_duplicated_sam_record___all_unique_query_names_contains_one_record(self):
        records = [query_name_1_mock, query_name_2_mock, query_name_1_mock]
        unique_sam_records_filter = UniqueSamRecordsFilter(records)
        assert unique_sam_records_filter.all_unique_query_names == [query_name_2_mock.query_name]

    def test___constructor___two_uniques_and_one_duplicated_sam_record___all_unique_query_names_contains_two_records(self):
        records = [query_name_3_mock, query_name_1_mock, query_name_2_mock, query_name_1_mock]
        unique_sam_records_filter = UniqueSamRecordsFilter(records)
        assert unique_sam_records_filter.all_unique_query_names == [query_name_2_mock.query_name, query_name_3_mock.query_name]

    def test___constructor___one_unique_and_two_duplicated_sam_record___all_unique_query_names_contains_one_record(self):
        records = [query_name_3_mock, query_name_2_mock, query_name_3_mock, query_name_2_mock, query_name_1_mock]
        unique_sam_records_filter = UniqueSamRecordsFilter(records)
        assert unique_sam_records_filter.all_unique_query_names == [query_name_1_mock.query_name]

    def test___constructor___two_unique_and_two_duplicated_sam_record___all_unique_query_names_contains_the_two_uniques(self):
        records = [query_name_3_mock, query_name_3_mock, query_name_2_mock, query_name_4_mock, query_name_3_mock, query_name_2_mock, query_name_1_mock]
        unique_sam_records_filter = UniqueSamRecordsFilter(records)
        assert unique_sam_records_filter.all_unique_query_names == [query_name_1_mock.query_name, query_name_4_mock.query_name]

    @patch.object(UniqueSamRecordsFilter, "all_unique_query_names", new_callable=PropertyMock,
                  return_value=[query_name_1_mock.query_name, query_name_2_mock.query_name, query_name_3_mock.query_name])
    def test___record_should_be_filtered_out___record_is_not_in_all_unique_query_names___should_be_filtered_out(self, *mocks):
        unique_sam_records_filter = UniqueSamRecordsFilter([])
        assert unique_sam_records_filter.record_should_be_filtered_out(query_name_4_mock)

    @patch.object(UniqueSamRecordsFilter, "all_unique_query_names", new_callable=PropertyMock,
                  return_value=[query_name_1_mock.query_name, query_name_2_mock.query_name, query_name_3_mock.query_name])
    def test___record_should_be_filtered_out___record_is_in_all_unique_query_names___should_NOT_be_filtered_out(self, *mocks):
        unique_sam_records_filter = UniqueSamRecordsFilter([])
        assert not unique_sam_records_filter.record_should_be_filtered_out(query_name_2_mock)