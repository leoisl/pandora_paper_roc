from unittest.mock import patch, MagicMock
from evaluate.filter import Filter


class TestFilter:
    def test_filterRecords_noRecordsReturnsEmpty(self):
        records = []
        filter = Filter()

        actual = filter.filter_records(records)
        expected = []

        assert actual == expected

    @patch.object(
        Filter, Filter.record_should_be_filtered_out.__name__, return_value=False
    )
    def test_filterRecords_oneRecordNotFilteredOutReturnsRecord(self, *mock):
        record = MagicMock()
        records = [record]
        filter = Filter()

        actual = filter.filter_records(records)
        expected = [record]

        assert actual == expected

    @patch.object(
        Filter, Filter.record_should_be_filtered_out.__name__, return_value=True
    )
    def test_filterRecords_oneRecordFilteredOutReturnsEmpty(self, *mock):
        record = MagicMock()
        records = [record]
        filter = Filter()

        actual = filter.filter_records(records)
        expected = []

        assert actual == expected

    @patch.object(
        Filter, Filter.record_should_be_filtered_out.__name__, side_effect=[True, False]
    )
    def test_filterRecords_twoRecordsOneFilteredOutOneNotFilteredOutReturnsOneRecord(
        self, *mock
    ):
        record_filtered_out = MagicMock()
        record_not_filtered_out = MagicMock()
        records = [record_filtered_out, record_not_filtered_out]
        filter = Filter()

        actual = filter.filter_records(records)
        expected = [record_not_filtered_out]

        assert actual == expected

    @patch.object(
        Filter, Filter.record_should_be_filtered_out.__name__, return_value=False
    )
    def test_filterRecords_twoRecordsNoneFilteredOutReturnsTwoRecords(self, *mock):
        record_1 = MagicMock()
        record_2 = MagicMock()
        records = [record_1, record_2]
        filter = Filter()

        actual = filter.filter_records(records)
        expected = records

        assert actual == expected

    @patch.object(
        Filter, Filter.record_should_be_filtered_out.__name__, return_value=True
    )
    def test_filterRecords_twoRecordsBothFilteredOutReturnsNoRecords(self, *mock):
        record_1 = MagicMock()
        record_2 = MagicMock()
        records = [record_1, record_2]
        filter = Filter()

        actual = filter.filter_records(records)
        expected = []

        assert actual == expected
