from pathlib import Path

import pytest

from evaluate.bwa import *

TEST_CASES = Path("tests/test_cases")
TEST_PANEL = TEST_CASES / "test_panel.fa"


class TestBwa:
    def test_align_refNotIndexed_raiseIndexError(self):
        with pytest.raises(IndexError) as excinfo:
            bwa = BWA()
            query = ">test\nTACGACACAGTGACGACATAAC"
            query_path = Path("tests/test_cases/TestBwa.1.fa")
            query_path.write_text(query)
            bwa.align(query_path)

        assert "indexed by BWA" in str(excinfo.value)

    def test_align_validQuery_returnValidSam(self):
        bwa = BWA()
        bwa.index(TEST_PANEL)
        query = ">test\nGACGTTAAATGCAAAAATCGCACGTCTTGAGCAGGATATAAAA"
        query_path = Path("tests/test_cases/TestBwa.2.fa")
        query_path.write_text(query)
        stdout, stderr = bwa.align(query_path)
        header, sam = bwa.parse_sam_string(stdout)
        expected = "test\t0\tC15154T\t55\t60\t43M\t*\t0\t0\tGACGTTAAATGCAAAAATCGCACGTCTTGAGCAGGATATAAAA\t*\tNM:i:0\tMD:Z:43\tAS:i:43\tXS:i:0"
        actual = sam[0].to_string()

        assert len(sam) == 1
        assert actual == expected

    def test_align_validQuery_returnValidHeader(self):
        bwa = BWA()
        bwa.index(TEST_PANEL)
        query = ">test\nGACGTTAAATGCAAAAATCGCACGTCTTGAGCAGGATATAAAA"
        query_path = Path("tests/test_cases/TestBwa.3.fa")
        query_path.write_text(query)
        stdout, stderr = bwa.align(query_path)
        header, sam = bwa.parse_sam_string(stdout)
        expected_start = (
            "@SQ\tSN:C15154T\tLN:201\n@SQ\tSN:T16509G\tLN:201\n@PG\tID:bwa\tPN:bwa\tVN:"
        )
        expected_end = f"CL:bwa mem -t 1 {TEST_PANEL} {query_path}\n"
        actual = str(header)

        # have to do it this way as the bwa version number is in the header and may vary
        assert actual.startswith(expected_start)
        assert actual.endswith(expected_end)

    def test_getOptions_defaultOneThread(self):
        bwa = BWA()

        actual = bwa.get_options()
        expected = ["-t", "1"]

        assert actual == expected

    def test_getOptions_setThreeThreads_returnThreeThreads(self):
        bwa = BWA(3)

        actual = bwa.get_options()
        expected = ["-t", "3"]

        assert actual == expected
