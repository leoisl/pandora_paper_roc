from evaluate.bwa import *
from pathlib import Path
import pytest


TEST_CASES = Path("tests/test_cases")
TEST_PANEL = TEST_CASES / "test_panel.fa"


class TestBwa:
    def test_align_refNotIndexed_raiseIndexError(self):
        with pytest.raises(IndexError) as excinfo:
            bwa = BWA()
            query = ">test\nTACGACACAGTGACGACATAAC"
            bwa.align(query)

        assert "indexed by BWA" in str(excinfo.value)

    def test_align_validQuery_returnValidSam(self):
        bwa = BWA()
        bwa.index(TEST_PANEL)
        query = ">test\nGACGTTAAATGCAAAAATCGCACGTCTTGAGCAGGATATAAAA"
        _, sam = bwa.align(query)
        expected = "test\t0\tC15154T\t55\t60\t43M\t*\t0\t0\tGACGTTAAATGCAAAAATCGCACGTCTTGAGCAGGATATAAAA\t*\tNM:i:0\tMD:Z:43\tAS:i:43\tXS:i:0"
        actual = sam[0].to_string()

        assert len(sam) == 1
        assert actual == expected

    def test_align_validQuery_returnValidHeader(self):
        bwa = BWA()
        bwa.index(TEST_PANEL)
        query = ">test\nGACGTTAAATGCAAAAATCGCACGTCTTGAGCAGGATATAAAA"
        header, _ = bwa.align(query)
        expected_start = (
            "@SQ\tSN:C15154T\tLN:201\n@SQ\tSN:T16509G\tLN:201\n@PG\tID:bwa\tPN:bwa\tVN:"
        )
        expected_end = f"CL:bwa mem -t 1 {TEST_PANEL} -\n"
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
