from evaluate.probe import *
import pytest


class TestInterval:
    def test_isNull_nullIntervalReturnsTrue(self):
        assert Interval(-1, -1).is_null()

    def test_isNull_nonNullIntervalReturnsFalse(self):
        assert not Interval(0, -1).is_null()

    def test_fromString_emptyStringReturnsNullInterval(self):
        string = ""

        actual = Interval.from_string(string)

        assert actual.is_null()

    def test_fromString_invalidStringReturnsNullInterval(self):
        string = "[10,20]"

        actual = Interval.from_string(string)

        assert actual.is_null()

    def test_fromString_validStringWithSpaceAfterCommaReturnsNullInterval(self):
        string = "[10, 20)"

        actual = Interval.from_string(string)

        assert actual.is_null()

    def test_fromString_validStringReturnsInterval(self):
        string = "[10,20)"

        actual = Interval.from_string(string)
        expected = Interval(10, 20)

        assert actual == expected


class TestProbeHeader:
    def test_equality_equalReturnsTrue(self):
        p1 = ProbeHeader(sample="foo")
        p2 = ProbeHeader(sample="foo")

        assert p1 == p2

    def test_equality_notEqualReturnsFalse(self):
        p1 = ProbeHeader(interval=Interval(2, 5))
        p2 = ProbeHeader(interval=Interval(2, 4))

        assert p1 != p2

    def test_fromString_emptyStringReturnsEmptyProbeHeader(self):
        string = ""

        actual = ProbeHeader.from_string(string)
        expected = ProbeHeader()

        assert actual == expected

    def test_fromString_allFieldsInStringReturnsProbeHeaderWithAllFields(self):
        string = ">CHROM=1;SAMPLE=CFT073;POS=1;INTERVAL=[0,72);SVTYPE=INDEL;MEAN_FWD_COVG=2;MEAN_REV_COVG=3;GT_CONF=10.9922;"

        actual = ProbeHeader.from_string(string)
        expected = ProbeHeader(
            sample="CFT073",
            chrom="1",
            pos=1,
            interval=Interval(0, 72),
            svtype="INDEL",
            mean_fwd_covg=2,
            mean_rev_covg=3,
            gt_conf=10.9922,
        )

        assert actual == expected

    def test_fromString_someFieldsInStringReturnsProbeHeaderWithSomeFields(self):
        string = ">CHROM=1;SAMPLE=CFT073;SVTYPE=INDEL;MEAN_FWD_COVG=2;MEAN_REV_COVG=3;GT_CONF=10.9922;"

        actual = ProbeHeader.from_string(string)
        expected = ProbeHeader(
            sample="CFT073",
            chrom="1",
            svtype="INDEL",
            mean_fwd_covg=2,
            mean_rev_covg=3,
            gt_conf=10.9922,
        )

        assert actual == expected
    def test_fromString_stringWithInvalidFieldReturnsProbeHeaderWithOnlyValidFields(self):
        string = ">CHROM=1;SAMPLE=CFT073;SVTYPE=INDEL;MEAN_FWD_COVG=2;INVALID=foo;MEAN_REV_COVG=3;GT_CONF=10.9922;"

        actual = ProbeHeader.from_string(string)
        expected = ProbeHeader(
            sample="CFT073",
            chrom="1",
            svtype="INDEL",
            mean_fwd_covg=2,
            mean_rev_covg=3,
            gt_conf=10.9922,
        )

        assert actual == expected
