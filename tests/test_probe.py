from evaluate.probe import *
import pytest


class TestInterval:
    def test_fromString_emptyStringRaisesRegexError(self):
        string = ""

        with pytest.raises(RegexError):
            Interval.from_string(string)

    def test_fromString_invalidStringRaisesRegexError(self):
        string = "[10,20]"

        with pytest.raises(RegexError):
            Interval.from_string(string)

    def test_fromString_validStringWithSpaceAfterCommaRaisesRegexError(self):
        string = "[10, 20)"

        with pytest.raises(RegexError):
            Interval.from_string(string)

    def test_fromString_validStringReturnsInterval(self):
        string = "[10,20)"

        actual = Interval.from_string(string)
        expected = Interval(10, 20)

        assert actual == expected


class TestProbe:
    def test_equality_equalReturnsTrue(self):
        p1 = Probe(sample="foo")
        p2 = Probe(sample="foo")

        assert p1 == p2

    def test_equality_notEqualReturnsFalse(self):
        p1 = Probe(interval=Interval(2, 5))
        p2 = Probe(interval=Interval(2, 4))

        assert p1 != p2

    def test_fromString_emptyStringReturnsEmptyProbe(self):
        string = ""

        actual = Probe.from_string(string)
        expected = Probe()

        assert actual == expected

    def test_fromString_headerOnlyStringReturnsProbeWithNoSequence(self):
        string = ">CHROM=1;SAMPLE=CFT073;POS=1;INTERVAL=[0,72);SVTYPE=INDEL;MEAN_FWD_COVG=2;MEAN_REV_COVG=3;GT_CONF=10.9922;"

        actual = Probe.from_string(string)
        expected = Probe(
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
