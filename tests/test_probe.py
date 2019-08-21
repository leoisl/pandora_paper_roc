from evaluate.probe import *


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

    def test_fromString_stringWithInvalidFieldReturnsProbeHeaderWithOnlyValidFields(
        self
    ):
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


class TestProbe:
    def test_equality_twoEqualProbesReturnsTrue(self):
        p1 = Probe(
            ProbeHeader(sample="foo", interval=Interval(1, 2)), full_sequence="bar"
        )
        p2 = Probe(
            ProbeHeader(sample="foo", interval=Interval(1, 2)), full_sequence="bar"
        )

        assert p1 == p2

    def test_equality_twoNonEqualProbesReturnsFalse(self):
        p1 = Probe(
            ProbeHeader(sample="foo", interval=Interval(1, 2)), full_sequence="barr"
        )
        p2 = Probe(
            ProbeHeader(sample="foo", interval=Interval(1, 2)), full_sequence="bar"
        )

        assert p1 != p2

    def test_getLeftFlank_emptyFullSequenceReturnsEmptyFlank(self):
        header = ProbeHeader(interval=Interval(4, 5))
        full_sequence = ""
        probe = Probe(header, full_sequence)

        actual = probe.get_left_flank()
        expected = ""

        assert actual == expected

    def test_getLeftFlank_noLeftFlankInFullSequenceReturnsEmptyFlank(self):
        header = ProbeHeader(interval=Interval(0, 5))
        full_sequence = "abcdefg"
        probe = Probe(header, full_sequence)

        actual = probe.get_left_flank()
        expected = ""

        assert actual == expected

    def test_getLeftFlank_singleBaseLeftFlankInFullSequenceReturnsLeftFlank(self):
        header = ProbeHeader(interval=Interval(1, 5))
        full_sequence = "abcdefg"
        probe = Probe(header, full_sequence)

        actual = probe.get_left_flank()
        expected = "a"

        assert actual == expected

    def test_getLeftFlank_multiBaseLeftFlankInFullSequenceReturnsLeftFlank(self):
        header = ProbeHeader(interval=Interval(3, 5))
        full_sequence = "abcdefg"
        probe = Probe(header, full_sequence)

        actual = probe.get_left_flank()
        expected = "abc"

        assert actual == expected

    def test_getRightFlank_emptyFullSequenceReturnsEmptyFlank(self):
        header = ProbeHeader(interval=Interval(4, 5))
        full_sequence = ""
        probe = Probe(header, full_sequence)

        actual = probe.get_right_flank()
        expected = ""

        assert actual == expected

    def test_getRightFlank_noRightFlankInFullSequenceReturnsEmptyFlank(self):
        header = ProbeHeader(interval=Interval(2, 7))
        full_sequence = "abcdefg"
        probe = Probe(header, full_sequence)

        actual = probe.get_right_flank()
        expected = ""

        assert actual == expected

    def test_getRightFlank_singleBaseRightFlankInFullSequenceReturnsRightFlank(self):
        header = ProbeHeader(interval=Interval(1, 6))
        full_sequence = "abcdefg"
        probe = Probe(header, full_sequence)

        actual = probe.get_right_flank()
        expected = "g"

        assert actual == expected

    def test_getRightFlank_multiBaseRightFlankInFullSequenceReturnsRightFlank(self):
        header = ProbeHeader(interval=Interval(3, 4))
        full_sequence = "abcdefg"
        probe = Probe(header, full_sequence)

        actual = probe.get_right_flank()
        expected = "efg"

        assert actual == expected
