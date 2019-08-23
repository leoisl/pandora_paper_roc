from evaluate.probe import *


class TestInterval:
    def test_str_nullIntervalReturnsEmptyString(self):
        interval = Interval()

        actual = str(interval)
        expected = ""

        assert actual == expected

    def test_str_nonNullIntervalReturnsExpectedString(self):
        interval = Interval(1, 5)

        actual = str(interval)
        expected = "[1,5)"

        assert actual == expected

    def test_len_isZeroLengthReturnsZero(self):
        interval = Interval(2, 2)

        actual = len(interval)
        expected = 0

        assert actual == expected

    def test_len_isTwoReturnsTwo(self):
        interval = Interval(4, 6)

        actual = len(interval)
        expected = 2

        assert actual == expected

    def test_len_isNullReturnsZero(self):
        interval = Interval()

        actual = len(interval)
        expected = 0

        assert actual == expected

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

    def test_str_emptyProbeHeaderReturnsEmptyString(self):
        header = ProbeHeader()

        actual = str(header)
        expected = ""

        assert actual == expected

    def test_str_singleVariableProbeHeaderReturnsStringWithOneField(self):
        header = ProbeHeader(mean_rev_covg=9)

        actual = str(header)
        expected = ">MEAN_REV_COVG=9;"

        assert actual == expected

    def test_str_allFieldsInProbeHeaderReturnsStringWithAllFields(self):
        header = ProbeHeader(
            mean_rev_covg=9,
            mean_fwd_covg=9,
            pos=3,
            gt_conf=2.2,
            chrom="4",
            sample="foo",
            interval=Interval(5, 6),
            svtype="SNP",
        )

        actual = str(header)
        expected = ">CHROM=4;SAMPLE=foo;POS=3;INTERVAL=[5,6);SVTYPE=SNP;MEAN_FWD_COVG=9;MEAN_REV_COVG=9;GT_CONF=2.2;"

        assert actual == expected

        assert ProbeHeader.from_string(actual) == header


class TestProbe:
    def test_str_emptyProbeReturnsEmptyString(self):
        probe = Probe()

        actual = str(probe)
        expected = ""

        assert actual == expected

    def test_str_emptyFullSequenceReturnsHeaderWithNewline(self):
        probe = Probe(ProbeHeader(chrom="3"))

        actual = str(probe)
        expected = ">CHROM=3;\n"

        assert actual == expected

    def test_str_fullProbeReturnsHeaderAndSequence(self):
        probe = Probe(ProbeHeader(pos=4, chrom="3"), full_sequence="foo")

        actual = str(probe)
        expected = ">CHROM=3;POS=4;\nfoo"

        assert actual == expected

        assert Probe.from_string(actual) == probe

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

        actual = probe.left_flank
        expected = ""

        assert actual == expected

    def test_getLeftFlank_noLeftFlankInFullSequenceReturnsEmptyFlank(self):
        header = ProbeHeader(interval=Interval(0, 5))
        full_sequence = "abcdefg"
        probe = Probe(header, full_sequence)

        actual = probe.left_flank
        expected = ""

        assert actual == expected

    def test_getLeftFlank_singleBaseLeftFlankInFullSequenceReturnsLeftFlank(self):
        header = ProbeHeader(interval=Interval(1, 5))
        full_sequence = "abcdefg"
        probe = Probe(header, full_sequence)

        actual = probe.left_flank
        expected = "a"

        assert actual == expected

    def test_getLeftFlank_multiBaseLeftFlankInFullSequenceReturnsLeftFlank(self):
        header = ProbeHeader(interval=Interval(3, 5))
        full_sequence = "abcdefg"
        probe = Probe(header, full_sequence)

        actual = probe.left_flank
        expected = "abc"

        assert actual == expected

    def test_getRightFlank_emptyFullSequenceReturnsEmptyFlank(self):
        header = ProbeHeader(interval=Interval(4, 5))
        full_sequence = ""
        probe = Probe(header, full_sequence)

        actual = probe.right_flank
        expected = ""

        assert actual == expected

    def test_getRightFlank_noRightFlankInFullSequenceReturnsEmptyFlank(self):
        header = ProbeHeader(interval=Interval(2, 7))
        full_sequence = "abcdefg"
        probe = Probe(header, full_sequence)

        actual = probe.right_flank
        expected = ""

        assert actual == expected

    def test_getRightFlank_singleBaseRightFlankInFullSequenceReturnsRightFlank(self):
        header = ProbeHeader(interval=Interval(1, 6))
        full_sequence = "abcdefg"
        probe = Probe(header, full_sequence)

        actual = probe.right_flank
        expected = "g"

        assert actual == expected

    def test_getRightFlank_multiBaseRightFlankInFullSequenceReturnsRightFlank(self):
        header = ProbeHeader(interval=Interval(3, 4))
        full_sequence = "abcdefg"
        probe = Probe(header, full_sequence)

        actual = probe.right_flank
        expected = "efg"

        assert actual == expected

    def test_getCoreSequence_emptyFullSequenceReturnsEmptyString(self):
        header = ProbeHeader(interval=Interval(3, 4))
        full_sequence = ""
        probe = Probe(header, full_sequence)

        actual = probe.core_sequence
        expected = ""

        assert actual == expected

    def test_getCoreSequence_intervalForDeletionReturnsEmptyString(self):
        header = ProbeHeader(interval=Interval(3, 3))
        full_sequence = "abcdefg"
        probe = Probe(header, full_sequence)

        actual = probe.core_sequence
        expected = ""

        assert actual == expected

    def test_getCoreSequence_intervalForSingleBaseReturnsSingleBase(self):
        header = ProbeHeader(interval=Interval(3, 4))
        full_sequence = "abcdefg"
        probe = Probe(header, full_sequence)

        actual = probe.core_sequence
        expected = "d"

        assert actual == expected

    def test_getCoreSequence_intervalForSingleBaseAtStartOfSequenceReturnsSingleBase(
        self
    ):
        header = ProbeHeader(interval=Interval(0, 1))
        full_sequence = "abcdefg"
        probe = Probe(header, full_sequence)

        actual = probe.core_sequence
        expected = "a"

        assert actual == expected

    def test_getCoreSequence_intervalForSingleBaseAtEndOfSequenceReturnsSingleBase(
        self
    ):
        header = ProbeHeader(interval=Interval(6, 7))
        full_sequence = "abcdefg"
        probe = Probe(header, full_sequence)

        actual = probe.core_sequence
        expected = "g"

        assert actual == expected

    def test_getCoreSequence_intervalForMultiBaseSequenceReturnsMultiBase(self):
        header = ProbeHeader(interval=Interval(2, 5))
        full_sequence = "abcdefg"
        probe = Probe(header, full_sequence)

        actual = probe.core_sequence
        expected = "cde"

        assert actual == expected

    def test_getCoreSequence_intervalOutOfRangeReturnsEmptyString(self):
        header = ProbeHeader(interval=Interval(8, 9))
        full_sequence = "abcdefg"
        probe = Probe(header, full_sequence)

        actual = probe.core_sequence
        expected = ""

        assert actual == expected

    def test_fromString_emptyStringReturnsEmptyProbe(self):
        string = ""

        actual = Probe.from_string(string)
        expected = Probe()

        assert actual == expected

    def test_fromString_headerOnlyStringReturnsProbeWithNoFullSequence(self):
        string = ">CHROM=1;INTERVAL=[3,5);"

        actual = Probe.from_string(string)
        expected = Probe(header=ProbeHeader(chrom="1", interval=Interval(3, 5)))

        assert actual == expected

    def test_fromString_headerAndEmptySequenceInStringReturnsProbeWithNoFullSequence(
        self
    ):
        string = ">CHROM=1;INTERVAL=[3,5);\n"

        actual = Probe.from_string(string)
        expected = Probe(header=ProbeHeader(chrom="1", interval=Interval(3, 5)))

        assert actual == expected

    def test_fromString_headerAndSequenceInStringReturnsFullProbe(self):
        string = ">CHROM=1;INTERVAL=[3,5);\nfoo"

        actual = Probe.from_string(string)
        expected = Probe(
            header=ProbeHeader(chrom="1", interval=Interval(3, 5)), full_sequence="foo"
        )

        assert actual == expected

    def test_gtConf(self):
        probe = Probe(header=ProbeHeader(gt_conf=5.5))

        actual = probe.gt_conf
        expected = 5.5

        assert actual == expected

    def test_interval(self):
        probe = Probe(header=ProbeHeader(interval=Interval(4, 6)))

        actual = probe.interval
        expected = Interval(4, 6)

        assert actual == expected