import uuid
from io import StringIO

import pandas.errors
import pytest

from evaluate.mummer import *

REF = Path("tests/test_cases/ref.fa")
QUERY = Path("tests/test_cases/query.fa")
DELTA = Path("tests/test_cases/out.delta")
DELTA1 = Path("tests/test_cases/out.delta1")


class TestNucmer:
    def test_init_withFakeFileRaisesFileNotFoundError(self):
        reference = Path("foo")
        query = Path("bar")
        with pytest.raises(NucmerError):
            nucmer = Nucmer(reference, query)

    def test_generateCommand_defaultArgsHasNoExtraParamsAndHasOutAsPrefix(self):
        reference = REF
        query = QUERY
        nucmer = Nucmer(reference, query)

        actual = nucmer.generate_command()
        expected = ["nucmer", "--prefix", "out", str(reference), str(query)]

        assert actual == expected

    def test_generateCommand_extraParamsIncluded(self):
        reference = REF
        query = QUERY
        extra_params = "--maxmatch --forward"
        nucmer = Nucmer(reference, query, extra_params=extra_params)

        actual = nucmer.generate_command()
        expected = [
            "nucmer",
            extra_params,
            "--prefix",
            "out",
            str(reference),
            str(query),
        ]

        assert actual == expected

    def test_run_twoTestCases_returnsOkExitCode(self):
        reference = REF
        query = QUERY
        prefix = f"/tmp/{str(uuid.uuid4())}"
        extra_params = "--maxmatch"
        nucmer = Nucmer(reference, query, prefix=prefix, extra_params=extra_params)
        result = nucmer.run()

        actual = result.returncode
        expected = 0

        assert actual == expected

        Path(prefix + ".delta").unlink()


class TestDeltaFilter:
    def test_init_withFakeFileRaisesFileNotFoundError(self):
        deltafile = Path("foo")

        with pytest.raises(DeltaFilterError):
            deltafilter = DeltaFilter(deltafile)

    def test_generateCommand_defaultArgsHasNoExtraParamsAndHasOutAsPrefix(self):
        deltafile = DELTA
        deltafilter = DeltaFilter(deltafile)

        actual = deltafilter.generate_command()
        expected = ["delta-filter", str(deltafile)]

        assert actual == expected

    def test_generateCommand_extraParamsIncluded(self):
        deltafile = DELTA
        extra_params = "-1"
        deltafilter = DeltaFilter(deltafile, extra_params=extra_params)

        actual = deltafilter.generate_command()
        expected = ["delta-filter", extra_params, str(deltafile)]

        assert actual == expected

    def test_run_realDeltaFileReturnsOkExitCode(self):
        deltafile = DELTA
        extra_params = "-1"
        deltafilter = DeltaFilter(deltafile, extra_params=extra_params)
        result = deltafilter.run()

        actual = result.returncode
        expected = 0

        assert actual == expected

        actual = result.stdout.decode()
        expected = """tests/test_cases/ref.fa tests/test_cases/query.fa
NUCMER
>ref query 85 84
1 85 1 84 2 2 0
39
0
"""

        assert actual == expected


class TestShowSnps:
    def test_init_withFakeFileRaisesFileNotFoundError(self):
        deltafile = Path("foo")

        with pytest.raises(ShowSnpsError):
            showsnps = ShowSnps(deltafile)

    def test_generateCommand_defaultArgsHasNoExtraParams(self):
        deltafile = DELTA1
        showsnps = ShowSnps(deltafile)

        actual = showsnps.generate_command()
        expected = ["show-snps", str(deltafile)]

        assert actual == expected

    def test_generateCommand_oppositeDefaultArgsAndAmbiguousMappingExtraParam(self):
        deltafile = DELTA1
        print_header = False
        indels = False
        context = 3
        extra_params = "-C"
        showsnps = ShowSnps(
            deltafile,
            context=context,
            print_header=print_header,
            indels=indels,
            extra_params=extra_params,
        )

        actual = showsnps.generate_command()
        expected = ["show-snps", "-C", f"-x {context}", "-H", "-I", str(deltafile)]

        assert actual == expected

    def test_run_realDeltaFileReturnsOkExitCodeAndExpectedOutput(self):
        deltafile = DELTA1
        context = 3
        extra_params = "-rlTC"
        showsnps = ShowSnps(deltafile, context=context, extra_params=extra_params)
        result = showsnps.run()

        actual = result.returncode
        expected = 0

        assert actual == expected

        actual = result.stdout.decode()
        expected = """tests/test_cases/ref.fa tests/test_cases/query.fa
NUCMER

[P1]\t[SUB]\t[SUB]\t[P2]\t[BUFF]\t[DIST]\t[LEN R]\t[LEN Q]\t[CTX R]\t[CTX Q]\t[FRM]\t[TAGS]
39\tG\t.\t38\t34\t38\t85\t84\tGTAGTAG\tGTA.TAG\t1\t1\tref\tquery
73\tT\tA\t72\t13\t13\t85\t84\tGGATTGA\tGGAATGA\t1\t1\tref\tquery
"""

        assert actual == expected

    def test_toDataframe_emptyFileReturnsEmpty(self):
        snps = StringIO()

        actual = ShowSnps.to_dataframe(snps)

        assert actual.empty

    def test_toDataframe_invalidInputRaisesError(self):
        snps = StringIO("foo\nbar\nsome\nrandom\ntext\n")

        with pytest.raises(ValueError):
            ShowSnps.to_dataframe(snps)

    def test_toDataFrame_validInputReturnCorrectDataframe(self):
        snps = StringIO(
            """/home/michael/Projects/pandora1_paper/tests/test_cases/ref.fa /home/michael/Projects/pandora1_paper/tests/test_cases/query.fa
NUCMER

[P1]\t[SUB]\t[SUB]\t[P2]\t[BUFF]\t[DIST]\t[LEN R]\t[LEN Q]\t[CTX R]\t[CTX Q]\t[FRM]\t[TAGS]
39\tG\t.\t38\t34\t38\t85\t84\tGTAGTAG\tGTA.TAG\t1\t1\tref\tquery
73\tT\tA\t72\t13\t13\t85\t84\tGGATTGA\tGGAATGA\t1\t1\tref\tquery
"""
        )

        actual = ShowSnps.to_dataframe(snps)
        expected = ShowSNPsDataframe(
            {
                "ref_pos": [39, 73],
                "ref_sub": ["G", "T"],
                "query_sub": [".", "A"],
                "query_pos": [38, 72],
                "nearest_mismatch": [34, 13],
                "nearest_end": [38, 13],
                "ref_len": [85, 85],
                "query_len": [84, 84],
                "ref_context": ["GTAGTAG", "GGATTGA"],
                "query_context": ["GTA.TAG", "GGAATGA"],
                "ref_strand": [1, 1],
                "query_strand": [1, 1],
                "ref_chrom": ["ref", "ref"],
                "query_chrom": ["query", "query"],
            }
        )

        assert actual.equals(expected)


class TestShowSNPsDataframe:
    ref_length = 20
    ref_pos = 3
    rev_ref_pos = 18
    query_length = 10
    query_pos = 8
    rev_query_pos = 3

    @staticmethod
    def create_showssnps_content(
        ref_strand_is_reversed: bool, query_strand_is_reversed: bool
    ) -> StringIO:
        ref_strand = -1 if ref_strand_is_reversed else 1
        query_strand = -1 if query_strand_is_reversed else 1
        return StringIO(
            f"whatever whatever\n"
            f"NUCMER\n"
            f"\n"
            f"[P1]	[SUB]	[SUB]	[P2]	[BUFF]	[DIST]	[LEN R]	[LEN Q]	[CTX R]	[CTX Q]	[FRM]	[TAGS]\n"
            f"{TestShowSNPsDataframe.ref_pos}	G	C	{TestShowSNPsDataframe.query_pos}	1	1	{TestShowSNPsDataframe.ref_length}"
            f"	{TestShowSNPsDataframe.query_length}	AA	AA	{ref_strand}	{query_strand}	ref	query"
        )

    @staticmethod
    def create_expected_showsnps_dataframe(
        ref_strand_is_reversed: bool, query_strand_is_reversed: bool
    ) -> ShowSNPsDataframe:
        if ref_strand_is_reversed:
            ref_pos = TestShowSNPsDataframe.rev_ref_pos
        else:
            ref_pos = TestShowSNPsDataframe.ref_pos

        if query_strand_is_reversed:
            query_pos = TestShowSNPsDataframe.rev_query_pos
        else:
            query_pos = TestShowSNPsDataframe.query_pos

        return ShowSNPsDataframe(
            {
                "ref_pos": [ref_pos],
                "ref_sub": ["G"],
                "query_sub": ["C"],
                "query_pos": [query_pos],
                "nearest_mismatch": [1],
                "nearest_end": [1],
                "ref_len": [TestShowSNPsDataframe.ref_length],
                "query_len": [TestShowSNPsDataframe.query_length],
                "ref_context": ["AA"],
                "query_context": ["AA"],
                "ref_strand": [1],
                "query_strand": [1],
                "ref_chrom": ["ref"],
                "query_chrom": ["query"],
            }
        )

    def test_translateToFWDStrand_queryStrandIsForwardRefStrandIsForwardReturnsPositionsWithNoChange(
        self
    ):
        showsnps_content = TestShowSNPsDataframe.create_showssnps_content(
            ref_strand_is_reversed=False, query_strand_is_reversed=False
        )
        df = ShowSnps.to_dataframe(showsnps_content)

        actual = df.translate_to_FWD_strand()
        expected = TestShowSNPsDataframe.create_expected_showsnps_dataframe(
            ref_strand_is_reversed=False, query_strand_is_reversed=False
        )

        assert actual.equals(expected)

    def test_translateToFWDStrand_queryStrandIsReverseRefStrandIsForwardReturnsPositionsWithQueryChanged(
        self
    ):
        showsnps_content = TestShowSNPsDataframe.create_showssnps_content(
            ref_strand_is_reversed=False, query_strand_is_reversed=True
        )
        df = ShowSnps.to_dataframe(showsnps_content)

        actual = df.translate_to_FWD_strand()
        expected = TestShowSNPsDataframe.create_expected_showsnps_dataframe(
            ref_strand_is_reversed=False, query_strand_is_reversed=True
        )

        assert actual.equals(expected)

    def test_translateToFWDStrand_queryStrandIsForwardRefStrandIsReverseReturnsPositionsWithRefChanged(
        self
    ):
        showsnps_content = TestShowSNPsDataframe.create_showssnps_content(
            ref_strand_is_reversed=True, query_strand_is_reversed=False
        )
        df = ShowSnps.to_dataframe(showsnps_content)

        actual = df.translate_to_FWD_strand()
        expected = TestShowSNPsDataframe.create_expected_showsnps_dataframe(
            ref_strand_is_reversed=True, query_strand_is_reversed=False
        )

        assert actual.equals(expected)

    def test_translateToFWDStrand_queryStrandIsReverseRefStrandIsReverseReturnsBothPositionsChanged(
        self
    ):
        showsnps_content = TestShowSNPsDataframe.create_showssnps_content(
            ref_strand_is_reversed=True, query_strand_is_reversed=True
        )
        df = ShowSnps.to_dataframe(showsnps_content)

        actual = df.translate_to_FWD_strand()
        expected = TestShowSNPsDataframe.create_expected_showsnps_dataframe(
            ref_strand_is_reversed=True, query_strand_is_reversed=True
        )

        assert actual.equals(expected)

    def test_makeTruthPanelFromSnpsDataframe_emptyDataframeReturnsEmptyPanel(self):
        df = ShowSnps.to_dataframe(
            StringIO(
                """/home/michael/Projects/pandora1_paper/tests/test_cases/ref.fa /home/michael/Projects/pandora1_paper/tests/test_cases/query.fa
    NUCMER

    [P1]\t[SUB]\t[SUB]\t[P2]\t[BUFF]\t[DIST]\t[LEN R]\t[LEN Q]\t[CTX R]\t[CTX Q]\t[FRM]\t[TAGS]
    """
            )
        )

        actual = df.get_probes()
        expected = ("", "")

        assert actual == expected

    def test_makeTruthPanelFromSnpsDataframe_invalidDataframeRaisesError(self):
        df = pd.DataFrame(
            {"ref_pos": [39, 73], "ref_sub": ["G", "T"], "query_len": [84, 84]}
        )

        with pytest.raises(AttributeError):
            df.get_probes()

    def test_makeTruthPanelFromSnpsDataframe_validDataframeReturnsTwoProbesets(self):
        df = ShowSnps.to_dataframe(
            StringIO(
                """/home/michael/Projects/pandora1_paper/tests/test_cases/ref.fa /home/michael/Projects/pandora1_paper/tests/test_cases/query.fa
    NUCMER

    [P1]\t[SUB]\t[SUB]\t[P2]\t[BUFF]\t[DIST]\t[LEN R]\t[LEN Q]\t[CTX R]\t[CTX Q]\t[FRM]\t[TAGS]
    39\tG\t.\t38\t34\t38\t85\t84\tGTAGTAG\tGTA.TAG\t1\t1\tref\tquery
    73\tT\tA\t72\t13\t13\t85\t84\tGGATTGA\tGGAATGA\t1\t1\tref\tquery
    """
            )
        )

        actual = df.get_probes()
        expected = (
            str(
                Probe(
                    ProbeHeader(chrom="ref", pos=39, interval=ProbeInterval(3, 4)),
                    full_sequence="GTAGTAG",
                )
            )
            + "\n"
            + str(
                Probe(
                    ProbeHeader(chrom="ref", pos=73, interval=ProbeInterval(3, 4)),
                    full_sequence="GGATTGA",
                )
            ),
            str(
                Probe(
                    ProbeHeader(chrom="query", pos=38, interval=ProbeInterval(3, 3)),
                    full_sequence="GTATAG",
                )
            )
            + "\n"
            + str(
                Probe(
                    ProbeHeader(chrom="query", pos=72, interval=ProbeInterval(3, 4)),
                    full_sequence="GGAATGA",
                )
            ),
        )

        assert actual == expected

    def test_makeTruthPanelFromSnpsDataframe_probeNearGeneStartReturnsTruncatedLeftFlank(
        self
    ):
        df = ShowSnps.to_dataframe(
            StringIO(
                """/home/michael/Projects/pandora1_paper/tests/test_cases/ref.fa /home/michael/Projects/pandora1_paper/tests/test_cases/query.fa
    NUCMER

    [P1]\t[SUB]\t[SUB]\t[P2]\t[BUFF]\t[DIST]\t[LEN R]\t[LEN Q]\t[CTX R]\t[CTX Q]\t[FRM]\t[TAGS]
    2\tG\t.\t3\t34\t38\t85\t84\t--AGTAG\t-TA.TAG\t1\t1\tref\tquery
    """
            )
        )

        actual = df.get_probes()
        expected_ref = str(
            Probe(
                ProbeHeader(chrom="ref", pos=2, interval=ProbeInterval(1, 2)),
                full_sequence="AGTAG",
            )
        )
        expected_query = str(
            Probe(
                ProbeHeader(chrom="query", pos=3, interval=ProbeInterval(2, 2)),
                full_sequence="TATAG",
            )
        )
        expected = (str(expected_ref), str(expected_query))

        assert actual == expected

    def test_makeTruthPanelFromSnpsDataframe_probeNearGeneEndReturnsTruncatedRightFlank(
        self
    ):
        df = ShowSnps.to_dataframe(
            StringIO(
                """/home/michael/Projects/pandora1_paper/tests/test_cases/ref.fa /home/michael/Projects/pandora1_paper/tests/test_cases/query.fa
    NUCMER

    [P1]\t[SUB]\t[SUB]\t[P2]\t[BUFF]\t[DIST]\t[LEN R]\t[LEN Q]\t[CTX R]\t[CTX Q]\t[FRM]\t[TAGS]
    12\tG\t.\t13\t34\t38\t85\t84\tAAAGTA-\tATA.T--\t1\t1\tref\tquery
    """
            )
        )

        actual = df.get_probes()
        expected_ref = str(
            Probe(
                ProbeHeader(chrom="ref", pos=12, interval=ProbeInterval(3, 4)),
                full_sequence="AAAGTA",
            )
        )
        expected_query = str(
            Probe(
                ProbeHeader(chrom="query", pos=13, interval=ProbeInterval(3, 3)),
                full_sequence="ATAT",
            )
        )
        expected = (str(expected_ref), str(expected_query))

        assert actual == expected

    def test_makeTruthPanelFromSnpsDataframe_probeAtGeneStartReturnsTruncatedLeftFlank(
        self
    ):
        df = ShowSnps.to_dataframe(
            StringIO(
                """/home/michael/Projects/pandora1_paper/tests/test_cases/ref.fa /home/michael/Projects/pandora1_paper/tests/test_cases/query.fa
    NUCMER

    [P1]\t[SUB]\t[SUB]\t[P2]\t[BUFF]\t[DIST]\t[LEN R]\t[LEN Q]\t[CTX R]\t[CTX Q]\t[FRM]\t[TAGS]
    1\tG\t.\t1\t34\t38\t85\t84\t---GTAG\t---.TAG\t1\t1\tref\tquery
    """
            )
        )

        actual = df.get_probes()
        expected_ref = str(
            Probe(
                ProbeHeader(chrom="ref", pos=1, interval=ProbeInterval(0, 1)),
                full_sequence="GTAG",
            )
        )
        expected_query = str(
            Probe(
                ProbeHeader(chrom="query", pos=1, interval=ProbeInterval(0, 0)),
                full_sequence="TAG",
            )
        )
        expected = (str(expected_ref), str(expected_query))

        assert actual == expected

    def test_makeTruthPanelFromSnpsDataframe_probeAtGeneEndReturnsTruncatedRightFlank(
        self
    ):
        df = ShowSnps.to_dataframe(
            StringIO(
                """/home/michael/Projects/pandora1_paper/tests/test_cases/ref.fa /home/michael/Projects/pandora1_paper/tests/test_cases/query.fa
    NUCMER

    [P1]\t[SUB]\t[SUB]\t[P2]\t[BUFF]\t[DIST]\t[LEN R]\t[LEN Q]\t[CTX R]\t[CTX Q]\t[FRM]\t[TAGS]
    10\tG\t.\t10\t34\t38\t85\t84\tAAAG---\tAAA.---\t1\t1\tref\tquery
    """
            )
        )

        actual = df.get_probes()
        expected_ref = str(
            Probe(
                ProbeHeader(chrom="ref", pos=10, interval=ProbeInterval(3, 4)),
                full_sequence="AAAG",
            )
        )
        expected_query = str(
            Probe(
                ProbeHeader(chrom="query", pos=10, interval=ProbeInterval(3, 3)),
                full_sequence="AAA",
            )
        )
        expected = (str(expected_ref), str(expected_query))

        assert actual == expected

    def test_makeTruthPanelFromSnpsDataframe_mergeConsecutiveRecordsForIndelInQuery(
        self
    ):
        df = ShowSnps.to_dataframe(
            StringIO(
                """ref query
    NUCMER

    [P1]\t[SUB]\t[SUB]\t[P2]\t[BUFF]\t[DIST]\t[LEN R]\t[LEN Q]\t[CTX R]\t[CTX Q]\t[FRM]\t[TAGS]
    39\tG\t.\t38\t34\t38\t85\t84\tGTAGTAG\tGTA.TAG\t1\t1\tref\tquery
    73\tT\t.\t72\t13\t13\t85\t84\tGGATTTG\tGGA.TGA\t1\t1\tref\tquery
    74\tT\t.\t72\t13\t13\t85\t84\tGATTTGA\tGGA.TGA\t1\t1\tref\tquery
    75\tT\t.\t72\t13\t13\t85\t84\tATTTGAA\tGGA.TGA\t1\t1\tref\tquery
    79\tT\tA\t78\t13\t13\t85\t84\tGGATTGA\tGGAATGA\t1\t1\tref\tquery
    """
            )
        )

        actual = df.get_probes()
        expected_ref = str(
            Probe(
                ProbeHeader(chrom="ref", pos=39, interval=ProbeInterval(3, 4)),
                full_sequence="GTAGTAG",
            )
        )
        expected_ref += "\n" + str(
            Probe(
                ProbeHeader(chrom="ref", pos=73, interval=ProbeInterval(3, 6)),
                full_sequence="GGATTTGAA",
            )
        )
        expected_ref += "\n" + str(
            Probe(
                ProbeHeader(chrom="ref", pos=79, interval=ProbeInterval(3, 4)),
                full_sequence="GGATTGA",
            )
        )
        expected_query = str(
            Probe(
                ProbeHeader(chrom="query", pos=38, interval=ProbeInterval(3, 3)),
                full_sequence="GTATAG",
            )
        )
        expected_query += "\n" + str(
            Probe(
                ProbeHeader(chrom="query", pos=72, interval=ProbeInterval(3, 3)),
                full_sequence="GGATGA",
            )
        )
        expected_query += "\n" + str(
            Probe(
                ProbeHeader(chrom="query", pos=78, interval=ProbeInterval(3, 4)),
                full_sequence="GGAATGA",
            )
        )
        expected = (str(expected_ref), str(expected_query))

        assert actual == expected

    def test_makeTruthPanelFromSnpsDataframe_mergeConsecutiveRecordsForIndelInRef(self):
        df = ShowSnps.to_dataframe(
            StringIO(
                """ref query
    NUCMER

    [P1]\t[SUB]\t[SUB]\t[P2]\t[BUFF]\t[DIST]\t[LEN R]\t[LEN Q]\t[CTX R]\t[CTX Q]\t[FRM]\t[TAGS]
    39\tG\t.\t38\t34\t38\t85\t84\tGTAGTAG\tGTA.TAG\t1\t1\tref\tquery
    72\t.\tT\t73\t13\t13\t85\t84\tGGA.TGA\tGGATTTG\t1\t1\tref\tquery
    72\t.\tT\t74\t13\t13\t85\t84\tGGA.TGA\tGATTTGA\t1\t1\tref\tquery
    72\t.\tT\t75\t13\t13\t85\t84\tGGA.TGA\tATTTGAA\t1\t1\tref\tquery
    79\tT\tA\t78\t13\t13\t85\t84\tGGATTGA\tGGAATGA\t1\t1\tref\tquery
    """
            )
        )

        actual = df.get_probes()
        expected_ref = str(
            Probe(
                ProbeHeader(chrom="ref", pos=39, interval=ProbeInterval(3, 4)),
                full_sequence="GTAGTAG",
            )
        )
        expected_ref += "\n" + str(
            Probe(
                ProbeHeader(chrom="ref", pos=72, interval=ProbeInterval(3, 3)),
                full_sequence="GGATGA",
            )
        )
        expected_ref += "\n" + str(
            Probe(
                ProbeHeader(chrom="ref", pos=79, interval=ProbeInterval(3, 4)),
                full_sequence="GGATTGA",
            )
        )
        expected_query = str(
            Probe(
                ProbeHeader(chrom="query", pos=38, interval=ProbeInterval(3, 3)),
                full_sequence="GTATAG",
            )
        )
        expected_query += "\n" + str(
            Probe(
                ProbeHeader(chrom="query", pos=73, interval=ProbeInterval(3, 6)),
                full_sequence="GGATTTGAA",
            )
        )
        expected_query += "\n" + str(
            Probe(
                ProbeHeader(chrom="query", pos=78, interval=ProbeInterval(3, 4)),
                full_sequence="GGAATGA",
            )
        )
        expected = (str(expected_ref), str(expected_query))

        assert actual == expected

    def test_makeTruthPanelFromSnpsDataframe_mergeConsecutiveRecordsForMnpInRef(self):
        df = ShowSnps.to_dataframe(
            StringIO(
                """ref query
    NUCMER

    [P1]\t[SUB]\t[SUB]\t[P2]\t[BUFF]\t[DIST]\t[LEN R]\t[LEN Q]\t[CTX R]\t[CTX Q]\t[FRM]\t[TAGS]
    39\tG\t.\t38\t34\t38\t85\t84\tGTAGTAG\tGTA.TAG\t1\t1\tref\tquery
    72\tA\tT\t73\t13\t13\t85\t84\tGGAAGCA\tGGATTTG\t1\t1\tref\tquery
    73\tG\tT\t74\t13\t13\t85\t84\tGAAGCAA\tGATTTGA\t1\t1\tref\tquery
    74\tC\tT\t75\t13\t13\t85\t84\tAAGCAAA\tATTTGAA\t1\t1\tref\tquery
    79\tT\tA\t78\t13\t13\t85\t84\tGGATTGA\tGGAATGA\t1\t1\tref\tquery
    """
            )
        )

        actual = df.get_probes()
        expected_ref = str(
            Probe(
                ProbeHeader(chrom="ref", pos=39, interval=ProbeInterval(3, 4)),
                full_sequence="GTAGTAG",
            )
        )
        expected_ref += "\n" + str(
            Probe(
                ProbeHeader(chrom="ref", pos=72, interval=ProbeInterval(3, 6)),
                full_sequence="GGAAGCAAA",
            )
        )
        expected_ref += "\n" + str(
            Probe(
                ProbeHeader(chrom="ref", pos=79, interval=ProbeInterval(3, 4)),
                full_sequence="GGATTGA",
            )
        )
        expected_query = str(
            Probe(
                ProbeHeader(chrom="query", pos=38, interval=ProbeInterval(3, 3)),
                full_sequence="GTATAG",
            )
        )
        expected_query += "\n" + str(
            Probe(
                ProbeHeader(chrom="query", pos=73, interval=ProbeInterval(3, 6)),
                full_sequence="GGATTTGAA",
            )
        )
        expected_query += "\n" + str(
            Probe(
                ProbeHeader(chrom="query", pos=78, interval=ProbeInterval(3, 4)),
                full_sequence="GGAATGA",
            )
        )
        expected = (str(expected_ref), str(expected_query))

        assert actual == expected

    def test_makePosZeroBased_emptyReturnsEmpty(self):
        df = ShowSNPsDataframe()

        actual = df.make_pos_zero_based()
        expected = df

        assert actual.equals(expected)

    def test_makePosZeroBased_oneValueReturnsValueMinusOne(self):
        df = ShowSNPsDataframe({"ref_pos": [1], "query_pos": [8]})

        actual = df.make_pos_zero_based()
        expected = ShowSNPsDataframe({"ref_pos": [0], "query_pos": [7]})

        assert actual.equals(expected)

    def test_makePosZeroBased_twoValuesReturnsAllValuesMinusOne(self):
        df = ShowSNPsDataframe({"ref_pos": [1, 89], "query_pos": [8, 10]})

        actual = df.make_pos_zero_based()
        expected = ShowSNPsDataframe({"ref_pos": [0, 88], "query_pos": [7, 9]})

        assert actual.equals(expected)
