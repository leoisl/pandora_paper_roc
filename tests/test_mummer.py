import pytest
import uuid
from pathlib import Path
from io import StringIO
import pandas as pd
import pandas.errors
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
        expected = """/home/michael/Projects/pandora1_paper/tests/test_cases/ref.fa /home/michael/Projects/pandora1_paper/tests/test_cases/query.fa
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
        expected = """/home/michael/Projects/pandora1_paper/tests/test_cases/ref.fa /home/michael/Projects/pandora1_paper/tests/test_cases/query.fa
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

        with pytest.raises(pandas.errors.ParserError):
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
        expected = pd.DataFrame(
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
                "ref_chrom": ["ref", "ref"],
                "query_chrom": ["query", "query"],
            }
        )

        assert actual.equals(expected)
