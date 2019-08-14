from pathlib import Path
from evaluate.mummer import NucmerError, ShowSnps
from evaluate.evaluate import (
    generate_mummer_snps,
    make_truth_panels_from_snps_dataframe,
)
from io import StringIO
import pandas as pd
import pytest

REF = Path("tests/test_cases/ref.fa")
QUERY = Path("tests/test_cases/query.fa")
SNPS = Path("tests/test_cases/out.snps")


def test_getMummerSnps_invalidQueryFileRaisesNucmerError():
    reference = REF
    query = Path("foo")
    with pytest.raises(NucmerError):
        generate_mummer_snps(reference, query)


def test_getMummerSnps_validSequenceFilesProducesExpectedSnpsFile():
    reference = REF
    query = QUERY
    prefix = Path("/tmp/ref_query")
    context = 3

    actual = generate_mummer_snps(reference, query, prefix, flank_width=context)
    expected = StringIO(
        """/home/michael/Projects/pandora1_paper/tests/test_cases/ref.fa /home/michael/Projects/pandora1_paper/tests/test_cases/query.fa
NUCMER

[P1]\t[SUB]\t[SUB]\t[P2]\t[BUFF]\t[DIST]\t[LEN R]\t[LEN Q]\t[CTX R]\t[CTX Q]\t[FRM]\t[TAGS]
39\tG\t.\t38\t34\t38\t85\t84\tGTAGTAG\tGTA.TAG\t1\t1\tref\tquery
73\tT\tA\t72\t13\t13\t85\t84\tGGATTGA\tGGAATGA\t1\t1\tref\tquery
"""
    )

    assert actual.readlines() == expected.readlines()


def test_makeTruthPanelFromSnpsDataframe_emptyDataframeReturnsEmptyPanel():
    df = pd.DataFrame()

    actual = make_truth_panels_from_snps_dataframe(df)
    expected = ("", "")

    assert actual == expected


def test_makeTruthPanelFromSnpsDataframe_invalidDataframeRaisesError():
    df = pd.DataFrame(
        {"ref_pos": [39, 73], "ref_sub": ["G", "T"], "query_len": [84, 84]}
    )

    with pytest.raises(AttributeError):
        make_truth_panels_from_snps_dataframe(df)


def test_makeTruthPanelFromSnpsDataframe_validDataframeReturnsTwoProbesets():
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

    actual = make_truth_panels_from_snps_dataframe(df)
    expected = (
        ">ref_POS=39_SUB=G\nGTAGTAG\n>ref_POS=73_SUB=T\nGGATTGA\n",
        ">query_POS=38_SUB=.\nGTATAG\n>query_POS=72_SUB=A\nGGAATGA\n",
    )

    assert actual == expected
