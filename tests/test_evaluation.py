from io import StringIO
from pathlib import Path

import pandas as pd
import pysam
import pytest

from evaluate.evaluation import (
    generate_mummer_snps,
    is_mapping_invalid,
    make_truth_panels,
)
from evaluate.mummer import NucmerError, ShowSnps

REF = Path("tests/test_cases/ref.fa")
QUERY = Path("tests/test_cases/query.fa")
SNPS = Path("tests/test_cases/out.snps")


def create_sam_header(name: str, length: int) -> pysam.AlignmentHeader:
    return pysam.AlignmentHeader.from_text(
        f"@SQ	SN:{name}	LN:{length}\n@PG	ID:bwa	PN:bwa	VN:0.7.17-r1188	CL:bwa mem -t 1 panel.fa -"
    )


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

    actual = generate_mummer_snps(
        reference, query, prefix, flank_width=context, indels=True
    )
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
    df = ShowSnps.to_dataframe(
        StringIO(
            """/home/michael/Projects/pandora1_paper/tests/test_cases/ref.fa /home/michael/Projects/pandora1_paper/tests/test_cases/query.fa
NUCMER

[P1]\t[SUB]\t[SUB]\t[P2]\t[BUFF]\t[DIST]\t[LEN R]\t[LEN Q]\t[CTX R]\t[CTX Q]\t[FRM]\t[TAGS]
"""
        )
    )

    actual = make_truth_panels(df)
    expected = ("", "")

    assert actual == expected


def test_makeTruthPanelFromSnpsDataframe_invalidDataframeRaisesError():
    df = pd.DataFrame(
        {"ref_pos": [39, 73], "ref_sub": ["G", "T"], "query_len": [84, 84]}
    )

    with pytest.raises(AttributeError):
        make_truth_panels(df)


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

    actual = make_truth_panels(df)
    expected = (
        ">ref_POS=39_SUB=G\nGTAGTAG\n>ref_POS=73_SUB=T\nGGATTGA\n",
        ">query_POS=38_SUB=.\nGTATAG\n>query_POS=72_SUB=A\nGGAATGA\n",
    )

    assert actual == expected


def test_makeTruthPanelFromSnpsDataframe_mergeConsecutiveRecordsForIndelInQuery():
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

    actual = make_truth_panels(df)
    expected_ref_probes = (
        ">ref_POS=39_SUB=G\nGTAGTAG\n"
        ">ref_POS=73_SUB=TTT\nGGATTTGAA\n"
        ">ref_POS=79_SUB=T\nGGATTGA\n"
    )
    expected_query_probes = (
        ">query_POS=38_SUB=.\nGTATAG\n"
        ">query_POS=72_SUB=...\nGGATGA\n"
        ">query_POS=78_SUB=A\nGGAATGA\n"
    )
    expected = (expected_ref_probes, expected_query_probes)

    assert actual == expected


def test_makeTruthPanelFromSnpsDataframe_mergeConsecutiveRecordsForIndelInRef():
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

    actual = make_truth_panels(df)
    expected_ref_probes = (
        ">ref_POS=39_SUB=G\nGTAGTAG\n"
        ">ref_POS=72_SUB=...\nGGATGA\n"
        ">ref_POS=79_SUB=T\nGGATTGA\n"
    )
    expected_query_probes = (
        ">query_POS=38_SUB=.\nGTATAG\n"
        ">query_POS=73_SUB=TTT\nGGATTTGAA\n"
        ">query_POS=78_SUB=A\nGGAATGA\n"
    )
    expected = (expected_ref_probes, expected_query_probes)

    assert actual == expected


def test_makeTruthPanelFromSnpsDataframe_mergeConsecutiveRecordsForMnpInRef():
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

    actual = make_truth_panels(df)
    expected_ref_probes = (
        ">ref_POS=39_SUB=G\nGTAGTAG\n"
        ">ref_POS=72_SUB=AGC\nGGAAGCAAA\n"
        ">ref_POS=79_SUB=T\nGGATTGA\n"
    )
    expected_query_probes = (
        ">query_POS=38_SUB=.\nGTATAG\n"
        ">query_POS=73_SUB=TTT\nGGATTTGAA\n"
        ">query_POS=78_SUB=A\nGGAATGA\n"
    )
    expected = (expected_ref_probes, expected_query_probes)

    assert actual == expected


def test_isMappingInvalid_unmappedEntry_returnTrue():
    header = create_sam_header("C15154T", 201)
    record = pysam.AlignedSegment.fromstring(
        "GC00004785_pos200_entry0	4	C15154T	124	48	69M	*	0	0	CAAATCGGAAGCTAACAGAGCCAATACGCGCCTTGACGCCCAGGACTATTTTGATTGCCTGCGCTGCTT	*	NM:i:0	MD:Z:69	AS:i:69	XS:i:53",
        header,
    )

    assert is_mapping_invalid(record)


def test_isMappingInvalid_mappedEntry_returnFalse():
    header = create_sam_header("C15154T", 201)
    record = pysam.AlignedSegment.fromstring(
        "GC00004785_pos200_entry0	0	C15154T	124	48	69M	*	0	0	CAAATCGGAAGCTAACAGAGCCAATACGCGCCTTGACGCCCAGGACTATTTTGATTGCCTGCGCTGCTT	*	NM:i:0	MD:Z:69	AS:i:69	XS:i:53",
        header,
    )
    is_mapping_valid = not is_mapping_invalid(record)

    assert is_mapping_valid


def test_isMappingInvalid_supplementaryEntry_returnTrue():
    header = create_sam_header("C15154T", 201)
    record = pysam.AlignedSegment.fromstring(
        "GC00004785_pos200_entry0	2048	C15154T	124	48	69M	*	0	0	CAAATCGGAAGCTAACAGAGCCAATACGCGCCTTGACGCCCAGGACTATTTTGATTGCCTGCGCTGCTT	*	NM:i:0	MD:Z:69	AS:i:69	XS:i:53",
        header,
    )

    assert is_mapping_invalid(record)


def test_isMappingInvalid_secondaryEntry_returnTrue():
    header = create_sam_header("C15154T", 201)
    record = pysam.AlignedSegment.fromstring(
        "GC00004785_pos200_entry0	256	C15154T	124	48	69M	*	0	0	CAAATCGGAAGCTAACAGAGCCAATACGCGCCTTGACGCCCAGGACTATTTTGATTGCCTGCGCTGCTT	*	NM:i:0	MD:Z:69	AS:i:69	XS:i:53",
        header,
    )

    assert is_mapping_invalid(record)
