from io import StringIO
from pathlib import Path

import pandas as pd
import pysam
import pytest

from evaluate.evaluation import *
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
        reference, query, prefix, flank_width=context, indels=True, print_header=False
    )
    expected = StringIO(
        """39\tG\t.\t38\t34\t38\t85\t84\tGTAGTAG\tGTA.TAG\t1\t1\tref\tquery
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
        str(
            Probe(
                ProbeHeader(chrom="ref", pos=39, interval=Interval(3, 4)),
                full_sequence="GTAGTAG",
            )
        )
        + "\n"
        + str(
            Probe(
                ProbeHeader(chrom="ref", pos=73, interval=Interval(3, 4)),
                full_sequence="GGATTGA",
            )
        ),
        str(
            Probe(
                ProbeHeader(chrom="query", pos=38, interval=Interval(3, 3)),
                full_sequence="GTATAG",
            )
        )
        + "\n"
        + str(
            Probe(
                ProbeHeader(chrom="query", pos=72, interval=Interval(3, 4)),
                full_sequence="GGAATGA",
            )
        ),
    )

    assert actual == expected


def test_makeTruthPanelFromSnpsDataframe_probeNearGeneStartReturnsTruncatedLeftFlank():
    df = ShowSnps.to_dataframe(
        StringIO(
            """/home/michael/Projects/pandora1_paper/tests/test_cases/ref.fa /home/michael/Projects/pandora1_paper/tests/test_cases/query.fa
NUCMER

[P1]\t[SUB]\t[SUB]\t[P2]\t[BUFF]\t[DIST]\t[LEN R]\t[LEN Q]\t[CTX R]\t[CTX Q]\t[FRM]\t[TAGS]
2\tG\t.\t3\t34\t38\t85\t84\t--AGTAG\t-TA.TAG\t1\t1\tref\tquery
"""
        )
    )

    actual = make_truth_panels(df)
    expected_ref = str(
        Probe(
            ProbeHeader(chrom="ref", pos=2, interval=Interval(1, 2)),
            full_sequence="AGTAG",
        )
    )
    expected_query = str(
        Probe(
            ProbeHeader(chrom="query", pos=3, interval=Interval(2, 2)),
            full_sequence="TATAG",
        )
    )
    expected = (str(expected_ref), str(expected_query))

    assert actual == expected


def test_makeTruthPanelFromSnpsDataframe_probeNearGeneEndReturnsTruncatedRightFlank():
    df = ShowSnps.to_dataframe(
        StringIO(
            """/home/michael/Projects/pandora1_paper/tests/test_cases/ref.fa /home/michael/Projects/pandora1_paper/tests/test_cases/query.fa
NUCMER

[P1]\t[SUB]\t[SUB]\t[P2]\t[BUFF]\t[DIST]\t[LEN R]\t[LEN Q]\t[CTX R]\t[CTX Q]\t[FRM]\t[TAGS]
12\tG\t.\t13\t34\t38\t85\t84\tAAAGTA-\tATA.T--\t1\t1\tref\tquery
"""
        )
    )

    actual = make_truth_panels(df)
    expected_ref = str(
        Probe(
            ProbeHeader(chrom="ref", pos=12, interval=Interval(3, 4)),
            full_sequence="AAAGTA",
        )
    )
    expected_query = str(
        Probe(
            ProbeHeader(chrom="query", pos=13, interval=Interval(3, 3)),
            full_sequence="ATAT",
        )
    )
    expected = (str(expected_ref), str(expected_query))

    assert actual == expected


def test_makeTruthPanelFromSnpsDataframe_probeAtGeneStartReturnsTruncatedLeftFlank():
    df = ShowSnps.to_dataframe(
        StringIO(
            """/home/michael/Projects/pandora1_paper/tests/test_cases/ref.fa /home/michael/Projects/pandora1_paper/tests/test_cases/query.fa
NUCMER

[P1]\t[SUB]\t[SUB]\t[P2]\t[BUFF]\t[DIST]\t[LEN R]\t[LEN Q]\t[CTX R]\t[CTX Q]\t[FRM]\t[TAGS]
1\tG\t.\t1\t34\t38\t85\t84\t---GTAG\t---.TAG\t1\t1\tref\tquery
"""
        )
    )

    actual = make_truth_panels(df)
    expected_ref = str(
        Probe(
            ProbeHeader(chrom="ref", pos=1, interval=Interval(0, 1)),
            full_sequence="GTAG",
        )
    )
    expected_query = str(
        Probe(
            ProbeHeader(chrom="query", pos=1, interval=Interval(0, 0)),
            full_sequence="TAG",
        )
    )
    expected = (str(expected_ref), str(expected_query))

    assert actual == expected


def test_makeTruthPanelFromSnpsDataframe_probeAtGeneEndReturnsTruncatedRightFlank():
    df = ShowSnps.to_dataframe(
        StringIO(
            """/home/michael/Projects/pandora1_paper/tests/test_cases/ref.fa /home/michael/Projects/pandora1_paper/tests/test_cases/query.fa
NUCMER

[P1]\t[SUB]\t[SUB]\t[P2]\t[BUFF]\t[DIST]\t[LEN R]\t[LEN Q]\t[CTX R]\t[CTX Q]\t[FRM]\t[TAGS]
10\tG\t.\t10\t34\t38\t85\t84\tAAAG---\tAAA.---\t1\t1\tref\tquery
"""
        )
    )

    actual = make_truth_panels(df)
    expected_ref = str(
        Probe(
            ProbeHeader(chrom="ref", pos=10, interval=Interval(3, 4)),
            full_sequence="AAAG",
        )
    )
    expected_query = str(
        Probe(
            ProbeHeader(chrom="query", pos=10, interval=Interval(3, 3)),
            full_sequence="AAA",
        )
    )
    expected = (str(expected_ref), str(expected_query))

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
    expected_ref = str(
        Probe(
            ProbeHeader(chrom="ref", pos=39, interval=Interval(3, 4)),
            full_sequence="GTAGTAG",
        )
    )
    expected_ref += "\n" + str(
        Probe(
            ProbeHeader(chrom="ref", pos=73, interval=Interval(3, 6)),
            full_sequence="GGATTTGAA",
        )
    )
    expected_ref += "\n" + str(
        Probe(
            ProbeHeader(chrom="ref", pos=79, interval=Interval(3, 4)),
            full_sequence="GGATTGA",
        )
    )
    expected_query = str(
        Probe(
            ProbeHeader(chrom="query", pos=38, interval=Interval(3, 3)),
            full_sequence="GTATAG",
        )
    )
    expected_query += "\n" + str(
        Probe(
            ProbeHeader(chrom="query", pos=72, interval=Interval(3, 3)),
            full_sequence="GGATGA",
        )
    )
    expected_query += "\n" + str(
        Probe(
            ProbeHeader(chrom="query", pos=78, interval=Interval(3, 4)),
            full_sequence="GGAATGA",
        )
    )
    expected = (str(expected_ref), str(expected_query))

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
    expected_ref = str(
        Probe(
            ProbeHeader(chrom="ref", pos=39, interval=Interval(3, 4)),
            full_sequence="GTAGTAG",
        )
    )
    expected_ref += "\n" + str(
        Probe(
            ProbeHeader(chrom="ref", pos=72, interval=Interval(3, 3)),
            full_sequence="GGATGA",
        )
    )
    expected_ref += "\n" + str(
        Probe(
            ProbeHeader(chrom="ref", pos=79, interval=Interval(3, 4)),
            full_sequence="GGATTGA",
        )
    )
    expected_query = str(
        Probe(
            ProbeHeader(chrom="query", pos=38, interval=Interval(3, 3)),
            full_sequence="GTATAG",
        )
    )
    expected_query += "\n" + str(
        Probe(
            ProbeHeader(chrom="query", pos=73, interval=Interval(3, 6)),
            full_sequence="GGATTTGAA",
        )
    )
    expected_query += "\n" + str(
        Probe(
            ProbeHeader(chrom="query", pos=78, interval=Interval(3, 4)),
            full_sequence="GGAATGA",
        )
    )
    expected = (str(expected_ref), str(expected_query))

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
    expected_ref = str(
        Probe(
            ProbeHeader(chrom="ref", pos=39, interval=Interval(3, 4)),
            full_sequence="GTAGTAG",
        )
    )
    expected_ref += "\n" + str(
        Probe(
            ProbeHeader(chrom="ref", pos=72, interval=Interval(3, 6)),
            full_sequence="GGAAGCAAA",
        )
    )
    expected_ref += "\n" + str(
        Probe(
            ProbeHeader(chrom="ref", pos=79, interval=Interval(3, 4)),
            full_sequence="GGATTGA",
        )
    )
    expected_query = str(
        Probe(
            ProbeHeader(chrom="query", pos=38, interval=Interval(3, 3)),
            full_sequence="GTATAG",
        )
    )
    expected_query += "\n" + str(
        Probe(
            ProbeHeader(chrom="query", pos=73, interval=Interval(3, 6)),
            full_sequence="GGATTTGAA",
        )
    )
    expected_query += "\n" + str(
        Probe(
            ProbeHeader(chrom="query", pos=78, interval=Interval(3, 4)),
            full_sequence="GGAATGA",
        )
    )
    expected = (str(expected_ref), str(expected_query))

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


def test_assessSamRecord_unmappedRecordReturnsUnmapped():
    header = create_sam_header(
        "GC00000422_2_SAMPLE=CFT073_POS=603_CALL_INTERVAL=[25,32)_SVTYPE=PH_SNPs_MEAN_FWD_COVG=23_MEAN_REV_COVG=13_GT_CONF=89.5987",
        57,
    )
    record = pysam.AlignedSegment.fromstring(
        ""
        "3_POS=14788_CALL_INTERVAL=[21,22)\t4\t*\t0\t0\t*\t*\t0\t0\tCGCGAAAGCCCTGACCATCTGCACCGTGTCTGACCACATCCGC\t*\tAS:i:0\tXS:i:0",
        header,
    )

    actual = assess_sam_record(record)
    expected = "unmapped"

    assert actual == expected


def test_assessSamRecord_incorrectSecondayAlignmentMismatchReturnsIncorrect():
    flag = 256
    cigar = "43M"
    nm = "NM:i:1"
    md = "MD:Z:21T21"
    mapq = 0
    pos = 5
    query_header = ProbeHeader(chrom="3", pos=14788, interval=Interval(21, 22))
    ref_header = ProbeHeader(
        chrom="GC00000422_2",
        sample="CFT073",
        pos=603,
        interval=Interval(25, 32),
        svtype="PH_SNPs",
        mean_fwd_covg=23,
        mean_rev_covg=13,
        gt_conf=89.5987,
    )
    sequence = "CGCGAAAGCCCTGACCATCTGCACCGTGTCTGACCACATCCGC"
    header = create_sam_header(str(ref_header), 57)
    record = pysam.AlignedSegment.fromstring(
        f"{query_header}\t{flag}\t{ref_header}\t{pos}\t{mapq}\t{cigar}\t*\t0\t0\t{sequence}\t*\t{nm}\t{md}\tAS:i:43\tXS:i:32",
        header,
    )

    actual = assess_sam_record(record)
    expected = "secondary_incorrect"

    assert actual == expected


def test_assessSamRecord_correctSecondayAlignmentReturnsCorrect():
    flag = 256
    cigar = "43M"
    nm = "NM:i:0"
    md = "MD:Z:43"
    mapq = 0
    pos = 5
    query_name = "3_POS=14788_CALL_INTERVAL=[21,22)"
    ref_name = "GC00000422_2_SAMPLE=CFT073_POS=603_CALL_INTERVAL=[25,32)_SVTYPE=PH_SNPs_MEAN_FWD_COVG=23_MEAN_REV_COVG=13_GT_CONF=89.5987"
    sequence = "CGCGAAAGCCCTGACCATCTGCACCGTGTCTGACCACATCCGC"
    header = create_sam_header(ref_name, 57)
    record = pysam.AlignedSegment.fromstring(
        f"{query_name}\t{flag}\t{ref_name}\t{pos}\t{mapq}\t{cigar}\t*\t0\t0\t{sequence}\t*\t{nm}\t{md}\tAS:i:43\tXS:i:32",
        header,
    )

    actual = assess_sam_record(record)
    expected = "secondary_correct"

    assert actual == expected


def test_assessSamRecord_incorrectPrimaryAlignmentMismatchReturnsIncorrect():
    header = create_sam_header(
        "GC00000422_2_SAMPLE=CFT073_POS=603_CALL_INTERVAL=[25,32)_SVTYPE=PH_SNPs_MEAN_FWD_COVG=23_MEAN_REV_COVG=13_GT_CONF=89.5987",
        57,
    )
    record = pysam.AlignedSegment.fromstring(
        "3_POS=14788_CALL_INTERVAL=[21,22)\t0\tGC00000422_2_SAMPLE=CFT073_POS=603_CALL_INTERVAL=[25,32)_SVTYPE=PH_SNPs_MEAN_FWD_COVG=23_MEAN_REV_COVG=13_GT_CONF=89.5987\t5\t60\t43M\t*\t0\t0\tCGCGAAAGCCCTGACCATCTGCACCGTGTCTGACCACATCCGC\t*\tNM:i:1\tMD:Z:21T21\tAS:i:43\tXS:i:32",
        header,
    )

    # query_probe_seq = TTGGCGCGAAAGCCCTGACCATCTGTACCGTGTCTGACCACATCCGCACTCACGAGC
    # truth_probe_seq =     CGCGAAAGCCCTGACCATCTGCACCGTGTCTGACCACATCCGC

    actual = assess_sam_record(record)
    expected = "incorrect"

    assert actual == expected


def test_assessSamRecord_correctPrimaryAlignmentMismatchInFlankReturnsCorrect():
    header = create_sam_header(
        "GC00000422_2_SAMPLE=CFT073_POS=603_CALL_INTERVAL=[25,32)_SVTYPE=PH_SNPs_MEAN_FWD_COVG=23_MEAN_REV_COVG=13_GT_CONF=89.5987",
        57,
    )
    record = pysam.AlignedSegment.fromstring(
        "3_POS=14788_CALL_INTERVAL=[21,22)\t0\tGC00000422_2_SAMPLE=CFT073_POS=603_CALL_INTERVAL=[25,32)_SVTYPE=PH_SNPs_MEAN_FWD_COVG=23_MEAN_REV_COVG=13_GT_CONF=89.5987\t5\t60\t43M\t*\t0\t0\tCGCGAAAGCCCTGACCATCTGCACCGTGTCTGACCACATCCGC\t*\tNM:i:1\tMD:Z:20T22\tAS:i:43\tXS:i:32",
        header,
    )

    # query_probe_seq = TTGGCGCGAAAGCCCTGACCATCTTCACCGTGTCTGACCACATCCGCACTCACGAGC
    # truth_probe_seq =     CGCGAAAGCCCTGACCATCTGCACCGTGTCTGACCACATCCGC

    actual = assess_sam_record(record)
    expected = "correct"

    assert actual == expected


def test_assessSamRecord_correctPrimaryAlignmentReturnsCorrect():
    header = create_sam_header(
        "GC00000422_2_SAMPLE=CFT073_POS=603_CALL_INTERVAL=[25,32)_SVTYPE=PH_SNPs_MEAN_FWD_COVG=23_MEAN_REV_COVG=13_GT_CONF=89.5987",
        57,
    )
    record = pysam.AlignedSegment.fromstring(
        "3_POS=14788_CALL_INTERVAL=[21,22)\t0\tGC00000422_2_SAMPLE=CFT073_POS=603_CALL_INTERVAL=[25,32)_SVTYPE=PH_SNPs_MEAN_FWD_COVG=23_MEAN_REV_COVG=13_GT_CONF=89.5987\t5\t60\t43M\t*\t0\t0\tCGCGAAAGCCCTGACCATCTGCACCGTGTCTGACCACATCCGC\t*\tNM:i:0\tMD:Z:43\tAS:i:43\tXS:i:32",
        header,
    )

    # query_probe_seq = TTGGCGCGAAAGCCCTGACCATCTGCACCGTGTCTGACCACATCCGCACTCACGAGC
    # truth_probe_seq =     CGCGAAAGCCCTGACCATCTGCACCGTGTCTGACCACATCCGC

    actual = assess_sam_record(record)
    expected = "correct"

    assert actual == expected
