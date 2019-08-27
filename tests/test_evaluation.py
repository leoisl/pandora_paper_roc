import pytest

from evaluate.evaluation import *
from evaluate.mummer import NucmerError
from evaluate.probe import Interval

TEST_DIR = Path("tests/test_cases")
REF = TEST_DIR / "ref.fa"
QUERY = TEST_DIR / "query.fa"
SNPS = TEST_DIR / "out.snps"


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
