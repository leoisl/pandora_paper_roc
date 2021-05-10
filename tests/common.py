import tempfile
from typing import Type
import pysam
from evaluate.classifier import Classifier
from evaluate.probe import ProbeHeader, ProbeInterval
import tempfile
from pathlib import Path
from evaluate.classification import AlignmentAssessment
import pandas as pd

def create_classifier_with_two_entries(cls: Type) -> Type[Classifier]:
    flag = 0
    cigar = "56M"
    nm = "NM:i:0"
    md = "MD:Z:56"
    mapq = 60
    pos = 1
    query_header = ProbeHeader(interval=ProbeInterval(12, 17))
    sequence = "AAAAAAAAAAACGGCTCGCATAGACACGACGACGACACGTACGATCGATCAGTCAT"
    ref_header = ProbeHeader(
        chrom="GC00000422_2",
        sample="CFT073",
        pos=603,
        interval=ProbeInterval(25, 32),
        svtype="PH_SNPs",
        gt_conf=89.5987,
    )
    header = create_sam_header(str(ref_header), 64)
    contents = str(header) + "\n"
    record1 = pysam.AlignedSegment.fromstring(
        f"{query_header}\t{flag}\t{ref_header}\t{pos}\t{mapq}\t{cigar}\t*\t0\t0\t{sequence}\t*\t{nm}\t{md}\tAS:i:43\tXS:i:32",
        header,
    )
    contents += record1.to_string() + "\n"

    flag = 2048
    cigar = "43M"
    nm = "NM:i:1"
    md = "MD:Z:21T21"
    mapq = 0
    pos = 5
    query_header = ProbeHeader(chrom="3", pos=14788, interval=ProbeInterval(21, 22))
    sequence = "CGCGAAAGCCCTGACCATCTGCACCGTGTCTGACCACATCCGC"
    header = create_sam_header(str(ref_header), 57)
    record2 = pysam.AlignedSegment.fromstring(
        f"{query_header}\t{flag}\t{ref_header}\t{pos}\t{mapq}\t{cigar}\t*\t0\t0\t{sequence}\t*\t{nm}\t{md}\tAS:i:43\tXS:i:32",
        header,
    )
    contents += record2.to_string() + "\n"
    sam = create_tmp_sam(contents)
    return cls(sam)


def create_tmp_sam(contents: str) -> pysam.AlignmentFile:
    with tempfile.NamedTemporaryFile(mode="r+") as tmp:
        tmp.write(contents)
        tmp.truncate()
        return pysam.AlignmentFile(tmp.name)


def create_unmapped_sam_record() -> pysam.AlignedSegment:
    header = create_sam_header(
        "GC00000422_2_SAMPLE=CFT073_POS=603_CALL_INTERVAL=[25,32)_SVTYPE=PH_SNPs_MEAN_FWD_COVG=23_MEAN_REV_COVG=13_GT_CONF=89.5987",
        57,
    )
    record = pysam.AlignedSegment.fromstring(
        ""
        "3_POS=14788_CALL_INTERVAL=[21,22)\t4\t*\t0\t0\t*\t*\t0\t0\tCGCGAAAGCCCTGACCATCTGCACCGTGTCTGACCACATCCGC\t*\tAS:i:0\tXS:i:0",
        header,
    )
    return record


def create_partially_mapped_sam_record() -> pysam.AlignedSegment:
    ref_name = "reference"
    ref_length = 59
    header = create_sam_header(ref_name, ref_length)
    flag = 0
    cigar = "30M38S"
    nm = "NM:i:0"
    md = "MD:Z:30"
    mapq = 60
    pos = 6
    query_name = "IV=[23,33);"
    sequence = "AAAAAAAAAAAAAAAAAAAAAAACGGCTCGCATAGACACGACGACGACACGTACGATCGATCAGTCAT"
    sam_string = f"{query_name}\t{flag}\t{ref_name}\t{pos}\t{mapq}\t{cigar}\t*\t0\t0\t{sequence}\t*\t{nm}\t{md}\tAS:i:0\tXS:i:0"
    record = pysam.AlignedSegment.fromstring(sam_string, header)
    return record


def create_incorrect_secondary_sam_record() -> pysam.AlignedSegment:
    flag = 256
    cigar = "43M"
    nm = "NM:i:1"
    md = "MD:Z:21T21"
    mapq = 0
    pos = 5
    query_header = ProbeHeader(chrom="3", pos=14788, interval=ProbeInterval(21, 22))
    ref_header = ProbeHeader(
        chrom="GC00000422_2",
        sample="CFT073",
        pos=603,
        interval=ProbeInterval(25, 32),
        svtype="PH_SNPs",
        gt_conf=89.5987,
    )
    sequence = "CGCGAAAGCCCTGACCATCTGCACCGTGTCTGACCACATCCGC"
    header = create_sam_header(str(ref_header), 57)
    record = pysam.AlignedSegment.fromstring(
        f"{query_header}\t{flag}\t{ref_header}\t{pos}\t{mapq}\t{cigar}\t*\t0\t0\t{sequence}\t*\t{nm}\t{md}\tAS:i:43\tXS:i:32",
        header,
    )
    return record


def create_correct_secondary_sam_record() -> pysam.AlignedSegment:
    flag = 256
    cigar = "43M"
    nm = "NM:i:1"
    md = "MD:Z:19T23"
    mapq = 0
    pos = 5
    query_header = ProbeHeader(chrom="3", pos=14788, interval=ProbeInterval(21, 22))
    ref_header = ProbeHeader(
        chrom="GC00000422_2",
        sample="CFT073",
        pos=603,
        interval=ProbeInterval(25, 32),
        svtype="PH_SNPs",
        gt_conf=89.5987,
    )
    sequence = "CGCGAAAGCCCTGACCATCTGCACCGTGTCTGACCACATCCGC"
    header = create_sam_header(str(ref_header), 57)
    record = pysam.AlignedSegment.fromstring(
        f"{query_header}\t{flag}\t{ref_header}\t{pos}\t{mapq}\t{cigar}\t*\t0\t0\t{sequence}\t*\t{nm}\t{md}\tAS:i:43\tXS:i:32",
        header,
    )
    return record


def create_incorrect_supplementary_sam_record() -> pysam.AlignedSegment:
    flag = 2048
    cigar = "43M"
    nm = "NM:i:1"
    md = "MD:Z:21T21"
    mapq = 0
    pos = 5
    query_header = ProbeHeader(chrom="3", pos=14788, interval=ProbeInterval(21, 22))
    ref_header = ProbeHeader(
        chrom="GC00000422_2",
        sample="CFT073",
        pos=603,
        interval=ProbeInterval(25, 32),
        svtype="PH_SNPs",
        gt_conf=89.5987,
    )
    sequence = "CGCGAAAGCCCTGACCATCTGCACCGTGTCTGACCACATCCGC"
    header = create_sam_header(str(ref_header), 57)
    record = pysam.AlignedSegment.fromstring(
        f"{query_header}\t{flag}\t{ref_header}\t{pos}\t{mapq}\t{cigar}\t*\t0\t0\t{sequence}\t*\t{nm}\t{md}\tAS:i:43\tXS:i:32",
        header,
    )
    return record


def create_correct_supplementary_sam_record() -> pysam.AlignedSegment:
    flag = 2048
    cigar = "43M"
    nm = "NM:i:1"
    md = "MD:Z:19T23"
    mapq = 0
    pos = 5
    query_header = ProbeHeader(chrom="3", pos=14788, interval=ProbeInterval(21, 22))
    ref_header = ProbeHeader(
        chrom="GC00000422_2",
        sample="CFT073",
        pos=603,
        interval=ProbeInterval(25, 32),
        svtype="PH_SNPs",
        gt_conf=89.5987,
    )
    sequence = "CGCGAAAGCCCTGACCATCTGCACCGTGTCTGACCACATCCGC"
    header = create_sam_header(str(ref_header), 57)
    record = pysam.AlignedSegment.fromstring(
        f"{query_header}\t{flag}\t{ref_header}\t{pos}\t{mapq}\t{cigar}\t*\t0\t0\t{sequence}\t*\t{nm}\t{md}\tAS:i:43\tXS:i:32",
        header,
    )
    return record


def create_correct_primary_sam_record() -> pysam.AlignedSegment:
    flag = 0
    cigar = "56M"
    nm = "NM:i:0"
    md = "MD:Z:56"
    mapq = 60
    pos = 1
    query_header = ProbeHeader(interval=ProbeInterval(12, 17))
    ref_header = ProbeHeader(
        chrom="GC00000422_2",
        sample="CFT073",
        pos=603,
        interval=ProbeInterval(25, 32),
        svtype="PH_SNPs",
        gt_conf=89.5987,
    )
    sequence = "AAAAAAAAAAACGGCTCGCATAGACACGACGACGACACGTACGATCGATCAGTCAT"
    header = create_sam_header(str(ref_header), 64)
    record = pysam.AlignedSegment.fromstring(
        f"{query_header}\t{flag}\t{ref_header}\t{pos}\t{mapq}\t{cigar}\t*\t0\t0\t{sequence}\t*\t{nm}\t{md}\tAS:i:43\tXS:i:32",
        header,
    )
    return record


def create_incorrect_primary_sam_record() -> pysam.AlignedSegment:
    flag = 0
    cigar = "56M"
    nm = "NM:i:1"
    md = "MD:Z:12T43"
    mapq = 60
    pos = 1
    query_header = ProbeHeader(interval=ProbeInterval(12, 13))
    ref_header = ProbeHeader(
        chrom="GC00000422_2",
        sample="CFT073",
        pos=603,
        interval=ProbeInterval(25, 32),
        svtype="PH_SNPs",
        gt_conf=89.5987,
    )
    sequence = "AAAAAAAAAAACGGCTCGCATAGACACGACGACGACACGTACGATCGATCAGTCAT"
    header = create_sam_header(str(ref_header), 64)
    record = pysam.AlignedSegment.fromstring(
        f"{query_header}\t{flag}\t{ref_header}\t{pos}\t{mapq}\t{cigar}\t*\t0\t0\t{sequence}\t*\t{nm}\t{md}\tAS:i:43\tXS:i:32",
        header,
    )
    return record


def create_sam_header(name: str, length: int) -> pysam.AlignmentHeader:
    return pysam.AlignmentHeader.from_text(
        f"@SQ	SN:{name}	LN:{length}\n@PG	ID:bwa	PN:bwa	VN:0.7.17-r1188	CL:bwa mem -t 1 panel.fa -"
    )


TEST_CASES = Path("tests/test_cases")
TEST_VCF = TEST_CASES / "test.vcf"
TEST_PANEL = TEST_CASES / "test_panel.fa"
TEST_REF_SEQ = TEST_CASES / "test_reference.fa"
TEST_TMP_PANEL = "/tmp/deleteme.fa"
TEST_MAKE_PROBE_VCF = TEST_CASES / "test_make_probe.vcf"
TEST_QUERY_VCF = TEST_CASES / "test_query.vcf"
TEST_QUERY_REF = TEST_CASES / "test_query.fa"


def retrieve_entry_from_test_vcf(idx: int) -> pysam.VariantRecord:
    with pysam.VariantFile(TEST_VCF) as vcf:
        for i, record in enumerate(vcf):
            if i == idx:
                return record
    raise IndexError("You asked for an index that is beyond the number in the test VCF")


def retrieve_entry_from_test_query_vcf(idx: int) -> pysam.VariantRecord:
    with pysam.VariantFile(TEST_QUERY_VCF) as vcf:
        for i, record in enumerate(vcf):
            if i == idx:
                return record
    raise IndexError("You asked for an index that is beyond the number in the test VCF")


def create_tmp_file(contents: str) -> Path:
    with tempfile.NamedTemporaryFile(mode="r+", delete=False) as tmp:
        tmp.write(contents)
        tmp.truncate()

    return Path(tmp.name)


def create_recall_report_row(
    truth_probe_header:str, classification: AlignmentAssessment, gt_conf: float = 0, sample: str = "sample1", with_gt_conf=False
) -> pd.Series:
    vcf_probe_header = ProbeHeader(gt_conf=gt_conf)
    data = {
        "sample": sample,
        "query_probe_header": str(truth_probe_header),
        "ref_probe_header": str(vcf_probe_header),
        "classification": classification.value,
        "good_eval": classification.value in ["primary_correct", "secondary_correct", "supplementary_correct"],
        "PVID": None,
        "NB_ALL": None,
        "ALL_ID": None,
        "NB_DIFF_ALL_SEQ": None,
        "ALL_SEQ_ID": None,
        "NB_OF_SAMPLES": None,
    }
    if with_gt_conf:
        data["GT_CONF"] = gt_conf

    return pd.Series(data=data)


def create_precision_report_row(
    classification: float, gt_conf: float = 0, sample: str = "sample1"
) -> pd.Series:
    ref_probe_header = ProbeHeader()
    pandora_probe_header = ProbeHeader(gt_conf=gt_conf)
    data = {
        "sample": sample,
        "query_probe_header": str(pandora_probe_header),
        "ref_probe_header": str(ref_probe_header),
        "classification": classification,
    }
    return pd.Series(data=data)
