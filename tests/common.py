import tempfile
from typing import Type
import pysam

from evaluate.classifier import Classifier, RecallClassifier
from evaluate.probe import ProbeHeader, ProbeInterval
from tests.test_evaluation import create_sam_header


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
        mean_fwd_covg=23,
        mean_rev_covg=13,
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
    query_name = "INTERVAL=[23,33);"
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
        mean_fwd_covg=23,
        mean_rev_covg=13,
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
        mean_fwd_covg=23,
        mean_rev_covg=13,
        gt_conf=89.5987,
    )
    sequence = "AAAAAAAAAAACGGCTCGCATAGACACGACGACGACACGTACGATCGATCAGTCAT"
    header = create_sam_header(str(ref_header), 64)
    record = pysam.AlignedSegment.fromstring(
        f"{query_header}\t{flag}\t{ref_header}\t{pos}\t{mapq}\t{cigar}\t*\t0\t0\t{sequence}\t*\t{nm}\t{md}\tAS:i:43\tXS:i:32",
        header,
    )
    return record
