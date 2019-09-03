import tempfile
from io import StringIO

import pandas as pd
import pysam

from evaluate.probe import Interval, ProbeHeader
from evaluate.recall import (
    RecallReporter,
    RecallClassifier,
    RecallClassification,
    RecallCalculator,
    StatisticalClassification,
)
from .test_evaluation import create_sam_header


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
    return record


def create_correct_secondary_sam_record() -> pysam.AlignedSegment:
    flag = 256
    cigar = "43M"
    nm = "NM:i:1"
    md = "MD:Z:19T23"
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
    return record


def create_incorrect_supplementary_sam_record() -> pysam.AlignedSegment:
    flag = 2048
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
    return record


def create_correct_supplementary_sam_record() -> pysam.AlignedSegment:
    flag = 2048
    cigar = "43M"
    nm = "NM:i:1"
    md = "MD:Z:19T23"
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
    return record


def create_correct_primary_sam_record() -> pysam.AlignedSegment:
    flag = 0
    cigar = "56M"
    nm = "NM:i:0"
    md = "MD:Z:56"
    mapq = 60
    pos = 1
    query_header = ProbeHeader(interval=Interval(12, 17))
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
    query_header = ProbeHeader(interval=Interval(12, 13))
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
    sequence = "AAAAAAAAAAACGGCTCGCATAGACACGACGACGACACGTACGATCGATCAGTCAT"
    header = create_sam_header(str(ref_header), 64)
    record = pysam.AlignedSegment.fromstring(
        f"{query_header}\t{flag}\t{ref_header}\t{pos}\t{mapq}\t{cigar}\t*\t0\t0\t{sequence}\t*\t{nm}\t{md}\tAS:i:43\tXS:i:32",
        header,
    )
    return record


def create_recall_classifier_with_two_entries() -> RecallClassifier:
    flag = 0
    cigar = "56M"
    nm = "NM:i:0"
    md = "MD:Z:56"
    mapq = 60
    pos = 1
    query_header = ProbeHeader(interval=Interval(12, 17))
    sequence = "AAAAAAAAAAACGGCTCGCATAGACACGACGACGACACGTACGATCGATCAGTCAT"
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
    query_header = ProbeHeader(chrom="3", pos=14788, interval=Interval(21, 22))
    sequence = "CGCGAAAGCCCTGACCATCTGCACCGTGTCTGACCACATCCGC"
    header = create_sam_header(str(ref_header), 57)
    record2 = pysam.AlignedSegment.fromstring(
        f"{query_header}\t{flag}\t{ref_header}\t{pos}\t{mapq}\t{cigar}\t*\t0\t0\t{sequence}\t*\t{nm}\t{md}\tAS:i:43\tXS:i:32",
        header,
    )
    contents += record2.to_string() + "\n"
    sam = create_tmp_sam(contents)
    return RecallClassifier(sam)


class TestRecallClassification:
    def test_equality_twoEqualReturnsTrue(self):
        c1 = RecallClassification()
        c1.truth_probe = Probe(ProbeHeader(chrom="2", pos=5))
        c2 = RecallClassification()
        c2.truth_probe = Probe(ProbeHeader(chrom="2", pos=5))

        assert c1 == c2

    def test_equality_twoNonEqualReturnsFalse(self):
        c1 = RecallClassification()
        c1.truth_probe = Probe(ProbeHeader(chrom="2", pos=5))
        c2 = RecallClassification()
        c2.truth_probe = Probe(ProbeHeader(chrom="3", pos=5))

        assert c1 != c2

    def test_wholeProbeMaps_unmappedRecordReturnsFalse(self):
        ref_name = "reference"
        ref_length = 55
        header = create_sam_header(ref_name, ref_length)
        flag = 4
        cigar = "*"
        nm = ""
        md = ""
        mapq = 0
        pos = 0
        query_name = "INTERVAL=[11,21);"
        ref_name = "*"
        sequence = "AAAAAAAAAAACGGCTCGCATAGACACGACGACGACACGTACGATCGATCAGTCAT"
        sam_string = f"{query_name}\t{flag}\t{ref_name}\t{pos}\t{mapq}\t{cigar}\t*\t0\t0\t{sequence}\t*\tAS:i:0\tXS:i:0"
        record = pysam.AlignedSegment.fromstring(sam_string, header)
        classification = RecallClassification(record=record)

        assert not classification._whole_probe_maps()

    def test_wholeProbeMaps_probeCompletelyMapsReturnsTrue(self):
        ref_name = "reference"
        ref_length = 55
        header = create_sam_header(ref_name, ref_length)
        flag = 0
        cigar = "9S47M"
        nm = "NM:i:0"
        md = "MD:Z:47"
        mapq = 60
        pos = 1
        query_name = "INTERVAL=[11,21);"
        sequence = "AAAAAAAAAAACGGCTCGCATAGACACGACGACGACACGTACGATCGATCAGTCAT"
        sam_string = f"{query_name}\t{flag}\t{ref_name}\t{pos}\t{mapq}\t{cigar}\t*\t0\t0\t{sequence}\t*\t{nm}\t{md}\tAS:i:0\tXS:i:0"
        record = pysam.AlignedSegment.fromstring(sam_string, header)
        classification = RecallClassification(record=record)

        assert classification._whole_probe_maps()

    def test_wholeProbeMaps_probeStartsAtFirstAlignmentPositionMapsReturnsTrue(self):
        ref_name = "reference"
        ref_length = 55
        header = create_sam_header(ref_name, ref_length)
        flag = 0
        cigar = "11S45M"
        nm = "NM:i:0"
        md = "MD:Z:45"
        mapq = 60
        pos = 1
        query_name = "INTERVAL=[11,21);"
        sequence = "AAAAAAAAAAACGGCTCGCATAGACACGACGACGACACGTACGATCGATCAGTCAT"
        sam_string = f"{query_name}\t{flag}\t{ref_name}\t{pos}\t{mapq}\t{cigar}\t*\t0\t0\t{sequence}\t*\t{nm}\t{md}\tAS:i:0\tXS:i:0"
        record = pysam.AlignedSegment.fromstring(sam_string, header)
        classification = RecallClassification(record=record)

        assert classification._whole_probe_maps()

    def test_wholeProbeMaps_probeStartOneBaseBeforeFirstAlignmentPositionMapsReturnsFalse(
        self
    ):
        ref_name = "reference"
        ref_length = 55
        header = create_sam_header(ref_name, ref_length)
        flag = 0
        cigar = "12S44M"
        nm = "NM:i:0"
        md = "MD:Z:44"
        mapq = 60
        pos = 1
        query_name = "INTERVAL=[11,21);"
        sequence = "AAAAAAAAAAACGGCTCGCATAGACACGACGACGACACGTACGATCGATCAGTCAT"
        sam_string = f"{query_name}\t{flag}\t{ref_name}\t{pos}\t{mapq}\t{cigar}\t*\t0\t0\t{sequence}\t*\t{nm}\t{md}\tAS:i:0\tXS:i:0"
        record = pysam.AlignedSegment.fromstring(sam_string, header)
        classification = RecallClassification(record=record)

        assert not classification._whole_probe_maps()

    def test_wholeProbeMaps_probeStartsThreeBaseBeforeFirstAlignmentPositionMapsReturnsFalse(
        self
    ):
        ref_name = "reference"
        ref_length = 55
        header = create_sam_header(ref_name, ref_length)
        flag = 0
        cigar = "14S42M"
        nm = "NM:i:0"
        md = "MD:Z:42"
        mapq = 60
        pos = 1
        query_name = "INTERVAL=[11,21);"
        sequence = "AAAAAAAAAAACGGCTCGCATAGACACGACGACGACACGTACGATCGATCAGTCAT"
        sam_string = f"{query_name}\t{flag}\t{ref_name}\t{pos}\t{mapq}\t{cigar}\t*\t0\t0\t{sequence}\t*\t{nm}\t{md}\tAS:i:0\tXS:i:0"
        record = pysam.AlignedSegment.fromstring(sam_string, header)
        classification = RecallClassification(record=record)

        assert not classification._whole_probe_maps()

    def test_wholeProbeMaps_ProbeEndsOneBaseBeforeAlignmentStartsReturnsFalse(self):
        ref_name = "reference"
        ref_length = 55
        header = create_sam_header(ref_name, ref_length)
        flag = 0
        cigar = "21S35M"
        nm = "NM:i:0"
        md = "MD:Z:35"
        mapq = 60
        pos = 1
        query_name = "INTERVAL=[11,21);"
        sequence = "AAAAAAAAAAACGGCTCGCATAGACACGACGACGACACGTACGATCGATCAGTCAT"
        sam_string = f"{query_name}\t{flag}\t{ref_name}\t{pos}\t{mapq}\t{cigar}\t*\t0\t0\t{sequence}\t*\t{nm}\t{md}\tAS:i:0\tXS:i:0"
        record = pysam.AlignedSegment.fromstring(sam_string, header)
        classification = RecallClassification(record=record)

        assert not classification._whole_probe_maps()

    def test_wholeProbeMaps_probeEndsThreeBasesBeforeAlignmentStartsReturnsFalse(self):
        ref_name = "reference"
        ref_length = 55
        header = create_sam_header(ref_name, ref_length)
        flag = 0
        cigar = "23S33M"
        nm = "NM:i:0"
        md = "MD:Z:33"
        mapq = 60
        pos = 1
        query_name = "INTERVAL=[11,21);"
        sequence = "AAAAAAAAAAACGGCTCGCATAGACACGACGACGACACGTACGATCGATCAGTCAT"
        sam_string = f"{query_name}\t{flag}\t{ref_name}\t{pos}\t{mapq}\t{cigar}\t*\t0\t0\t{sequence}\t*\t{nm}\t{md}\tAS:i:0\tXS:i:0"
        record = pysam.AlignedSegment.fromstring(sam_string, header)
        classification = RecallClassification(record=record)

        assert not classification._whole_probe_maps()

    def test_wholeProbeMaps_probeEndsOnLastBaseOfAlignmentReturnsTrue(self):
        ref_name = "reference"
        ref_length = 59
        header = create_sam_header(ref_name, ref_length)
        flag = 0
        cigar = "33M35S"
        nm = "NM:i:0"
        md = "MD:Z:33"
        mapq = 60
        pos = 6
        query_name = "INTERVAL=[23,33);"
        sequence = (
            "AAAAAAAAAAAAAAAAAAAAAAACGGCTCGCATAGACACGACGACGACACGTACGATCGATCAGTCAT"
        )
        sam_string = f"{query_name}\t{flag}\t{ref_name}\t{pos}\t{mapq}\t{cigar}\t*\t0\t0\t{sequence}\t*\t{nm}\t{md}\tAS:i:0\tXS:i:0"
        record = pysam.AlignedSegment.fromstring(sam_string, header)
        classification = RecallClassification(record=record)

        assert classification._whole_probe_maps()

    def test_wholeProbeMaps_probeEndsOneBaseAfterLastBaseOfAlignmentReturnsFalse(self):
        ref_name = "reference"
        ref_length = 59
        header = create_sam_header(ref_name, ref_length)
        flag = 0
        cigar = "32M36S"
        nm = "NM:i:0"
        md = "MD:Z:32"
        mapq = 60
        pos = 6
        query_name = "INTERVAL=[23,33);"
        sequence = (
            "AAAAAAAAAAAAAAAAAAAAAAACGGCTCGCATAGACACGACGACGACACGTACGATCGATCAGTCAT"
        )
        sam_string = f"{query_name}\t{flag}\t{ref_name}\t{pos}\t{mapq}\t{cigar}\t*\t0\t0\t{sequence}\t*\t{nm}\t{md}\tAS:i:0\tXS:i:0"
        record = pysam.AlignedSegment.fromstring(sam_string, header)
        classification = RecallClassification(record=record)

        assert not classification._whole_probe_maps()

    def test_wholeProbeMaps_probeEndsThreeBasesAfterLastBaseOfAlignmentReturnsFalse(
        self
    ):
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
        sequence = (
            "AAAAAAAAAAAAAAAAAAAAAAACGGCTCGCATAGACACGACGACGACACGTACGATCGATCAGTCAT"
        )
        sam_string = f"{query_name}\t{flag}\t{ref_name}\t{pos}\t{mapq}\t{cigar}\t*\t0\t0\t{sequence}\t*\t{nm}\t{md}\tAS:i:0\tXS:i:0"
        record = pysam.AlignedSegment.fromstring(sam_string, header)
        classification = RecallClassification(record=record)

        assert not classification._whole_probe_maps()

    def test_wholeProbeMaps_probeStartsAfterLastBaseOfAlignmentReturnsFalse(self):
        ref_name = "reference"
        ref_length = 59
        header = create_sam_header(ref_name, ref_length)
        flag = 0
        cigar = "22M46S"
        nm = "NM:i:0"
        md = "MD:Z:22"
        mapq = 60
        pos = 6
        query_name = "INTERVAL=[23,33);"
        sequence = (
            "AAAAAAAAAAAAAAAAAAAAAAACGGCTCGCATAGACACGACGACGACACGTACGATCGATCAGTCAT"
        )
        sam_string = f"{query_name}\t{flag}\t{ref_name}\t{pos}\t{mapq}\t{cigar}\t*\t0\t0\t{sequence}\t*\t{nm}\t{md}\tAS:i:0\tXS:i:0"
        record = pysam.AlignedSegment.fromstring(sam_string, header)
        classification = RecallClassification(record=record)

        assert not classification._whole_probe_maps()

    def test_isCorrect_probeIsSnpAndIsMismatchReturnsFalse(self):
        ref_name = "reference"
        ref_length = 64
        header = create_sam_header(ref_name, ref_length)
        flag = 0
        cigar = "56M"
        nm = "NM:i:1"
        md = "MD:Z:11A44"
        mapq = 60
        pos = 6
        query_name = "INTERVAL=[11,12);"
        sequence = "AAAAAAAAAAACGGCTCGCATAGACACGACGACGACACGTACGATCGATCAGTCAT"
        sam_string = f"{query_name}\t{flag}\t{ref_name}\t{pos}\t{mapq}\t{cigar}\t*\t0\t0\t{sequence}\t*\t{nm}\t{md}\tAS:i:0\tXS:i:0"
        record = pysam.AlignedSegment.fromstring(sam_string, header)
        probe = Probe(
            header=ProbeHeader.from_string(query_name), full_sequence=sequence
        )
        classification = RecallClassification(record=record)

        assert not classification.is_correct()

    def test_isCorrect_probeIsSnpAndIsHasMismatchToLeftReturnsTrue(self):
        ref_name = "reference"
        ref_length = 64
        header = create_sam_header(ref_name, ref_length)
        flag = 0
        cigar = "56M"
        nm = "NM:i:1"
        md = "MD:Z:10C45"
        mapq = 60
        pos = 6
        query_name = "INTERVAL=[11,12);"
        sequence = "AAAAAAAAAAACGGCTCGCATAGACACGACGACGACACGTACGATCGATCAGTCAT"
        sam_string = f"{query_name}\t{flag}\t{ref_name}\t{pos}\t{mapq}\t{cigar}\t*\t0\t0\t{sequence}\t*\t{nm}\t{md}\tAS:i:0\tXS:i:0"
        record = pysam.AlignedSegment.fromstring(sam_string, header)
        probe = Probe(
            header=ProbeHeader.from_string(query_name), full_sequence=sequence
        )
        classification = RecallClassification(record=record)

        assert classification.is_correct()

    def test_isCorrect_probeIsSnpAndIsHasMismatchToRightReturnsTrue(self):
        ref_name = "reference"
        ref_length = 64
        header = create_sam_header(ref_name, ref_length)
        flag = 0
        cigar = "56M"
        nm = "NM:i:1"
        md = "MD:Z:12C43"
        mapq = 60
        pos = 6
        query_name = "INTERVAL=[11,12);"
        sequence = "AAAAAAAAAAACGGCTCGCATAGACACGACGACGACACGTACGATCGATCAGTCAT"
        sam_string = f"{query_name}\t{flag}\t{ref_name}\t{pos}\t{mapq}\t{cigar}\t*\t0\t0\t{sequence}\t*\t{nm}\t{md}\tAS:i:0\tXS:i:0"
        record = pysam.AlignedSegment.fromstring(sam_string, header)
        probe = Probe(
            header=ProbeHeader.from_string(query_name), full_sequence=sequence
        )
        classification = RecallClassification(record=record)

        assert classification.is_correct()

    def test_isCorrect_probesMatchPerfectlyReturnsTrue(self):
        ref_name = "reference"
        ref_length = 59
        header = create_sam_header(ref_name, ref_length)
        flag = 0
        cigar = "33M35S"
        nm = "NM:i:0"
        md = "MD:Z:33"
        mapq = 60
        pos = 6
        query_name = "INTERVAL=[23,33);"
        sequence = (
            "AAAAAAAAAAAAAAAAAAAAAAACGGCTCGCATAGACACGACGACGACACGTACGATCGATCAGTCAT"
        )
        sam_string = f"{query_name}\t{flag}\t{ref_name}\t{pos}\t{mapq}\t{cigar}\t*\t0\t0\t{sequence}\t*\t{nm}\t{md}\tAS:i:0\tXS:i:0"
        record = pysam.AlignedSegment.fromstring(sam_string, header)
        probe = Probe(
            header=ProbeHeader.from_string(query_name), full_sequence=sequence
        )
        classification = RecallClassification(record=record)

        assert classification.is_correct()

    def test_isCorrect_probeIsDeletionBasesEitherSideMatchReturnsTrue(self):
        ref_name = "reference"
        ref_length = 64
        header = create_sam_header(ref_name, ref_length)
        flag = 0
        cigar = "56M"
        nm = "NM:i:0"
        md = "MD:Z:56"
        mapq = 60
        pos = 6
        query_name = "INTERVAL=[12,12);"
        sequence = "AAAAAAAAAAACGGCTCGCATAGACACGACGACGACACGTACGATCGATCAGTCAT"
        sam_string = f"{query_name}\t{flag}\t{ref_name}\t{pos}\t{mapq}\t{cigar}\t*\t0\t0\t{sequence}\t*\t{nm}\t{md}\tAS:i:0\tXS:i:0"
        record = pysam.AlignedSegment.fromstring(sam_string, header)
        probe = Probe(
            header=ProbeHeader.from_string(query_name), full_sequence=sequence
        )
        classification = RecallClassification(record=record)

        assert classification.is_correct()

    def test_isCorrect_probeIsDeletionBaseToLeftIsMismatchReturnsFalse(self):
        ref_name = "reference"
        ref_length = 64
        header = create_sam_header(ref_name, ref_length)
        flag = 0
        cigar = "56M"
        nm = "NM:i:1"
        md = "MD:Z:11T44"
        mapq = 60
        pos = 6
        query_name = "INTERVAL=[12,12);"
        sequence = "AAAAAAAAAAACGGCTCGCATAGACACGACGACGACACGTACGATCGATCAGTCAT"
        sam_string = f"{query_name}\t{flag}\t{ref_name}\t{pos}\t{mapq}\t{cigar}\t*\t0\t0\t{sequence}\t*\t{nm}\t{md}\tAS:i:0\tXS:i:0"
        record = pysam.AlignedSegment.fromstring(sam_string, header)
        probe = Probe(
            header=ProbeHeader.from_string(query_name), full_sequence=sequence
        )
        classification = RecallClassification(record=record)

        assert not classification.is_correct()

    def test_isCorrect_probeIsDeletionBaseToRightIsMismatchReturnsFalse(self):
        ref_name = "reference"
        ref_length = 64
        header = create_sam_header(ref_name, ref_length)
        flag = 0
        cigar = "56M"
        nm = "NM:i:1"
        md = "MD:Z:12T43"
        mapq = 60
        pos = 6
        query_name = "INTERVAL=[12,12);"
        sequence = "AAAAAAAAAAACGGCTCGCATAGACACGACGACGACACGTACGATCGATCAGTCAT"
        sam_string = f"{query_name}\t{flag}\t{ref_name}\t{pos}\t{mapq}\t{cigar}\t*\t0\t0\t{sequence}\t*\t{nm}\t{md}\tAS:i:0\tXS:i:0"
        record = pysam.AlignedSegment.fromstring(sam_string, header)
        probe = Probe(
            header=ProbeHeader.from_string(query_name), full_sequence=sequence
        )
        classification = RecallClassification(record=record)

        assert not classification.is_correct()

    def test_isCorrect_mismatchInFirstCoreProbeBaseReturnsFalse(self):
        ref_name = "reference"
        ref_length = 64
        header = create_sam_header(ref_name, ref_length)
        flag = 0
        cigar = "56M"
        nm = "NM:i:1"
        md = "MD:Z:11A44"
        mapq = 60
        pos = 1
        query_name = "INTERVAL=[11,21);"
        sequence = "AAAAAAAAAAACGGCTCGCATAGACACGACGACGACACGTACGATCGATCAGTCAT"
        sam_string = f"{query_name}\t{flag}\t{ref_name}\t{pos}\t{mapq}\t{cigar}\t*\t0\t0\t{sequence}\t*\t{nm}\t{md}\tAS:i:0\tXS:i:0"
        record = pysam.AlignedSegment.fromstring(sam_string, header)
        probe = Probe(
            header=ProbeHeader.from_string(query_name), full_sequence=sequence
        )
        classification = RecallClassification(record=record)

        assert not classification.is_correct()

    def test_isCorrect_mismatchInBaseBeforeCoreProbeStartsReturnsTrue(self):
        ref_name = "reference"
        ref_length = 64
        header = create_sam_header(ref_name, ref_length)
        flag = 0
        cigar = "56M"
        nm = "NM:i:1"
        md = "MD:Z:10T45"
        mapq = 60
        pos = 1
        query_name = "INTERVAL=[11,21);"
        sequence = "AAAAAAAAAAACGGCTCGCATAGACACGACGACGACACGTACGATCGATCAGTCAT"
        sam_string = f"{query_name}\t{flag}\t{ref_name}\t{pos}\t{mapq}\t{cigar}\t*\t0\t0\t{sequence}\t*\t{nm}\t{md}\tAS:i:0\tXS:i:0"
        record = pysam.AlignedSegment.fromstring(sam_string, header)
        probe = Probe(
            header=ProbeHeader.from_string(query_name), full_sequence=sequence
        )
        classification = RecallClassification(record=record)

        assert classification.is_correct()

    def test_isCorrect_mismatchInLastCoreProbeBaseReturnsFalse(self):
        ref_name = "reference"
        ref_length = 64
        header = create_sam_header(ref_name, ref_length)
        flag = 0
        cigar = "56M"
        nm = "NM:i:1"
        md = "MD:Z:20C35"
        mapq = 60
        pos = 1
        query_name = "INTERVAL=[11,21);"
        sequence = "AAAAAAAAAAACGGCTCGCATAGACACGACGACGACACGTACGATCGATCAGTCAT"
        sam_string = f"{query_name}\t{flag}\t{ref_name}\t{pos}\t{mapq}\t{cigar}\t*\t0\t0\t{sequence}\t*\t{nm}\t{md}\tAS:i:0\tXS:i:0"
        record = pysam.AlignedSegment.fromstring(sam_string, header)
        probe = Probe(
            header=ProbeHeader.from_string(query_name), full_sequence=sequence
        )
        classification = RecallClassification(record=record)

        assert not classification.is_correct()

    def test_isCorrect_mismatchInBaseAfterCoreProbeEndsReturnsTrue(self):
        ref_name = "reference"
        ref_length = 64
        header = create_sam_header(ref_name, ref_length)
        flag = 0
        cigar = "56M"
        nm = "NM:i:1"
        md = "MD:Z:21C34"
        mapq = 60
        pos = 1
        query_name = "INTERVAL=[11,21);"
        sequence = "AAAAAAAAAAACGGCTCGCATAGACACGACGACGACACGTACGATCGATCAGTCAT"
        sam_string = f"{query_name}\t{flag}\t{ref_name}\t{pos}\t{mapq}\t{cigar}\t*\t0\t0\t{sequence}\t*\t{nm}\t{md}\tAS:i:0\tXS:i:0"
        record = pysam.AlignedSegment.fromstring(sam_string, header)
        probe = Probe(
            header=ProbeHeader.from_string(query_name), full_sequence=sequence
        )
        classification = RecallClassification(record=record)

        assert classification.is_correct()

    def test_isCorrect_insertionInRefAfterFirstProbeCoreBaseReturnsFalse(self):
        ref_name = "reference"
        ref_length = 64
        header = create_sam_header(ref_name, ref_length)
        flag = 0
        cigar = "12M1D44M"
        nm = "NM:i:1"
        md = "MD:Z:12^G44"
        mapq = 60
        pos = 1
        query_name = "INTERVAL=[11,21);"
        sequence = "AAAAAAAAAAACGGCTCGCATAGACACGACGACGACACGTACGATCGATCAGTCAT"
        sam_string = f"{query_name}\t{flag}\t{ref_name}\t{pos}\t{mapq}\t{cigar}\t*\t0\t0\t{sequence}\t*\t{nm}\t{md}\tAS:i:0\tXS:i:0"
        record = pysam.AlignedSegment.fromstring(sam_string, header)
        probe = Probe(
            header=ProbeHeader.from_string(query_name), full_sequence=sequence
        )
        classification = RecallClassification(record=record)

        assert not classification.is_correct()

    def test_isCorrect_deletionInRefOfFirstProbeCoreBaseReturnsFalse(self):
        ref_name = "reference"
        ref_length = 63
        header = create_sam_header(ref_name, ref_length)
        flag = 0
        cigar = "11M1I44M"
        nm = "NM:i:1"
        md = "MD:Z:55"
        mapq = 60
        pos = 1
        query_name = "INTERVAL=[11,21);"
        sequence = "AAAAAAAAAAACGGCTCGCATAGACACGACGACGACACGTACGATCGATCAGTCAT"
        sam_string = f"{query_name}\t{flag}\t{ref_name}\t{pos}\t{mapq}\t{cigar}\t*\t0\t0\t{sequence}\t*\t{nm}\t{md}\tAS:i:0\tXS:i:0"
        record = pysam.AlignedSegment.fromstring(sam_string, header)
        probe = Probe(
            header=ProbeHeader.from_string(query_name), full_sequence=sequence
        )
        classification = RecallClassification(record=record)

        assert not classification.is_correct()

    def test_isCorrect_insertionInRefBeforeLastProbeCoreBaseReturnsFalse(self):
        ref_name = "reference"
        ref_length = 65
        header = create_sam_header(ref_name, ref_length)
        flag = 0
        cigar = "20M1D36M"
        nm = "NM:i:1"
        md = "MD:Z:20^C36"
        mapq = 60
        pos = 1
        query_name = "INTERVAL=[11,21);"
        sequence = "AAAAAAAAAAACGGCTCGCATAGACACGACGACGACACGTACGATCGATCAGTCAT"
        sam_string = f"{query_name}\t{flag}\t{ref_name}\t{pos}\t{mapq}\t{cigar}\t*\t0\t0\t{sequence}\t*\t{nm}\t{md}\tAS:i:0\tXS:i:0"
        record = pysam.AlignedSegment.fromstring(sam_string, header)
        probe = Probe(
            header=ProbeHeader.from_string(query_name), full_sequence=sequence
        )
        classification = RecallClassification(record=record)

        assert not classification.is_correct()

    def test_isCorrect_insertionInRefAfterLastProbeCoreBaseReturnsTrue(self):
        ref_name = "reference"
        ref_length = 65
        header = create_sam_header(ref_name, ref_length)
        flag = 0
        cigar = "21M1D35M"
        nm = "NM:i:1"
        md = "MD:Z:21^G35"
        mapq = 60
        pos = 1
        query_name = "INTERVAL=[11,21);"
        sequence = "AAAAAAAAAAACGGCTCGCATAGACACGACGACGACACGTACGATCGATCAGTCAT"
        sam_string = f"{query_name}\t{flag}\t{ref_name}\t{pos}\t{mapq}\t{cigar}\t*\t0\t0\t{sequence}\t*\t{nm}\t{md}\tAS:i:0\tXS:i:0"
        record = pysam.AlignedSegment.fromstring(sam_string, header)
        probe = Probe(
            header=ProbeHeader.from_string(query_name), full_sequence=sequence
        )
        classification = RecallClassification(record=record)

        assert classification.is_correct()

    def test_isCorrect_deletionInRefOfLastProbeCoreBaseReturnsFalse(self):
        ref_name = "reference"
        ref_length = 63
        header = create_sam_header(ref_name, ref_length)
        flag = 0
        cigar = "20M1I35M"
        nm = "NM:i:1"
        md = "MD:Z:55"
        mapq = 60
        pos = 1
        query_name = "INTERVAL=[11,21);"
        sequence = "AAAAAAAAAAACGGCTCGCATAGACACGACGACGACACGTACGATCGATCAGTCAT"
        sam_string = f"{query_name}\t{flag}\t{ref_name}\t{pos}\t{mapq}\t{cigar}\t*\t0\t0\t{sequence}\t*\t{nm}\t{md}\tAS:i:0\tXS:i:0"
        record = pysam.AlignedSegment.fromstring(sam_string, header)
        probe = Probe(
            header=ProbeHeader.from_string(query_name), full_sequence=sequence
        )
        classification = RecallClassification(record=record)

        assert not classification.is_correct()

    def test_assessment_unmappedRecordReturnsUnmapped(self):
        record = create_unmapped_sam_record()
        classification = RecallClassification(record=record)

        actual = classification.assessment()
        expected = "unmapped"

        assert actual == expected

    def test_assessment_recordWithCoreProbeOnlyPartiallyMappedReturnsPartiallyMapped(
        self
    ):
        record = create_partially_mapped_sam_record()
        classification = RecallClassification(record=record)

        actual = classification.assessment()
        expected = "partially_mapped"

        assert actual == expected

    def test_assessment_incorrectSecondayAlignmentMismatchReturnsIncorrect(self):
        record = create_incorrect_secondary_sam_record()
        classification = RecallClassification(record=record)

        actual = classification.assessment()
        expected = "secondary_incorrect"

        assert actual == expected

    def test_assessment_correctSecondayAlignmentReturnsCorrect(self):
        record = create_correct_secondary_sam_record()
        classification = RecallClassification(record=record)

        actual = classification.assessment()
        expected = "secondary_correct"

        assert actual == expected

    def test_assessment_incorrectSupplementaryAlignmentMismatchReturnsIncorrect(self):
        record = create_incorrect_supplementary_sam_record()
        classification = RecallClassification(record=record)

        actual = classification.assessment()
        expected = "supplementary_incorrect"

        assert actual == expected

    def test_assessment_correctSupplementaryAlignmentReturnsCorrect(self):
        record = create_correct_supplementary_sam_record()
        classification = RecallClassification(record=record)

        actual = classification.assessment()
        expected = "supplementary_correct"

        assert actual == expected

    def test_assessment_correctPrimaryAlignmenReturnsCorrect(self):
        record = create_correct_primary_sam_record()
        classification = RecallClassification(record=record)

        actual = classification.assessment()
        expected = "correct"

        assert actual == expected

    def test_assessment_incorrectPrimaryAlignmentMismatchReturnsIncorrect(self):
        record = create_incorrect_primary_sam_record()
        classification = RecallClassification(record=record)

        actual = classification.assessment()
        expected = "incorrect"

        assert actual == expected


class TestRecallClassifier:
    def test_classify_noAlignmentFileReturnsEmpty(self):
        classifier = RecallClassifier()

        actual = classifier.classify()
        expected = []

        assert actual == expected

    def test_classify(self):
        classifier = create_recall_classifier_with_two_entries()

        actual = classifier.classify()
        expected = [
            RecallClassification(record=create_correct_primary_sam_record()),
            RecallClassification(record=create_incorrect_supplementary_sam_record()),
        ]

        assert actual == expected

        expected_classifications = ["correct", "supplementary_incorrect"]
        assert [x.assessment() for x in actual] == expected_classifications


class TestRecallReporter:
    def test_generateReport_noClassifierReturnsEmpty(self):
        sample = "sample"
        classifier = RecallClassifier(name=sample)
        reporter = RecallReporter(classifiers=[classifier])

        actual = reporter.generate_report()
        expected = pd.DataFrame(
            [],
            columns=[
                "sample",
                "truth_probe_header",
                "vcf_probe_header",
                "classification",
            ],
        )

        assert actual.equals(expected)

    def test_generateReport_twoClassificationsReturnsDataframeWithTwoEntries(self):
        classifier = create_recall_classifier_with_two_entries()
        sample = "sample"
        classifier.name = sample
        reporter = RecallReporter(classifiers=[classifier])

        actual = reporter.generate_report()
        expected_data = []
        for assessment, record in [
            ("correct", create_correct_primary_sam_record()),
            ("supplementary_incorrect", create_incorrect_supplementary_sam_record()),
        ]:
            expected_data.append(
                [sample, record.query_name, record.reference_name, assessment]
            )
        expected = pd.DataFrame(
            expected_data,
            columns=[
                "sample",
                "truth_probe_header",
                "vcf_probe_header",
                "classification",
            ],
        )

        assert actual.equals(expected)

    def test_save_emptyReporterReturnsHeadersOnly(self):
        delim = "\t"
        reporter = RecallReporter(classifiers=[RecallClassifier()], delim=delim)
        fh = StringIO(newline="")
        reporter.save(fh)

        fh.seek(0)
        actual = fh.read()
        expected = delim.join(reporter.columns) + "\n"

        assert actual == expected

    def test_save_reporterWithTwoClassificationsWritesHeadersAndTwoRows(self):
        primary_correct_record = create_correct_primary_sam_record()
        suppl_incorrect_record = create_incorrect_supplementary_sam_record()
        delim = "\t"
        classifier = create_recall_classifier_with_two_entries()
        sample = "sample"
        classifier.name = sample
        reporter = RecallReporter(classifiers=[classifier], delim=delim)

        fh = StringIO(newline="")
        reporter.save(fh)

        fh.seek(0)
        actual = fh.read()
        expected_data = []
        for assessment, record in [
            ("correct", primary_correct_record),
            ("supplementary_incorrect", suppl_incorrect_record),
        ]:
            expected_data.append(
                [sample, record.query_name, record.reference_name, assessment]
            )
        expected = StringIO(newline="")
        pd.DataFrame(
            expected_data,
            columns=[
                "sample",
                "truth_probe_header",
                "vcf_probe_header",
                "classification",
            ],
        ).to_csv(expected, sep=delim, header=True, index=False)
        expected.seek(0)
        expected = expected.read()

        assert actual == expected

    def test_save_reporterWithTwoClassificationsWritesHeadersAndTwoRowsWithCommaDelim(
        self
    ):
        primary_correct_record = create_correct_primary_sam_record()
        suppl_incorrect_record = create_incorrect_supplementary_sam_record()
        delim = ","
        classifier = create_recall_classifier_with_two_entries()
        sample = "sample"
        classifier.name = sample
        reporter = RecallReporter(classifiers=[classifier], delim=delim)

        fh = StringIO(newline="")
        reporter.save(fh)

        fh.seek(0)
        actual = fh.read()
        expected_data = []
        for assessment, record in [
            ("correct", primary_correct_record),
            ("supplementary_incorrect", suppl_incorrect_record),
        ]:
            expected_data.append(
                [sample, record.query_name, record.reference_name, assessment]
            )
        expected = StringIO(newline="")
        pd.DataFrame(
            expected_data,
            columns=[
                "sample",
                "truth_probe_header",
                "vcf_probe_header",
                "classification",
            ],
        ).to_csv(expected, sep=delim, header=True, index=False)
        expected.seek(0)
        expected = expected.read()

        assert actual == expected

    def test_save_reporterWithTwoClassifiersWritesTwoSamplesWithTwoRows(self):
        primary_correct_record = create_correct_primary_sam_record()
        suppl_incorrect_record = create_incorrect_supplementary_sam_record()
        delim = ","
        classifier1 = create_recall_classifier_with_two_entries()
        sample = "sample"
        classifier1.name = sample
        classifier2 = create_recall_classifier_with_two_entries()
        sample2 = "sample2"
        classifier2.name = sample2

        reporter = RecallReporter(classifiers=[classifier1, classifier2], delim=delim)

        fh = StringIO(newline="")
        reporter.save(fh)

        fh.seek(0)
        actual = fh.read()
        expected_data = []
        for s in [sample, sample2]:
            for assessment, record in [
                ("correct", primary_correct_record),
                ("supplementary_incorrect", suppl_incorrect_record),
            ]:
                expected_data.append(
                    [s, record.query_name, record.reference_name, assessment]
                )
        expected = StringIO(newline="")
        pd.DataFrame(
            expected_data,
            columns=[
                "sample",
                "truth_probe_header",
                "vcf_probe_header",
                "classification",
            ],
        ).to_csv(expected, sep=delim, header=True, index=False)
        expected.seek(0)
        expected = expected.read()

        assert actual == expected


def create_report_row(
    classification: str, gt_conf: float = 0, sample: str = "sample1"
) -> pd.Series:
    truth_probe_header = ProbeHeader()
    vcf_probe_header = ProbeHeader(gt_conf=gt_conf)
    data = {
        "sample": sample,
        "truth_probe_header": str(truth_probe_header),
        "vcf_probe_header": str(vcf_probe_header),
        "classification": classification,
    }
    return pd.Series(data=data)


class TestRecallCalculator:
    def test_statisticalClassification_unmappedReturnsFalseNegative(self):
        classification = "unmapped"
        row = create_report_row(classification)

        actual = RecallCalculator.statistical_classification(row)
        expected = StatisticalClassification.FALSE_NEGATIVE

        assert actual == expected

    def test_statisticalClassification_partiallyMappedReturnsFalseNegative(self):
        classification = "partially_mapped"
        row = create_report_row(classification)

        actual = RecallCalculator.statistical_classification(row)
        expected = StatisticalClassification.FALSE_NEGATIVE

        assert actual == expected


    def test_statisticalClassification_incorrectAboveConfReturnsFalsePositive(self):
        classification = "incorrect"
        gt_conf = 10
        row = create_report_row(classification, gt_conf)
        threshold = 5

        actual = RecallCalculator.statistical_classification(row, conf_threshold=threshold)
        expected = StatisticalClassification.FALSE_POSITIVE

        assert actual == expected

    def test_statisticalClassification_incorrectBelowConfReturnsFalseNegative(self):
        classification = "incorrect"
        gt_conf = 10
        row = create_report_row(classification, gt_conf)
        threshold = 50

        actual = RecallCalculator.statistical_classification(row, conf_threshold=threshold)
        expected = StatisticalClassification.FALSE_NEGATIVE

        assert actual == expected

    def test_statisticalClassification_correctAboveConfReturnsTruePositive(self):
        classification = "correct"
        gt_conf = 10
        row = create_report_row(classification, gt_conf)
        threshold = 5

        actual = RecallCalculator.statistical_classification(row, conf_threshold=threshold)
        expected = StatisticalClassification.TRUE_POSITIVE

        assert actual == expected

    def test_statisticalClassification_correctBelowConfReturnsFalseNegative(self):
        classification = "correct"
        gt_conf = 10
        row = create_report_row(classification, gt_conf)
        threshold = 50

        actual = RecallCalculator.statistical_classification(row, conf_threshold=threshold)
        expected = StatisticalClassification.FALSE_NEGATIVE

        assert actual == expected

    def test_statisticalClassification_secondaryIncorrectAboveConfReturnsFalsePositive(
        self
    ):
        classification = "secondary_incorrect"
        gt_conf = 10
        row = create_report_row(classification, gt_conf)
        threshold = 5

        actual = RecallCalculator.statistical_classification(row, conf_threshold=threshold)
        expected = StatisticalClassification.FALSE_POSITIVE

        assert actual == expected

    def test_statisticalClassification_secondaryIncorrectBelowConfReturnsFalseNegative(
        self
    ):
        classification = "secondary_incorrect"
        gt_conf = 10
        row = create_report_row(classification, gt_conf)
        threshold = 50

        actual = RecallCalculator.statistical_classification(row, conf_threshold=threshold)
        expected = StatisticalClassification.FALSE_NEGATIVE

        assert actual == expected

    def test_statisticalClassification_secondaryCorrectAboveConfReturnsTruePositive(
        self
    ):
        classification = "secondary_correct"
        gt_conf = 10
        row = create_report_row(classification, gt_conf)
        threshold = 5

        actual = RecallCalculator.statistical_classification(row, conf_threshold=threshold)
        expected = StatisticalClassification.TRUE_POSITIVE

        assert actual == expected

    def test_statisticalClassification_secondaryCorrectBelowConfReturnsFalseNegative(
        self
    ):
        classification = "secondary_correct"
        gt_conf = 10
        row = create_report_row(classification, gt_conf)
        threshold = 50

        actual = RecallCalculator.statistical_classification(row, conf_threshold=threshold)
        expected = StatisticalClassification.FALSE_NEGATIVE

        assert actual == expected

    def test_statisticalClassification_supplementaryIncorrectAboveConfReturnsFalsePositive(
        self
    ):
        classification = "supplementary_incorrect"
        gt_conf = 10
        row = create_report_row(classification, gt_conf)
        threshold = 5

        actual = RecallCalculator.statistical_classification(row, conf_threshold=threshold)
        expected = StatisticalClassification.FALSE_POSITIVE

        assert actual == expected

    def test_statisticalClassification_supplementaryIncorrectBelowConfReturnsFalseNegative(
        self
    ):
        classification = "supplementary_incorrect"
        gt_conf = 10
        row = create_report_row(classification, gt_conf)
        threshold = 50

        actual = RecallCalculator.statistical_classification(row, conf_threshold=threshold)
        expected = StatisticalClassification.FALSE_NEGATIVE

        assert actual == expected

    def test_statisticalClassification_supplementaryCorrectAboveConfReturnsTruePositive(
        self
    ):
        classification = "supplementary_correct"
        gt_conf = 10
        row = create_report_row(classification, gt_conf)
        threshold = 5

        actual = RecallCalculator.statistical_classification(row, conf_threshold=threshold)
        expected = StatisticalClassification.TRUE_POSITIVE

        assert actual == expected

    def test_statisticalClassification_supplementaryCorrectBelowConfReturnsFalseNegative(
        self
    ):
        classification = "supplementary_correct"
        gt_conf = 10
        row = create_report_row(classification, gt_conf)
        threshold = 50

        actual = RecallCalculator.statistical_classification(row, conf_threshold=threshold)
        expected = StatisticalClassification.FALSE_NEGATIVE

        assert actual == expected
