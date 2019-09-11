from tests.test_evaluation import create_sam_header
from tests.common import (
    create_unmapped_sam_record,
    create_partially_mapped_sam_record,
    create_incorrect_secondary_sam_record,
    create_correct_secondary_sam_record,
    create_incorrect_supplementary_sam_record,
    create_correct_supplementary_sam_record,
    create_correct_primary_sam_record,
    create_incorrect_primary_sam_record,
)
from evaluate.classification import *
import pytest


class TestClassification:
    def test_equality_twoEqualReturnsTrue(self):
        c1 = Classification()
        c1.query_probe = Probe(ProbeHeader(chrom="2", pos=5))
        c2 = Classification()
        c2.query_probe = Probe(ProbeHeader(chrom="2", pos=5))

        assert c1 == c2

    def test_equality_twoNonEqualReturnsFalse(self):
        c1 = Classification()
        c1.query_probe = Probe(ProbeHeader(chrom="2", pos=5))
        c2 = Classification()
        c2.query_probe = Probe(ProbeHeader(chrom="3", pos=5))

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
        classification = Classification(record=record)

        assert not classification._whole_query_probe_maps()

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
        classification = Classification(record=record)

        assert classification._whole_query_probe_maps()

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
        classification = Classification(record=record)

        assert classification._whole_query_probe_maps()

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
        classification = Classification(record=record)

        assert not classification._whole_query_probe_maps()

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
        classification = Classification(record=record)

        assert not classification._whole_query_probe_maps()

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
        classification = Classification(record=record)

        assert not classification._whole_query_probe_maps()

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
        classification = Classification(record=record)

        assert not classification._whole_query_probe_maps()

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
        classification = Classification(record=record)

        assert classification._whole_query_probe_maps()

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
        classification = Classification(record=record)

        assert not classification._whole_query_probe_maps()

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
        classification = Classification(record=record)

        assert not classification._whole_query_probe_maps()

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
        classification = Classification(record=record)

        assert not classification._whole_query_probe_maps()

    def test_assessment_raisesNotImplementedError(self):
        classification = Classification()
        with pytest.raises(NotImplementedError):
            classification.assessment()


class TestRecallClassification:
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
        expected = AlignmentAssessment.UNMAPPED

        assert actual == expected

    def test_assessment_recordWithCoreProbeOnlyPartiallyMappedReturnsPartiallyMapped(
        self
    ):
        record = create_partially_mapped_sam_record()
        classification = RecallClassification(record=record)

        actual = classification.assessment()
        expected = AlignmentAssessment.PARTIALLY_MAPPED

        assert actual == expected

    def test_assessment_incorrectSecondayAlignmentMismatchReturnsIncorrect(self):
        record = create_incorrect_secondary_sam_record()
        classification = RecallClassification(record=record)

        actual = classification.assessment()
        expected = AlignmentAssessment.SECONDARY_INCORRECT

        assert actual == expected

    def test_assessment_correctSecondayAlignmentReturnsCorrect(self):
        record = create_correct_secondary_sam_record()
        classification = RecallClassification(record=record)

        actual = classification.assessment()
        expected = AlignmentAssessment.SECONDARY_CORRECT

        assert actual == expected

    def test_assessment_incorrectSupplementaryAlignmentMismatchReturnsIncorrect(self):
        record = create_incorrect_supplementary_sam_record()
        classification = RecallClassification(record=record)

        actual = classification.assessment()
        expected = AlignmentAssessment.SUPPLEMENTARY_INCORRECT

        assert actual == expected

    def test_assessment_correctSupplementaryAlignmentReturnsCorrect(self):
        record = create_correct_supplementary_sam_record()
        classification = RecallClassification(record=record)

        actual = classification.assessment()
        expected = AlignmentAssessment.SUPPLEMENTARY_CORRECT

        assert actual == expected

    def test_assessment_correctPrimaryAlignmenReturnsCorrect(self):
        record = create_correct_primary_sam_record()
        classification = RecallClassification(record=record)

        actual = classification.assessment()
        expected = AlignmentAssessment.PRIMARY_CORRECT

        assert actual == expected

    def test_assessment_incorrectPrimaryAlignmentMismatchReturnsIncorrect(self):
        record = create_incorrect_primary_sam_record()
        classification = RecallClassification(record=record)

        actual = classification.assessment()
        expected = AlignmentAssessment.PRIMARY_INCORRECT

        assert actual == expected


class TestPrecisionClassification:
    def test_assessment_unmappedRecordReturnsZero(self):
        record = create_unmapped_sam_record()
        classification = PrecisionClassification(record=record)

        actual = classification.assessment()
        expected = 0.0

        assert actual == expected

    def test_assessment_recordWithCoreProbeOnlyPartiallyMappedReturnsZero(self):
        record = create_partially_mapped_sam_record()
        classification = PrecisionClassification(record=record)

        actual = classification.assessment()
        expected = 0.0

        assert actual == expected

    def test_assessment_probeIsSnpAndIsMismatch(self):
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
        classification = PrecisionClassification(record=record)

        expected = 0.0
        actual = classification.assessment()
        assert actual == expected

    def test_assessment_probeIsSnpAndIsHasMismatchToLeft(self):
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
        classification = PrecisionClassification(record=record)

        expected = 1.0
        actual = classification.assessment()
        assert actual == expected

    def test_assessment_probeIsSnpAndIsHasMismatchToRight(self):
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
        classification = PrecisionClassification(record=record)

        expected = 1.0
        actual = classification.assessment()
        assert actual == expected

    def test_assessment_probesMatchPerfectly(self):
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
        classification = PrecisionClassification(record=record)

        expected = 1.0
        actual = classification.assessment()
        assert actual == expected

    def test_assessment_probeIsDeletionBasesEitherSideMatch(self):
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
        classification = PrecisionClassification(record=record)

        expected = 1.0
        actual = classification.assessment()
        assert actual == expected

    def test_assessment_probeIsDeletionBaseToLeftIsMismatch(self):
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
        classification = PrecisionClassification(record=record)

        expected = 0.5
        actual = classification.assessment()
        assert actual == expected

    def test_assessment_probeIsDeletionBaseToRightIsMismatch(self):
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
        classification = PrecisionClassification(record=record)

        expected = 0.5
        actual = classification.assessment()
        assert actual == expected

    def test_assessment_mismatchInFirstCoreProbeBase(self):
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
        classification = PrecisionClassification(record=record)

        expected = 0.9
        actual = classification.assessment()
        assert actual == expected

    def test_assessment_mismatchInBaseBeforeCoreProbeStarts(self):
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
        classification = PrecisionClassification(record=record)

        expected = 1.0
        actual = classification.assessment()
        assert actual == expected

    def test_assessment_mismatchInLastCoreProbeBase(self):
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
        classification = PrecisionClassification(record=record)

        expected = 0.9
        actual = classification.assessment()
        assert actual == expected

    def test_assessment_mismatchInBaseAfterCoreProbeEnds(self):
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
        classification = PrecisionClassification(record=record)

        expected = 1.0
        actual = classification.assessment()
        assert actual == expected

    def test_assessment_insertionInRefAfterFirstProbeCoreBase(self):
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
        classification = PrecisionClassification(record=record)

        expected = 10 / 11
        actual = classification.assessment()
        assert actual == expected

    def test_assessment_deletionInRefOfFirstProbeCoreBase(self):
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
        classification = PrecisionClassification(record=record)

        expected = 9 / 10
        actual = classification.assessment()
        assert actual == expected

    def test_assessment_insertionInRefBeforeLastProbeCoreBase(self):
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
        classification = PrecisionClassification(record=record)

        expected = 10 / 11
        actual = classification.assessment()
        assert actual == expected

    def test_assessment_insertionInRefAfterLastProbeCoreBase(self):
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
        classification = PrecisionClassification(record=record)

        expected = 1.0
        actual = classification.assessment()
        assert actual == expected

    def test_assessment_deletionInRefOfLastProbeCoreBase(self):
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
        classification = PrecisionClassification(record=record)

        expected = 9 / 10
        actual = classification.assessment()
        assert actual == expected

    def test_assessment_mismatch4ProbeBases(self):
        ref_name = "reference"
        ref_length = 64
        header = create_sam_header(ref_name, ref_length)
        flag = 0
        cigar = "56M"
        nm = "NM:i:1"
        md = "MD:Z:12A2AGC38"
        mapq = 60
        pos = 1
        query_name = "INTERVAL=[11,21);"
        sequence = "AAAAAAAAAAACGGCTCGCATAGACACGACGACGACACGTACGATCGATCAGTCAT"
        sam_string = f"{query_name}\t{flag}\t{ref_name}\t{pos}\t{mapq}\t{cigar}\t*\t0\t0\t{sequence}\t*\t{nm}\t{md}\tAS:i:0\tXS:i:0"
        record = pysam.AlignedSegment.fromstring(sam_string, header)
        probe = Probe(
            header=ProbeHeader.from_string(query_name), full_sequence=sequence
        )
        classification = PrecisionClassification(record=record)

        expected = 0.6
        actual = classification.assessment()
        assert actual == expected

    def test_assessment_mismatch1insertion1deletion1ProbeBases(self):
        ref_name = "reference"
        ref_length = 56
        header = create_sam_header(ref_name, ref_length)
        flag = 0
        cigar = "11M1I8M1D36M"
        nm = "NM:i:3"
        md = "MD:Z:15G3^T36"
        mapq = 60
        pos = 1
        query_name = "INTERVAL=[11,21);"
        sequence = "AAAAAAAAAAACGGCTCGCATAGACACGACGACGACACGTACGATCGATCAGTCAT"
        sam_string = f"{query_name}\t{flag}\t{ref_name}\t{pos}\t{mapq}\t{cigar}\t*\t0\t0\t{sequence}\t*\t{nm}\t{md}\tAS:i:0\tXS:i:0"
        record = pysam.AlignedSegment.fromstring(sam_string, header)
        probe = Probe(
            header=ProbeHeader.from_string(query_name), full_sequence=sequence
        )
        classification = PrecisionClassification(record=record)

        expected = 8 / 11
        actual = classification.assessment()
        assert actual == expected
