from evaluate.filter import Filter
from .vcf import VCF


class StrandBiasFilter(Filter):
    def __init__(self, strand_bias_threshold: float):
        self._strand_bias_threshold = strand_bias_threshold

    @property
    def strand_bias_threshold(self) -> float:
        return self._strand_bias_threshold

    def record_should_be_filtered_out(self, record: VCF) -> bool:
        if record.mean_coverage == 0:
            return True
        strand_ratio = record.mean_coverage_forward / record.mean_coverage
        bad_strand_ratio = (
            strand_ratio <= self.strand_bias_threshold
            or strand_ratio >= 1.0 - self.strand_bias_threshold
        )
        return bad_strand_ratio
