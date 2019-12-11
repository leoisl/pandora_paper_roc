from evaluate.filter import Filter
from .vcf import PandoraVCF


class StrandBiasFilter(Filter):
    def __init__(self, strand_bias_threshold: float):
        self._strand_bias_threshold = strand_bias_threshold

    @property
    def strand_bias_threshold(self) -> float:
        return self._strand_bias_threshold

    def record_should_be_filtered_out(self, record: PandoraVCF) -> bool:
        if record.coverage == 0:
            return True
        strand_ratio = record._mean_coverage_forward / record.coverage
        bad_strand_ratio = (
            strand_ratio < self.strand_bias_threshold
            or strand_ratio > 1.0 - self.strand_bias_threshold
        )
        return bad_strand_ratio
