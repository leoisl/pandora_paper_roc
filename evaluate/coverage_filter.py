from evaluate.filter import Filter
from .vcf import VCF


class CoverageFilter(Filter):
    def __init__(self, coverage_threshold: float):
        self._coverage_threshold = coverage_threshold

    @property
    def coverage_threshold(self) -> float:
        return self._coverage_threshold

    def record_should_be_filtered_out(self, record: VCF) -> bool:
        return record.mean_coverage < self.coverage_threshold
