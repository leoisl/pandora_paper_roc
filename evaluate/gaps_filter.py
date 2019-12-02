from evaluate.filter import Filter
from .vcf import PandoraVCF


class GapsFilter(Filter):
    def __init__(self, gaps_threshold: float):
        self._gaps_threshold = gaps_threshold

    @property
    def gaps_threshold(self) -> float:
        return self._gaps_threshold

    def record_should_be_filtered_out(self, record: PandoraVCF) -> bool:
        return record._gaps > self.gaps_threshold
