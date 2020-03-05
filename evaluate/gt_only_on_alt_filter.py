from evaluate.filter import Filter
from .vcf import VCF


class GT_Only_on_Alt_filter(Filter):
    def record_should_be_filtered_out(self, record: VCF) -> bool:
        return record.genotype == 0 or record.genotype == "."
