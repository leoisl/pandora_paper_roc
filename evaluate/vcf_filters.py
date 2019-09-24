from .vcf import VCF
from collections import UserList
from .coverage_filter import CoverageFilter
from .strand_bias_filter import StrandBiasFilter


class VCF_Filters(UserList):
    def record_should_be_filtered_out(self, vcf_record: VCF) -> bool:
        return any(
            vcf_filter.record_should_be_filtered_out(vcf_record) for vcf_filter in self
        )

    @staticmethod
    def get_list_of_all_VCF_Filters(coverage_threshold: float, strand_bias_threshold: float) -> "VCF_Filters":
        vcf_filters = VCF_Filters()
        vcf_filters.append(CoverageFilter(coverage_threshold))
        vcf_filters.append(StrandBiasFilter(strand_bias_threshold))
        return vcf_filters
