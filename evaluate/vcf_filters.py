from .vcf import VCF
from collections import UserList
from .coverage_filter import CoverageFilter
from .strand_bias_filter import StrandBiasFilter
from .gaps_filter import GapsFilter


class VCF_Filters(UserList):
    def record_should_be_filtered_out(self, vcf_record: VCF) -> bool:
        return any(
            vcf_filter.record_should_be_filtered_out(vcf_record) for vcf_filter in self
        )

    @staticmethod
    def get_all_VCF_Filters(
        coverage_threshold: str, strand_bias_threshold: str, gaps_threshold: str
    ) -> "VCF_Filters":
        vcf_filters = VCF_Filters()

        if coverage_threshold != "Not_App":
            vcf_filters.append(CoverageFilter(float(coverage_threshold)))

        if strand_bias_threshold != "Not_App":
            vcf_filters.append(StrandBiasFilter(float(strand_bias_threshold)))

        if gaps_threshold != "Not_App":
            vcf_filters.append(GapsFilter(float(gaps_threshold)))

        return vcf_filters
