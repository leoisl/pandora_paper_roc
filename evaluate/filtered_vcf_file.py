from .vcf_file import VCFFile
import pysam
from .vcf_filters import VCF_Filters
from typing import Dict, List
from .vcf import VCF
from collections import defaultdict


class FilteredVCFFile(VCFFile):
    def __init__(self, pysam_variant_file: pysam.VariantFile, filters: VCF_Filters, VCF_creator_method, isolated_variants, flank_length):
        VCFFile.__init__(self, pysam_variant_file, VCF_creator_method, isolated_variants, flank_length)
        self._sample_to_gene_to_VCFs = FilteredVCFFile._filter_records(
            self.sample_to_gene_to_VCFs, filters
        )

    @staticmethod
    def _filter_records(
        sample_to_gene_to_VCFs_all_records: Dict[str, Dict[str, List[VCF]]],
        filters: VCF_Filters,
    ) -> Dict[str, Dict[str, List[VCF]]]:
        sample_to_gene_to_VCFs_filtered_records = defaultdict(lambda: defaultdict(list))
        for (
            sample,
            gene_to_VCFs_of_a_sample,
        ) in sample_to_gene_to_VCFs_all_records.items():
            for gene, vcfs in gene_to_VCFs_of_a_sample.items():
                filtered_VCFs = [
                    vcf
                    for vcf in vcfs
                    if not filters.record_should_be_filtered_out(vcf)
                ]
                if len(filtered_VCFs) > 0:
                    sample_to_gene_to_VCFs_filtered_records[sample][gene] = filtered_VCFs
        return sample_to_gene_to_VCFs_filtered_records
