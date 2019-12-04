import pysam
from typing import List, Dict
from .vcf import VCF, BuggedVCFError, NullVCFError, VCFFactory
from collections import defaultdict


class VCFFile:
    def __init__(self, pysam_variant_file: pysam.VariantFile, VCF_creator_method, isolated_variants, flank_length):
        self._sample_to_gene_to_VCFs = defaultdict(lambda: defaultdict(list))
        self._sample_to_gene_to_positions_covered_by_added_records = defaultdict(lambda: defaultdict(set))
        for variant_record in pysam_variant_file:
            for sample in variant_record.samples:
                gene = variant_record.chrom
                try:
                    vcf = VCF_creator_method(variant_record, sample)

                    if isolated_variants and self.record_overlaps_an_already_covered_position(vcf, sample, gene, flank_length):
                        print(f"Removing record {vcf}")
                        continue

                    self._sample_to_gene_to_VCFs[sample][gene].append(vcf)
                    self.add_positions_covered_by_variant(vcf, sample, gene, flank_length)
                except NullVCFError:
                    pass
                except BuggedVCFError:
                    assert False, f"We found a bugged VCF: {variant_record}"

    @property
    def sample_to_gene_to_VCFs(self) -> Dict[str, Dict[str, List[VCF]]]:
        return self._sample_to_gene_to_VCFs

    def record_overlaps_an_already_covered_position(self, vcf_record, sample, gene, flank_length) -> bool:
        for position_in_variant in vcf_record.positions_covered_by_variant_with_flanks(flank_length):
            if position_in_variant in self._sample_to_gene_to_positions_covered_by_added_records[sample][gene]:
                return True
        return False

    def add_positions_covered_by_variant(self, vcf_record, sample, gene, flank_length):
        for position_in_variant in vcf_record.positions_covered_by_variant_with_flanks(flank_length):
            self._sample_to_gene_to_positions_covered_by_added_records[sample][gene].add(position_in_variant)

    def get_VCF_records_given_sample_and_gene(
        self, sample: str, gene_name: str
    ) -> List[VCF]:
        return self.sample_to_gene_to_VCFs[sample][gene_name]
