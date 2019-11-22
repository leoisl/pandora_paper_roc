import pysam
from typing import List, Dict
from .vcf import VCF, BuggedVCFError, NullVCFError
from collections import defaultdict


class VCFFile:
    def __init__(self, pysam_variant_file: pysam.VariantFile):
        self._sample_to_gene_to_VCFs = defaultdict(lambda: defaultdict(list))
        for variant_record in pysam_variant_file:
            for sample in variant_record.samples:
                gene = variant_record.chrom
                try:
                    vcf = VCF.from_VariantRecord_and_Sample(variant_record, sample)
                    self._sample_to_gene_to_VCFs[sample][gene].append(vcf)
                except NullVCFError:
                    pass
                except BuggedVCFError:
                    assert False, f"We found a bugged VCF: {vcf.variant}"

    @property
    def sample_to_gene_to_VCFs(self) -> Dict[str, Dict[str, List[VCF]]]:
        return self._sample_to_gene_to_VCFs

    def get_VCF_records_given_sample_and_gene(
        self, sample: str, gene_name: str
    ) -> List[VCF]:
        return self.sample_to_gene_to_VCFs[sample][gene_name]
