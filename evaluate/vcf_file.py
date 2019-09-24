from pathlib import Path
import pysam
from typing import List
from .vcf import VCF
from collections import defaultdict


class VCFFile:
    def __init__(self, vcf_filepath: Path):
        self._sample_to_gene_to_VCFs = defaultdict(lambda: defaultdict(list))
        with pysam.VariantFile(vcf_filepath) as pysam_variant_file:
            for variant_record in pysam_variant_file:
                for sample in variant_record.samples:
                    gene = variant_record.chrom
                    self._sample_to_gene_to_VCFs[sample][gene].append(
                        VCF(variant_record, sample)
                    )

    @property
    def sample_to_gene_to_VCFs(self):
        return self._sample_to_gene_to_VCFs

    def get_VCF_records_given_sample_and_gene(
        self, sample: str, gene_name: str
    ) -> List[VCF]:
        return self.sample_to_gene_to_VCFs[sample][gene_name]
