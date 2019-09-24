from pathlib import Path
import pysam
from typing import List, Iterable
from .vcf import VCF

class VCFFile:
    def __init__(self, pysam_variant_file: pysam.VariantFile):
        self._pysam_variant_file = pysam_variant_file
        self.sample_to_gene_to_VCF = {}
        with op

    # @property
    # def pysam_variant_file(self):
    #     return self._pysam_variant_file

    def subset_samples(self, samples: List[str]):
        self.pysam_variant_file.subset_samples(samples)

    def fetch(self, gene_name: str) -> Iterable:
        return self.pysam_variant_file.fetch(contig=gene_name)

    def get_VCF_records_given_sample_and_gene(self, sample: str, gene_name: str) -> List[VCF]:
        self.subset_samples([sample])
        variants_iterator = self.fetch(gene_name)
        vcf_records = [VCF(variant, sample) for variant in variants_iterator]
        return vcf_records


    @staticmethod
    def from_path(filepath : Path) -> "VCFFile":
        indexed_vcf = pysam.VariantFile(Path(
            pysam.tabix_index(str(filepath), preset="vcf", keep_original=True, force=True)
        ))
        return VCFFile(indexed_vcf)