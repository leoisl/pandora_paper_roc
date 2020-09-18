import pysam
from typing import List, Dict, TextIO
from .vcf import VCF, NullVCFError, VCFFactory
from collections import defaultdict


class VCFFile:
    def __init__(self, pysam_variant_file: pysam.VariantFile, VCF_creator_method):
        self._header = pysam_variant_file.header
        self._sample_to_gene_to_VCFs = defaultdict(lambda: defaultdict(list))
        for variant_record in pysam_variant_file:
            for sample in variant_record.samples:
                gene = variant_record.chrom
                try:
                    vcf = VCF_creator_method(variant_record, sample)
                    self._sample_to_gene_to_VCFs[sample][gene].append(vcf)
                except NullVCFError:
                    pass

    @property
    def sample_to_gene_to_VCFs(self) -> Dict[str, Dict[str, List[VCF]]]:
        return self._sample_to_gene_to_VCFs

    @property
    def header(self) -> str:
        return str(self._header)

    def get_VCF_records_given_sample_and_gene(
        self, sample: str, gene_name: str
    ) -> List[VCF]:
        return self.sample_to_gene_to_VCFs[sample][gene_name]

    def write(self, filehandler: TextIO):
        filehandler.write(self.header)
        for (sample, gene_to_VCFs_of_a_sample) in self.sample_to_gene_to_VCFs.items():
            for gene, vcfs in gene_to_VCFs_of_a_sample.items():
                for vcf in vcfs:
                    filehandler.write(str(vcf))
