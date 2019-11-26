import pysam
from typing import Iterable

class BuggedVCFError(Exception):
    pass

class NullVCFError(Exception):
    pass

class VCF:
    @staticmethod
    def from_VariantRecord_and_Sample(variant: pysam.VariantRecord = None, sample: str = None) -> "VCF":
        vcf = VCF()
        vcf.variant = variant
        vcf.sample = sample

        if vcf.is_null_call:
            raise NullVCFError()

        if vcf.has_genotype_bug:
            raise BuggedVCFError()

        return vcf

    def __eq__(self, other):
        return self.__dict__ == other.__dict__

    @property
    def genotype(self) -> int:
        data_from_sample = self.variant.samples[self.sample]
        all_gts = data_from_sample.get("GT")
        assert len(all_gts) == 1
        return all_gts[0]

    @property
    def is_null_call(self) -> bool:
        return self.genotype is None

    @property
    def has_genotype_bug(self) -> bool:
        genotype = self.genotype
        genotype_called_wrongly = genotype not in self.highest_likelihood_indexes
        return genotype_called_wrongly

    @property
    def genotype_confidence(self) -> float:
        data_from_sample = self.variant.samples[self.sample]
        return float(data_from_sample.get("GT_CONF"))

    @property
    def called_variant_sequence(self) -> str:
        genotype = self.genotype
        if genotype == 0:
            return self.variant.ref
        else:
            return self.variant.alleles[genotype]

    @property
    def called_variant_length(self) -> int:
        return len(self.called_variant_sequence)

    @property
    def svtype(self) -> str:
        return self.variant.info["SVTYPE"]

    @property
    def mean_coverage_forward(self) -> int:
        genotype = self.genotype
        return int(self.variant.samples[self.sample]["MEAN_FWD_COVG"][genotype])

    @property
    def mean_coverage_reverse(self) -> int:
        genotype = self.genotype
        return int(self.variant.samples[self.sample]["MEAN_REV_COVG"][genotype])

    @property
    def mean_coverage(self) -> int:
        return self.mean_coverage_forward + self.mean_coverage_reverse

    @property
    def pos(self) -> int:
        return int(self.variant.pos)

    @property
    def start(self) -> int:
        return int(self.variant.start)

    @property
    def stop(self) -> int:
        return int(self.variant.stop)

    @property
    def rlen(self) -> int:
        return int(self.variant.rlen)

    @property
    def chrom(self) -> str:
        return self.variant.chrom

    @property
    def likelihoods(self) -> Iterable[float]:
        return [
            float(likelihood)
            for likelihood in self.variant.samples[self.sample]["LIKELIHOOD"]
        ]

    @property
    def highest_likelihood_indexes(self) -> Iterable[int]:
        highest_likelihood_indexes = [
            index
            for index, likelihood in enumerate(self.likelihoods)
            if likelihood == max(self.likelihoods)
        ]
        return highest_likelihood_indexes

    @property
    def gaps(self) -> float:
        return float(self.variant.samples[self.sample]["GAPS"][self.genotype])
