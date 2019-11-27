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

        # TODO : I don't think we should worry about checking if a VCF is bugged or not
        # TODO : in principle, we should not receive bugged VCFs...
        # if vcf.has_genotype_bug:
        #     raise BuggedVCFError()

        return vcf

    def __eq__(self, other):
        return self.__dict__ == other.__dict__

    @property
    def genotype(self) -> int:
        data_from_sample = self.variant.samples[self.sample]
        all_gts = data_from_sample.get("GT")
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
        gt_conf = data_from_sample.get("GT_CONF")
        if gt_conf is not None:
            return float(gt_conf)
        else:
            # TODO: this is required due to Snippy, put this in a hierarchy (since we will need to do this for nanopolish/medaka also)
            return float(self.variant.qual)

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
        # TODO: this is required due to Snippy, put this in a hierarchy (since we will need to do this for nanopolish/medaka also)
        try:
            return self.variant.info["SVTYPE"]
        except KeyError:
            return self.variant.info["TYPE"]

    @property
    def mean_coverage_forward(self) -> int:
        genotype = self.genotype
        # TODO: this is required due to Snippy, put this in a hierarchy (since we will need to do this for nanopolish/medaka also)
        try:
            return int(self.variant.samples[self.sample]["MEAN_FWD_COVG"][genotype])
        except KeyError:
            return 0

    @property
    def mean_coverage_reverse(self) -> int:
        genotype = self.genotype
        # TODO: this is required due to Snippy, put this in a hierarchy (since we will need to do this for nanopolish/medaka also)
        try:
            return int(self.variant.samples[self.sample]["MEAN_REV_COVG"][genotype])
        except KeyError:
            return 0

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
