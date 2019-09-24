import pysam
from typing import Iterable


class VCF:
    def __init__(self, variant: pysam.VariantRecord = None, sample: str = None):
        self.variant = variant
        self.sample = sample

    def __eq__(self, other):
        return (
            self.genotype,
            self.genotype_confidence,
            self.variant_sequence,
            self.svtype,
            self.mean_coverage_forward,
            self.mean_coverage_reverse,
            self.pos,
            self.chrom,
        ) == (
            other.genotype,
            other.genotype_confidence,
            other.variant_sequence,
            other.svtype,
            other.mean_coverage_forward,
            other.mean_coverage_reverse,
            other.pos,
            other.chrom,
        )

    @property
    def genotype(self) -> int:
        return self.variant.samples[self.sample]["GT"][0]

    @property
    def is_invalid_vcf_entry(self) -> bool:
        genotype = self.genotype
        no_genotype_called = genotype is None
        if no_genotype_called:
            return True

        genotype_called_wrongly = genotype not in self.highest_likelihood_indexes
        if genotype_called_wrongly:
            return True

        return False

    @property
    def genotype_confidence(self) -> float:
        data_from_sample = self.variant.samples[self.sample]
        return float(data_from_sample.get("GT_CONF", 0))

    @property
    def variant_sequence(self) -> str:
        genotype = self.genotype

        if genotype is None:
            return self.variant.ref
        else:
            return self.variant.alleles[genotype]

    @property
    def variant_length(self) -> int:
        return len(self.variant_sequence)

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
    def rlen(self) -> int:
        return int(self.variant.rlen)

    @property
    def chrom(self) -> str:
        return self.variant.chrom

    @property
    def likelihoods(self) -> Iterable[float]:
        return [float(likelihood) for likelihood in self.variant.samples[self.sample]["LIKELIHOOD"]]

    @property
    def highest_likelihood_indexes(self) -> Iterable[int]:
        highest_likelihood_indexes = [index for index, likelihood in enumerate(self.likelihoods) if
                                      likelihood == max(self.likelihoods)]
        return highest_likelihood_indexes
