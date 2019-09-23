import pysam


class VCF:
    def __init__(self, variant: pysam.VariantRecord = None, sample: str = None):
        self.variant = variant
        self.sample = sample

    @property
    def genotype(self) -> int:
        return self.variant.samples[self.sample]["GT"][0]

    @property
    def is_invalid_vcf_entry(self) -> bool:
        genotype = self.genotype
        return genotype is None

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
        return self.variant.pos

    @property
    def start(self) -> int:
        return self.variant.start

    @property
    def rlen(self) -> int:
        return self.variant.rlen

    @property
    def chrom(self) -> str:
        return self.variant.chrom
