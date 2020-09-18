import pysam
from typing import Iterable, List
from abc import ABC, abstractmethod

class NullVCFError(Exception):
    pass


class VCF(ABC):
    def __init__(self, variant: pysam.VariantRecord = None, sample: str = None):
        self.variant = variant
        self.sample = sample

        if self.variant is not None and self.sample is not None:
            if self.is_null_call:
                raise NullVCFError()

    @property
    def is_null_call(self) -> bool:
        return self.genotype is None

    @property
    @abstractmethod
    def genotype(self) -> int:
        pass


    @property
    @abstractmethod
    def genotype_confidence(self) -> float:
        pass

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
    @abstractmethod
    def svtype(self) -> str:
        pass


    @property
    @abstractmethod
    def coverage(self) -> int:
        pass

    @property
    def pos(self) -> int:
        return int(self.variant.pos)

    @property
    def ref(self) -> str:
        return self.variant.ref

    @property
    def ref_length(self) -> int:
        return len(self.ref)

    @property
    def start(self) -> int:
        return int(self.variant.start)

    @property
    def stop(self) -> int:
        return int(self.variant.stop)

    def positions_covered_by_variant_with_flanks(self, flank_length) -> List[int]:
        return list(range(self.start-flank_length, self.stop+flank_length))

    @property
    def rlen(self) -> int:
        return int(self.variant.rlen)

    @property
    def chrom(self) -> str:
        return self.variant.chrom

    def __str__(self) -> str:
        return str(self.variant)




class PandoraVCF(VCF):
    def __eq__(self, other):
        return self.__dict__ == other.__dict__

    ####################################################################################################################
    # Overriding methods
    @property
    def genotype(self) -> int:
        data_from_sample = self.variant.samples[self.sample]
        all_gts = data_from_sample.get("GT")
        assert len(all_gts) == 1
        return all_gts[0]


    @property
    def genotype_confidence(self) -> float:
        data_from_sample = self.variant.samples[self.sample]
        return float(data_from_sample.get("GT_CONF"))

    @property
    def svtype(self) -> str:
        return self.variant.info["SVTYPE"]

    @property
    def coverage(self) -> int:
        return self._mean_coverage_forward + self._mean_coverage_reverse

    ####################################################################################################################


    @property
    def _mean_coverage_forward(self) -> int:
        genotype = self.genotype
        return int(self.variant.samples[self.sample]["MEAN_FWD_COVG"][genotype])

    @property
    def _mean_coverage_reverse(self) -> int:
        genotype = self.genotype
        return int(self.variant.samples[self.sample]["MEAN_REV_COVG"][genotype])

    @property
    def _likelihoods(self) -> Iterable[float]:
        return [
            float(likelihood)
            for likelihood in self.variant.samples[self.sample]["LIKELIHOOD"]
        ]

    @property
    def _highest_likelihood_indexes(self) -> Iterable[int]:
        highest_likelihood_indexes = [
            index
            for index, likelihood in enumerate(self._likelihoods)
            if likelihood == max(self._likelihoods)
        ]
        return highest_likelihood_indexes

    @property
    def _gaps(self) -> float:
        return float(self.variant.samples[self.sample]["GAPS"][self.genotype])




class SnippyVCF(VCF):
    def __eq__(self, other):
        return self.__dict__ == other.__dict__

    ####################################################################################################################
    # Overriding methods
    @property
    def genotype(self) -> int:
        data_from_sample = self.variant.samples[self.sample]
        all_gts = data_from_sample.get("GT")
        assert all_gts == (1, 1)
        return 1


    @property
    def genotype_confidence(self) -> float:
        return float(self.variant.qual)

    @property
    def svtype(self) -> str:
        return self.variant.info["TYPE"][0]

    @property
    def coverage(self) -> int:
        # return int(self.variant.info["AO"][0])
        return 1000  # mocking snippy coverage, as we don't want to filter it by coverage

    ####################################################################################################################


class SamtoolsVCF(VCF):
    def __eq__(self, other):
        return self.__dict__ == other.__dict__

    ####################################################################################################################
    # Overriding methods
    @property
    def genotype(self) -> int:
        data_from_sample = self.variant.samples[self.sample]
        GT = data_from_sample.get("GT")
        assert GT == (1, )
        return 1


    @property
    def genotype_confidence(self) -> float:
        return float(self.variant.qual)

    @property
    def svtype(self) -> str:
        if "INDEL" in self.variant.info:
            return "INDEL"
        else:
            return "SNP"

    @property
    def coverage(self) -> int:
        return int(self.variant.info["DP"])

    ####################################################################################################################


class MedakaVCF(VCF):
    def __eq__(self, other):
        return self.__dict__ == other.__dict__

    ####################################################################################################################
    # Overriding methods
    @property
    def genotype(self) -> int:
        data_from_sample = self.variant.samples[self.sample]
        GT = data_from_sample.get("GT")
        assert GT == (1,)
        return 1

    @property
    def genotype_confidence(self) -> float:
        return float(self.variant.qual)

    @property
    def svtype(self) -> str:
        return "NA"

    @property
    def coverage(self) -> int:
        return 1000 # unknown

    ####################################################################################################################


class NanopolishVCF(VCF):
    def __eq__(self, other):
        return self.__dict__ == other.__dict__

    ####################################################################################################################
    # Overriding methods
    @property
    def genotype(self) -> int:
        data_from_sample = self.variant.samples[self.sample]
        GT = data_from_sample.get("GT")
        assert GT == (1,)
        return 1

    @property
    def genotype_confidence(self) -> float:
        return float(self.variant.qual)

    @property
    def svtype(self) -> str:
        return "NA"

    @property
    def coverage(self) -> int:
        return int(self.variant.info["BaseCalledReadsWithVariant"])

    ####################################################################################################################


class VCFFactory:
    @staticmethod
    def create_Pandora_VCF_from_VariantRecord_and_Sample(variant: pysam.VariantRecord = None, sample: str = None) -> PandoraVCF:
        vcf = PandoraVCF(variant, sample)
        return vcf

    @staticmethod
    def create_Snippy_VCF_from_VariantRecord_and_Sample(variant: pysam.VariantRecord = None, sample: str = None) -> SnippyVCF:
        vcf = SnippyVCF(variant, sample)
        return vcf

    @staticmethod
    def create_Samtools_VCF_from_VariantRecord_and_Sample(variant: pysam.VariantRecord = None,
                                                          sample: str = None) -> SamtoolsVCF:
        vcf = SamtoolsVCF(variant, sample)
        return vcf

    @staticmethod
    def create_Medaka_VCF_from_VariantRecord_and_Sample(variant: pysam.VariantRecord = None,
                                                          sample: str = None) -> MedakaVCF:
        vcf = MedakaVCF(variant, sample)
        return vcf


    @staticmethod
    def create_Nanopolish_VCF_from_VariantRecord_and_Sample(variant: pysam.VariantRecord = None,
                                                          sample: str = None) -> NanopolishVCF:
        vcf = NanopolishVCF(variant, sample)
        return vcf
