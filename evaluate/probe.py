import re
from typing import NamedTuple, Type


class RegexError(Exception):
    pass


class Interval(NamedTuple):
    start: int = 0
    end: int = 0

    @staticmethod
    def from_string(string: str) -> "Interval":
        regex = re.compile(r"\[(\d+),(\d+)\)")
        match = regex.search(string)
        if not match:
            raise RegexError(f"Could not parse interval string: {string}")

        return Interval(int(match.group(1)), int(match.group(2)))


class Probe:
    def __init__(
        self,
        sample: str = "",
        chrom: str = "",
        pos: int = 0,
        interval: Interval = Interval(),
        svtype: str = "",
        mean_fwd_covg: int = 0,
        mean_rev_covg: int = 0,
        gt_conf: float = 0,
        full_sequence: str = "",
    ):
        self.chrom = chrom
        self.sample = sample
        self.pos = pos
        self.interval = interval
        self.svtype = svtype
        self.mean_fwd_covg = mean_fwd_covg
        self.mean_rev_covg = mean_rev_covg
        self.gt_conf = gt_conf
        self.full_sequence = full_sequence

    def __eq__(self, other: "Probe"):
        return (
            self.chrom == other.chrom
            and self.sample == other.sample
            and self.pos == other.pos
            and self.interval == other.interval
            and self.svtype == other.svtype
            and self.mean_fwd_covg == other.mean_fwd_covg
            and self.mean_rev_covg == other.mean_rev_covg
            and self.gt_conf == other.gt_conf
            and self.full_sequence == other.full_sequence
        )

    @staticmethod
    def from_string(string: str) -> "Probe":
        if not string:
            return Probe()

        fields = string.split("\n")
        header = fields[0]
        probe = Probe.from_header(header)

        if len(fields) > 1 and fields[1]:
            probe.full_sequence = fields[1].rstrip()

        return probe

    @staticmethod
    def from_header(string: str) -> "Probe":
        def parse_field_from_header(
            field: str, header: str, return_type: Type=str, delim: str = ";"
        ):
            regex = re.compile(f"{field}=(.+?){delim}")
            match = regex.search(header)
            if match:
                return return_type(match.group(1))
            else:
                return return_type()

        chrom = parse_field_from_header("CHROM", string)
        sample = parse_field_from_header("SAMPLE", string)
        pos = parse_field_from_header("POS", string, return_type=int)
        svtype = parse_field_from_header("SVTYPE", string)
        mean_fwd_covg = parse_field_from_header("MEAN_FWD_COVG", string, return_type=int)
        mean_rev_covg = parse_field_from_header("MEAN_REV_COVG", string, return_type=int)
        gt_conf = parse_field_from_header("GT_CONF", string, return_type=float)
        interval = Interval.from_string(parse_field_from_header("INTERVAL", string))

        return Probe(
            sample, chrom, pos, interval, svtype, mean_fwd_covg, mean_rev_covg, gt_conf
        )
