import re
from typing import NamedTuple, Type

DELIM = ";"


class RegexError(Exception):
    pass


class ProbeInterval(NamedTuple):
    start: int = -1
    end: int = -1

    def __bool__(self) -> bool:
        return not self.is_null()

    def __str__(self) -> str:
        if self.is_null():
            return ""

        return f"[{self.start},{self.end})"

    def __len__(self) -> int:
        return self.end - self.start

    def is_null(self) -> bool:
        return self.start == -1 and self.end == -1

    interval_regex = re.compile(r"\[(\d+),(\d+)\)")

    @staticmethod
    def from_string(string: str) -> "ProbeInterval":
        match = ProbeInterval.interval_regex.search(string)
        if not match:
            return ProbeInterval()

        return ProbeInterval(int(match.group(1)), int(match.group(2)))


class ProbeHeader:
    def __init__(
        self,
        sample: str = "",
        chrom: str = "",
        pos: int = 0,
        interval: ProbeInterval = ProbeInterval(),
        svtype: str = "",
        gt_conf: float = None,
    ):
        self.chrom = chrom
        self.sample = sample
        self.pos = pos
        self.interval = interval
        self.svtype = svtype
        self.gt_conf = gt_conf

    def __eq__(self, other: "ProbeHeader") -> bool:
        return (
            self.chrom == other.chrom
            and self.sample == other.sample
            and self.pos == other.pos
            and self.interval == other.interval
            and self.svtype == other.svtype
            and self.gt_conf == other.gt_conf
        )

    def __str__(self) -> str:
        list_of_key_values_to_add = [f"{k.upper()}={str(v)}" for k, v in vars(self).items() if k != "gt_conf" and v]
        # forcibly add gt conf
        if self.gt_conf is not None:
            list_of_key_values_to_add.append(f"GT_CONF={round(self.gt_conf, 1):.1f}")

        contents = DELIM.join(list_of_key_values_to_add)

        if not contents:
            return ""

        return f">{contents}{DELIM}"

    @staticmethod
    def from_string(string: str) -> "ProbeHeader":
        def parse_field_from_header(
            field: str, header: str, return_type: Type, value_to_return_if_not_found, delim: str = DELIM,
        ):
            regex = re.compile(f"{field}=(.+?){delim}")
            match = regex.search(header)
            if match:
                return return_type(match.group(1))
            else:
                return value_to_return_if_not_found

        chrom = parse_field_from_header("CHROM", string, str, "")
        sample = parse_field_from_header("SAMPLE", string, str, "")
        pos = parse_field_from_header("POS", string, int, 0)
        svtype = parse_field_from_header("SVTYPE", string, str, "")
        gt_conf = parse_field_from_header("GT_CONF", string, float, None)
        interval = ProbeInterval.from_string(
            parse_field_from_header("INTERVAL", string, str, "")
        )

        return ProbeHeader(
            sample, chrom, pos, interval, svtype, gt_conf
        )


class Probe:
    def __init__(self, header: ProbeHeader = ProbeHeader(), full_sequence: str = ""):
        self.header = header
        self.full_sequence = full_sequence

    def __eq__(self, other: "Probe"):
        return self.header == other.header and self.full_sequence == other.full_sequence

    def __str__(self) -> str:
        header = str(self.header)

        if not header:
            return ""

        return f"{header}\n{self.full_sequence}"

    @property
    def left_flank(self) -> str:
        end = self.header.interval.start

        return self.full_sequence[:end]

    @property
    def right_flank(self) -> str:
        start = self.header.interval.end

        return self.full_sequence[start:]

    @property
    def core_sequence(self) -> str:
        return self.full_sequence[slice(*self.interval)]

    @property
    def interval(self) -> ProbeInterval:
        return self.header.interval

    @property
    def gt_conf(self) -> float:
        return self.header.gt_conf

    @property
    def chrom(self) -> str:
        return self.header.chrom

    @property
    def pos(self) -> int:
        return self.header.pos

    @property
    def is_deletion(self) -> bool:
        return len(self.interval) == 0

    @staticmethod
    def from_string(string: str) -> "Probe":
        if not string:
            return Probe()

        fields = string.split("\n")
        header = ProbeHeader.from_string(fields[0])
        probe = Probe(header=header)

        if len(fields) > 1 and fields[1]:
            probe.full_sequence = fields[1].rstrip()

        return probe
