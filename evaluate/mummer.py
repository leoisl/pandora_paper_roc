"""This file holds wrappers for running mummer commands."""
import logging
import subprocess
from pathlib import Path
from typing import List, TextIO, Tuple

import pandas as pd
from .probe import Probe, ProbeHeader, ProbeInterval
from .utils import arg_ranges


class NucmerError(Exception):
    pass


class DeltaFilterError(Exception):
    pass


class ShowSnpsError(Exception):
    pass


class Nucmer:
    def __init__(
        self, reference: Path, query: Path, prefix: str = "out", extra_params: str = ""
    ):
        if not reference.is_file():
            raise NucmerError(f"Reference file {str(reference)} does not exist.")
        elif not query.is_file():
            raise NucmerError(f"Query file {str(query)} does not exist.")
        if not Path(prefix).parent.is_dir():
            raise NucmerError(f"Prefix {str(Path(prefix))} parent does not exist.")

        self.reference = str(reference)
        self.query = str(query)
        self.prefix = prefix
        self.extra_params = extra_params

    def generate_command(self) -> List[str]:
        command = ["nucmer", "--prefix", self.prefix, self.reference, self.query]

        if self.extra_params:
            command.insert(1, self.extra_params)

        return command

    def run(self) -> subprocess.CompletedProcess:
        """The output file is written to <prefix>.delta"""
        nucmer_command = self.generate_command()

        logging.info(f"Running nucmer with command:\n{' '.join(nucmer_command)}")

        return subprocess.run(
            nucmer_command, stderr=subprocess.PIPE, stdout=subprocess.PIPE
        )


class DeltaFilter:
    def __init__(self, deltafile: Path, extra_params: str = ""):
        if not deltafile.is_file():
            raise DeltaFilterError(f"deltafile {str(deltafile)} does not exist.")

        self.deltafile = str(deltafile)
        self.extra_params = extra_params

    def generate_command(self) -> List[str]:
        command = ["delta-filter", self.deltafile]

        if self.extra_params:
            command.insert(1, self.extra_params)

        return command

    def run(self) -> subprocess.CompletedProcess:
        """The output file can be found in the stdout of the returned
        CompletedProcess."""
        deltafilter_command = self.generate_command()

        logging.info(
            f"Running delta-filter with command:\n{' '.join(deltafilter_command)}"
        )

        return subprocess.run(
            deltafilter_command, stderr=subprocess.PIPE, stdout=subprocess.PIPE
        )


class ShowSnps:
    def __init__(
        self,
        deltafile: Path,
        context: int = 0,
        print_header: bool = True,
        indels: bool = True,
        extra_params: str = "",
    ):
        if not deltafile.is_file():
            raise ShowSnpsError(f"deltafile {str(deltafile)} does not exist.")

        self.deltafile = str(deltafile)
        self.context = context
        self.print_header = print_header
        self.indels = indels
        self.extra_params = extra_params

    def generate_command(self) -> List[str]:
        command = ["show-snps", self.deltafile]

        if not self.indels:
            command.insert(1, "-I")

        if not self.print_header:
            command.insert(1, "-H")

        if self.context > 0:
            command.insert(1, f"-x {self.context}")

        if self.extra_params:
            command.insert(1, self.extra_params)

        return command

    def run(self) -> subprocess.CompletedProcess:
        """The output file can be found in the stdout of the returned
        CompletedProcess."""
        showsnps_command = self.generate_command()

        logging.info(f"Running show-snps with command:\n{' '.join(showsnps_command)}")

        return subprocess.run(
            showsnps_command, stdout=subprocess.PIPE, stderr=subprocess.PIPE
        )

    @staticmethod
    def to_dataframe(snps: TextIO) -> "ShowSNPsDataframe":
        """Note: this method is not general. i.e it is only setup at the moment to
        parse a show-snps file where the options used were -rlTC and -x"""
        cols = {
            "ref_pos": int,  #  P1,
            "ref_sub": str,  #  SUB,
            "query_sub": str,  #  SUB
            "query_pos": int,  #  P2
            "nearest_mismatch": int,  #  BUFF
            "nearest_end": int,  #  DIST
            "ref_len": int,  #  LEN R
            "query_len": int,  #  LEN Q
            "ref_context": str,  # CTX R
            "query_context": str,  #  CTX Q
            "ref_strand": int,
            "query_strand": int,
            "ref_chrom": str,
            "query_chrom": str,
        }
        names = list(cols.keys())
        return ShowSNPsDataframe(
            pd.read_csv(
                snps, sep="\t", skiprows=4, index_col=False, names=names, dtype=cols
            )
        )


class ShowSNPsDataframe(pd.DataFrame):
    @property
    def _constructor(self):
        return ShowSNPsDataframe

    def translate_to_FWD_strand(self) -> "ShowSNPsDataframe":
        def fix_position(position: int, strand_aln: int, length: int) -> int:
            if strand_aln == 1:
                return position
            else:
                return length - position + 1

        def translate_to_FWD_strand_core(line: pd.Series) -> pd.Series:
            line.ref_pos = fix_position(line.ref_pos, line.ref_strand, line.ref_len)
            line.ref_strand = 1
            line.query_pos = fix_position(
                line.query_pos, line.query_strand, line.query_len
            )
            line.query_strand = 1
            return line

        return self.apply(translate_to_FWD_strand_core, axis=1)

    def get_probes(self, id_prefix: str = "ID") -> Tuple[str, str]:
        ref_probes = []
        query_probes = []
        merged_indices = arg_ranges(self.ref_pos.tolist())

        for interval_index, interval in enumerate(merged_indices):
            ref_probe, query_probe = self.probes_for_interval(interval, variation_id = f"{id_prefix}_{interval_index}")
            ref_probes.append(str(ref_probe))
            query_probes.append(str(query_probe))

        return (
            "\n".join(probe for probe in ref_probes),
            "\n".join(probe for probe in query_probes),
        )

    def probes_for_interval(self, interval: Tuple[int, int], variation_id: str) -> Tuple[Probe, ...]:
        probes = []
        probe_prefixes = ["ref", "query"]
        consecutive_positions = self.iloc[slice(*interval)]
        first_row = consecutive_positions.iloc[0]
        flank_width = int((len(first_row[f"{probe_prefixes[0]}_context"]) - 1) / 2)

        for allele, prefix in enumerate(probe_prefixes):
            core_sequence = "".join(consecutive_positions[f"{prefix}_sub"]).replace(
                ".", ""
            )
            left_flank = first_row[f"{prefix}_context"][:flank_width].replace("-", "")
            right_flank = consecutive_positions.iloc[-1][f"{prefix}_context"][
                flank_width + 1 :
            ].replace("-", "")
            call_start_idx = len(left_flank)
            call_end_idx = call_start_idx + len(core_sequence)
            header = ProbeHeader(
                chrom=first_row[f"{prefix}_chrom"],
                pos=first_row[f"{prefix}_pos"],
                interval=ProbeInterval(call_start_idx, call_end_idx),
                variation_id=f"{variation_id}_allele_{allele}"
            )
            full_sequence = left_flank + core_sequence + right_flank
            probes.append(Probe(header=header, full_sequence=full_sequence))

        return tuple(probes)

    def make_pos_zero_based(self) -> "ShowSNPsDataframe":
        df = self.copy(deep=True)
        if df.empty:
            return df
        df["ref_pos"] = df["ref_pos"] - 1
        df["query_pos"] = df["query_pos"] - 1
        return df
