"""This file holds wrappers for running mummer commands."""
from pathlib import Path
import logging
import subprocess
from typing import List


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
