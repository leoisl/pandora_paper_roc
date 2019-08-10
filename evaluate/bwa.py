import subprocess
import pysam


class BWA:
    def __init__(self, threads=1):
        self.threads = threads
        self.reference = ""

    def index(self, reference: str):
        self.reference = reference

        completed_process = subprocess.run(
            ["bwa", "index", self.reference],
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
        )
        completed_process.check_returncode()

    def align(self, query: str) -> tuple:
        options = self.get_options()
        self.alignment = subprocess.run(
            ["bwa", "mem", *options, str(self.reference), "-"],
            input=query.encode(),
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
        )

        if self.alignment.returncode != 0:
            if b"fail to locate the index" in self.alignment.stderr:
                raise IndexError("Reference must be indexed by BWA before alignment.")
            else:
                self.alignment.check_returncode()

        return self._get_samfile()

    def get_options(self):
        options = []
        options.extend(["-t", str(self.threads)])

        return options

    def _get_samfile(self):
        header = ""
        sam_lines = []
        for line in self.alignment.stdout.decode().split("\n"):
            if line.startswith("@"):
                header += line + "\n"
            else:
                sam_lines.append(line)

        header = pysam.AlignmentHeader.from_text(header)

        return (
            header,
            [pysam.AlignedSegment.fromstring(sam, header) for sam in sam_lines if sam],
        )
