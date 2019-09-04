from pathlib import Path
from typing import Tuple, List

import pysam

from .bwa import BWA


def map_probes_to_truth(
    probes: Path, truth: Path, output: Path = Path(), threads: int = 1
) -> Tuple[pysam.VariantHeader, List[pysam.AlignedSegment]]:
    return BWA.map_query_to_ref(probes, truth, output, threads)
