from pathlib import Path
from evaluate.mummer import NucmerError
from evaluate.evaluate import generate_mummer_snps
import pytest

cft = Path("tests/test_cases/CFT073.ref.fa")


def test_getMummerSnps_invalidQueryFileRaisesNucmerError():
    reference = cft
    query = Path("foo")
    with pytest.raises(NucmerError):
        result = generate_mummer_snps(reference, query)
