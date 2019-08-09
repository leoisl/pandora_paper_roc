from pathlib import Path
from evaluate.utils import *

def test_stripExtensions_noExtensionsReturnsSameAsInput():
    path = Path("path/to/foo")

    actual = strip_extensions(path)
    expected = path

    assert actual == expected


def test_stripExtensions_oneExtensionReturnsInputWithoutExtension():
    path = Path("path/to/foo.fa")

    actual = strip_extensions(path)
    expected = Path("path/to/foo")

    assert actual == expected

def test_stripExtensions_twoExtensionsReturnsInputWithoutAnyExtensions():
    path = Path("path/to/foo.fa.gz")

    actual = strip_extensions(path)
    expected = Path("path/to/foo")

    assert actual == expected