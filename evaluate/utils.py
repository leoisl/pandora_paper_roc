from pathlib import Path
from typing import List, Tuple


def strip_extensions(path: Path) -> Path:
    return Path(str(path).replace("".join(path.suffixes), ""))


def arg_ranges(array: List[int]) -> List[Tuple[int, int]]:
    if not array:
        return []

    argranges = []
    index_cache = [0]
    previous = array[0]

    for index in range(1, len(array)):
        current = array[index]
        if consecutive(previous, current):
            index_cache.append(index)
        else:
            argranges.append(index_cache)
            index_cache = [index]
        previous = current

    argranges.append(index_cache)

    return collapse_ranges(argranges)


def consecutive(x: int, y: int) -> bool:
    diff = x - y
    return -1 <= diff <= 1


def collapse_ranges(ranges: List[List[int]]) -> List[Tuple[int, int]]:
    return [(xs[0], xs[-1] + 1) for xs in ranges if xs]
