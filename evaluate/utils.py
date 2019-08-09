from pathlib import Path


def strip_extensions(path: Path) -> Path:
    return Path(str(path).replace("".join(path.suffixes), ""))
