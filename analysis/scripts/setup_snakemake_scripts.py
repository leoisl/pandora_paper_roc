def setup_snakemake_scripts(log, log_level="DEBUG"):
    from pathlib import Path
    import sys

    sys.path.append(str(Path().absolute()))
    import logging

    logging.basicConfig(
        filename=log,
        filemode="w",
        level=log_level,
        format="[%(asctime)s]:%(levelname)s: %(message)s",
        datefmt="%m/%d/%Y %I:%M:%S %p",
    )