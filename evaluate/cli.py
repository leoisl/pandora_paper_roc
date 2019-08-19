"""Assess the recall of a VCF with respect to two samples.
"""
import argparse
import logging
import sys
from pathlib import Path

GT_MAX = 300
GT_STEP = 5
TRUTH_FLANK_WIDTH = 21
QUERY_FLANK_WIDTH = 45
TEMP_DIR = Path.cwd() / "tmp"
LOGGING_LEVELS = {
    0: "NOTSET",
    1: "CRITICAL",
    2: "ERROR",
    3: "WARNING",
    4: "INFO",
    5: "DEBUG",
}
DEFAULT_LOG_LVL = 4


def cli() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "-1",
        "--query1",
        help="One of two assembly/truth sequences to compare.",
        required=True,
        type=Path,
    )
    parser.add_argument(
        "-2",
        "--query2",
        help="One of two assembly/truth sequences to compare.",
        required=True,
        type=Path,
    )
    parser.add_argument(
        "-v", "--vcf", required=True, help="""Path to the VCF to evaluate.""", type=Path
    )
    parser.add_argument(
        "-r",
        "--vcf-ref",
        required=True,
        help="""Path to the sequence(s) that the VCF is with respect to.""",
        type=Path,
    )
    parser.add_argument(
        "-o",
        "--output",
        type=str,
        help="""Path to output results to. Default is STDOUT [-].""",
        default="-",
    )
    parser.add_argument(
        "-t",
        "--max-gt-threshold",
        type=int,
        help=f"Maximum genotype confidence threshold [{GT_MAX}]",
        default=GT_MAX,
    )
    parser.add_argument(
        "-s",
        "--gt-step",
        type=int,
        help=f"Steps between successive genotype confidence thresholds [{GT_STEP}]",
        default=GT_STEP,
    )
    parser.add_argument(
        "--truth-flank",
        type=int,
        help=(
            "Width of flank to add either side of the truth probes "
            f"[{TRUTH_FLANK_WIDTH}]"
        ),
        default=TRUTH_FLANK_WIDTH,
    )
    parser.add_argument(
        "--query-flank",
        type=int,
        help=(
            "Width of flank to add either side of the query/VCF probes "
            f"[{QUERY_FLANK_WIDTH}]"
        ),
        default=QUERY_FLANK_WIDTH,
    )
    parser.add_argument("--indels", help="Include indels.", action="store_true")
    parser.add_argument(
        "--temp",
        help=f"Location for storing temporary files [{TEMP_DIR}]",
        type=Path,
        default=TEMP_DIR,
    )
    parser.add_argument(
        "--threads", help="Number of threads to use for BWA. [1]", type=int, default=1
    )
    parser.add_argument(
        "--log_level",
        help=(
            "Level of logging. 0 is none, 5 is for debugging. "
            f"Default is {LOGGING_LEVELS[DEFAULT_LOG_LVL]} [{DEFAULT_LOG_LVL}] "
            "which will report info, warnings, errors, and critical information."
        ),
        default=DEFAULT_LOG_LVL,
        type=int,
        choices=range(6),
    )
    args = parser.parse_args()

    log_level = LOGGING_LEVELS[args.log_level]
    logging.basicConfig(
        level=log_level,
        format="[%(asctime)s]:%(levelname)s: %(message)s",
        datefmt="%m/%d/%Y %I:%M:%S %p",
    )

    if args.output == "-":
        args.output = sys.stdout
    else:
        p = Path(args.output)
        if not p.parent.is_dir():
            raise NotADirectoryError(
                "Directory specified for output file does not exist: {}".format(
                    p.parent
                )
            )
        args.output = p.open("w")

    args.temp.mkdir(exist_ok=True)

    logging.debug(args)

    return args
