#!/usr/bin/env python3
"""
bigwig_compare.py
-----------------
Wrapper around deeptools bigwigCompare.

Examples
--------
# Log2 ratio
python bigwig_compare.py \
    -b1 treatment.bw \
    -b2 control.bw \
    -o treatment_vs_control.bw \
    --operation log2

# Subtract
python bigwig_compare.py \
    -b1 sample1.bw \
    -b2 sample2.bw \
    -o difference.bw \
    --operation subtract

# Mean
python bigwig_compare.py \
    -b1 rep1.bw \
    -b2 rep2.bw \
    -o mean.bw \
    --operation mean
"""

import argparse
import subprocess
import sys
import shutil
from pathlib import Path


def check_deeptools() -> None:
    """Verify that deeptools (bigwigCompare) is available on PATH."""
    if shutil.which("bigwigCompare") is None:
        sys.exit(
            "[ERROR] 'bigwigCompare' not found on PATH.\n"
            "Install deeptools with:\n"
            "  conda install -c bioconda deeptools"
        )


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Compare two BigWig files using deeptools bigwigCompare.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )

    parser.add_argument(
        "-b1", "--bigwig1",
        required=True,
        help="First BigWig file."
    )

    parser.add_argument(
        "-b2", "--bigwig2",
        required=True,
        help="Second BigWig file."
    )

    parser.add_argument(
        "-o", "--output",
        required=True,
        help="Output file."
    )

    parser.add_argument(
        "--operation",
        default="log2",
        choices=[
            "log2",
            "ratio",
            "subtract",
            "add",
            "mean",
            "reciprocal_ratio",
            "first",
            "second"
        ],
        help="Operation to perform."
    )

    parser.add_argument(
        "--binSize",
        type=int,
        default=10,
        help="Bin size in bp."
    )

    parser.add_argument(
        "--numberOfProcessors", "-p",
        type=int,
        default=4,
        help="Number of processors."
    )

    parser.add_argument(
        "--pseudocount",
        nargs=2,
        type=float,
        default=[1.0, 1.0],
        metavar=("PC1", "PC2"),
        help="Pseudocounts added before ratio/log2 calculations."
    )

    parser.add_argument(
        "--skipNonCoveredRegions",
        action="store_true",
        help="Skip regions without coverage."
    )

    parser.add_argument(
        "--blackListFileName",
        default=None,
        help="Blacklist BED file."
    )

    parser.add_argument(
        "--skip-check",
        action="store_true",
        help="Skip deeptools availability check."
    )

    parser.add_argument(
        "--verbose", "-v",
        action="store_true",
        help="Print command before execution."
    )

    return parser.parse_args()


def validate_inputs(files):
    missing = [f for f in files if not Path(f).is_file()]
    if missing:
        sys.exit(
            "[ERROR] Missing input file(s):\n" +
            "\n".join(missing)
        )


def build_command(args):

    cmd = [
        "bigwigCompare",
        "--bigwig1", str(Path(args.bigwig1).resolve()),
        "--bigwig2", str(Path(args.bigwig2).resolve()),
        "--operation", args.operation,
        "--outFileName", str(Path(args.output).resolve()),
        "--binSize", str(args.binSize),
        "--numberOfProcessors", str(args.numberOfProcessors),
    ]

    if args.operation in ["log2", "ratio", "reciprocal_ratio"]:
        cmd.extend(
            [
                "--pseudocount",
                str(args.pseudocount[0]),
                str(args.pseudocount[1])
            ]
        )

    if args.skipNonCoveredRegions:
        cmd.append("--skipNonCoveredRegions")

    if args.blackListFileName:
        cmd.extend(
            [
                "--blackListFileName",
                str(Path(args.blackListFileName).resolve())
            ]
        )

    return cmd


def main():

    args = parse_args()

    if not args.skip_check:
        check_deeptools()

    validate_inputs([args.bigwig1, args.bigwig2])

    Path(args.output).parent.mkdir(
        parents=True,
        exist_ok=True
    )

    cmd = build_command(args)

    print(
        f"[INFO] Running bigwigCompare "
        f"({args.operation})"
    )

    if args.verbose:
        print("[CMD] " + " ".join(cmd))

    result = subprocess.run(
        cmd,
        capture_output=False,
        text=True
    )

    if result.returncode != 0:
        sys.exit(
            f"[ERROR] bigwigCompare exited "
            f"with code {result.returncode}"
        )

    print(
        f"[DONE] Output written to: "
        f"{args.output}"
    )


if __name__ == "__main__":
    main()