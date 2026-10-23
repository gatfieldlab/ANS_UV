#!/usr/bin/env python3
"""
bigwig_mean.py
--------------
Compute the mean of multiple BigWig files (already in CPM or any
comparable unit) using deeptools bigwigAverage.

Usage examples
--------------
# Average all .bw files in a directory:
python bigwig_mean.py -i /path/to/bw/*.bw -o mean_output.bw

# Average a specific list of files:
python bigwig_mean.py -i rep1.bw rep2.bw rep3.bw -o mean_output.bw

# Specify bin size and number of cores:
python bigwig_mean.py -i *.bw -o mean_output.bw --binSize 10 --numberOfProcessors 8

# Skip the deeptools check (useful in environments where it is known to be installed):
python bigwig_mean.py -i *.bw -o mean_output.bw --skip-check
"""

import argparse
import subprocess
import sys
import shutil
from pathlib import Path


def check_deeptools() -> None:
    """Verify that deeptools (bigwigAverage) is available on PATH."""
    if shutil.which("bigwigAverage") is None:
        sys.exit(
            "[ERROR] 'bigwigAverage' not found on PATH.\n"
            "Install deeptools with:  pip install deeptools\n"
            "or via conda:            conda install -c bioconda deeptools"
        )


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Compute the mean of BigWig files using deeptools bigwigAverage.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument(
        "-i", "--input",
        nargs="+",
        required=True,
        metavar="FILE",
        help="Input BigWig files (space-separated or glob-expanded by the shell).",
    )
    parser.add_argument(
        "-o", "--output",
        required=True,
        metavar="FILE",
        help="Output BigWig file for the mean track.",
    )
    parser.add_argument(
        "--binSize",
        type=int,
        default=10,
        metavar="INT",
        help="Bin size in bases used for averaging.",
    )
    parser.add_argument(
        "--numberOfProcessors", "-p",
        type=int,
        default=4,
        metavar="INT",
        help="Number of parallel processors.",
    )
    parser.add_argument(
        "--outFileFormat",
        choices=["bigwig", "bedgraph"],
        default="bigwig",
        help="Output file format.",
    )
    parser.add_argument(
        "--skipNonCoveredRegions",
        action="store_true",
        help="If set, bins without coverage in ALL files are skipped (set to NaN).",
    )
    parser.add_argument(
        "--blackListFileName",
        metavar="BED",
        default=None,
        help="Optional BED file of regions to exclude (e.g. ENCODE blacklist).",
    )
    parser.add_argument(
        "--skip-check",
        action="store_true",
        help="Skip the deeptools installation check.",
    )
    parser.add_argument(
        "--verbose", "-v",
        action="store_true",
        help="Print the exact command before running it.",
    )
    return parser.parse_args()


def build_command(args: argparse.Namespace) -> list[str]:
    """Build the bigwigAverage command list."""
    input_files = [str(Path(f).resolve()) for f in args.input]
    output_file = str(Path(args.output).resolve())

    cmd = [
        "bigwigAverage",
        "--bigwigs", *input_files,
        "--outFileName", output_file,
        "--binSize", str(args.binSize),
        "--numberOfProcessors", str(args.numberOfProcessors),
        "--outFileFormat", args.outFileFormat,
    ]

    if args.skipNonCoveredRegions:
        cmd.append("--skipNonCoveredRegions")

    if args.blackListFileName:
        cmd.extend(["--blackListFileName", str(Path(args.blackListFileName).resolve())])

    return cmd


def validate_inputs(input_files: list[str]) -> None:
    """Check that all input files exist and are readable."""
    missing = [f for f in input_files if not Path(f).is_file()]
    if missing:
        sys.exit(
            "[ERROR] The following input files were not found:\n"
            + "\n".join(f"  {f}" for f in missing)
        )


def main() -> None:
    args = parse_args()

    if not args.skip_check:
        check_deeptools()

    validate_inputs(args.input)

    # Ensure output directory exists
    out_path = Path(args.output)
    out_path.parent.mkdir(parents=True, exist_ok=True)

    cmd = build_command(args)

    print(f"[INFO] Averaging {len(args.input)} BigWig file(s) → {args.output}")
    if args.verbose:
        print("[CMD] " + " ".join(cmd))

    result = subprocess.run(cmd, capture_output=False, text=True)

    if result.returncode != 0:
        sys.exit(f"[ERROR] bigwigAverage exited with code {result.returncode}.")

    print(f"[DONE] Mean BigWig written to: {args.output}")


if __name__ == "__main__":
    main()