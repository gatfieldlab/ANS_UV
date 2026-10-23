#!/usr/bin/env python3
"""
compute_matrix.py

A comprehensive Python wrapper exposing (essentially) every option available
in deepTools' `computeMatrix` command line tool, for both the
`scale-regions` and `reference-point` sub-commands.

It builds the appropriate `computeMatrix` command line from the supplied
arguments and runs it via subprocess, giving you a scriptable, reusable,
and version-controllable way to launch computeMatrix jobs instead of
hand-writing (or copy-pasting) long CLI invocations.

Requirements
------------
deepTools must be installed and `computeMatrix` must be on PATH:
    pip install deeptools
    # or
    conda install -c bioconda deeptools

Basic usage
-----------
The tool has two modes, matching computeMatrix itself:

1) reference-point: values are computed relative to a single point in each
   region (TSS, TES, or center) — typical for signal enrichment around a
   fixed feature.

    python compute_matrix.py reference-point \\
        -S sample1.bw sample2.bw \\
        -R regions1.bed regions2.bed \\
        -o matrix.gz \\
        --referencePoint TSS \\
        --beforeRegionStartLength 1000 \\
        --afterRegionStartLength 1000 \\
        --binSize 10

2) scale-regions: every region is stretched or shrunk to a common length —
   typical for comparing signal across gene bodies of varying size.

    python compute_matrix.py scale-regions \\
        -S sample1.bw -R genes.bed -o matrix.gz \\
        --regionBodyLength 1000 \\
        --beforeRegionStartLength 500 \\
        --afterRegionStartLength 500

-S accepts one or more bigWig files and -R one or more BED/GTF files;
list several values to compute a single matrix across multiple
samples/region sets at once, exactly as in the underlying computeMatrix
tool.

Run `python compute_matrix.py reference-point --help` or
`... scale-regions --help` to see every available option and its default.
"""

import argparse
import shlex
import subprocess
import sys


def add_common_arguments(parser):
    """Options shared by both scale-regions and reference-point modes."""

    required = parser.add_argument_group("Required arguments")
    required.add_argument(
        "-S", "--scoreFileName", nargs="+", required=True,
        help="One or more bigWig file(s) containing the scores to be plotted.",
    )
    required.add_argument(
        "-R", "--regionsFileName", nargs="+", required=True,
        help="One or more BED/GTF file(s) containing regions of interest.",
    )
    required.add_argument(
        "-o", "--outFileName", required=True,
        help="Output file name for the gzipped matrix (consumed by plotHeatmap/plotProfile).",
    )

    output = parser.add_argument_group("Output options")
    output.add_argument(
        "--outFileNameMatrix",
        help="Optional: tab-separated file with the raw values underlying the matrix.",
    )
    output.add_argument(
        "--outFileSortedRegions",
        help="Optional: BED file with regions in the order/filtering used for the matrix.",
    )

    opt = parser.add_argument_group("General options")
    opt.add_argument(
        "-bs", "--binSize", type=int, default=10,
        help="Length, in bases, of the non-overlapping bins used to average scores. Default: 10.",
    )
    opt.add_argument(
        "--sortRegions", choices=["descend", "ascend", "no", "keep"], default="keep",
        help="Whether/how to sort regions. 'keep' preserves the input order. Default: keep.",
    )
    opt.add_argument(
        "--sortUsing",
        choices=["mean", "median", "max", "min", "sum", "region_length"],
        default="mean",
        help="Statistic used for sorting when --sortRegions is not 'no'/'keep'. Default: mean.",
    )
    opt.add_argument(
        "--sortUsingSamples", type=int, nargs="+",
        help="1-indexed list of sample numbers to base sorting on, e.g. 1 3.",
    )
    opt.add_argument(
        "--averageTypeBins",
        choices=["mean", "median", "min", "max", "std", "sum"],
        default="mean",
        help="Method used to summarize multiple values falling in one bin. Default: mean.",
    )
    opt.add_argument(
        "--missingDataAsZero", action="store_true",
        help="Treat missing (NA) data points as zero instead of excluding them.",
    )
    opt.add_argument(
        "--skipZeros", action="store_true",
        help="Skip regions that have all zero/NaN values across all samples.",
    )
    opt.add_argument(
        "--minThreshold", type=float,
        help="Exclude regions where any bin value is below this threshold.",
    )
    opt.add_argument(
        "--maxThreshold", type=float,
        help="Exclude regions where any bin value is above this threshold.",
    )
    opt.add_argument(
        "-bl", "--blackListFileName",
        help="BED file of blacklisted regions to exclude from the analysis.",
    )
    opt.add_argument(
        "--samplesLabel", nargs="+",
        help="Labels for each sample, in the same order as -S. Defaults to file names.",
    )
    opt.add_argument(
        "--smartLabels", action="store_true",
        help="Infer sample/region labels from file names (strip path and extension).",
    )
    opt.add_argument(
        "-p", "--numberOfProcessors", default="max",
        help="Number of processors to use: an int, 'max', or 'max/2'. Default: max.",
    )
    opt.add_argument("-q", "--quiet", action="store_true", help="Suppress progress messages.")
    opt.add_argument("--verbose", action="store_true", help="Verbose output.")

    gtf = parser.add_argument_group("GTF / gene-model options (relevant only if -R is a GTF)")
    gtf.add_argument(
        "--metagene", action="store_true",
        help="Compute the matrix over joined exons only (metagene model), for GTF input.",
    )
    gtf.add_argument(
        "--transcriptID", default="transcript",
        help="GTF feature type to treat as a transcript. Default: transcript.",
    )
    gtf.add_argument(
        "--exonID", default="exon",
        help="GTF feature type to treat as an exon. Default: exon.",
    )
    gtf.add_argument(
        "--transcript_id_designator", default="transcript_id",
        help="GTF attribute key linking exons to their transcript. Default: transcript_id.",
    )


def build_reference_point_parser(subparsers):
    p = subparsers.add_parser(
        "reference-point",
        help="Compute the matrix relative to a single reference point (TSS, TES, or center).",
    )
    add_common_arguments(p)
    rp = p.add_argument_group("reference-point specific options")
    rp.add_argument(
        "--referencePoint", choices=["TSS", "TES", "center"], default="TSS",
        help="The reference point within each region. Default: TSS.",
    )
    rp.add_argument(
        "-b", "--beforeRegionStartLength", type=int, default=0,
        help="Distance (bp) upstream of the reference point to include. Default: 0.",
    )
    rp.add_argument(
        "-a", "--afterRegionStartLength", type=int, default=0,
        help="Distance (bp) downstream of the reference point to include. Default: 0.",
    )
    rp.add_argument(
        "--nanAfterEnd", action="store_true",
        help="Set values past the end of a region to NaN (useful when -a extends beyond it).",
    )
    return p


def build_scale_regions_parser(subparsers):
    p = subparsers.add_parser(
        "scale-regions",
        help="Compute the matrix by scaling all regions to a common length.",
    )
    add_common_arguments(p)
    sr = p.add_argument_group("scale-regions specific options")
    sr.add_argument(
        "-m", "--regionBodyLength", type=int, default=1000,
        help="Length (bp) all regions are scaled to. Default: 1000.",
    )
    sr.add_argument(
        "-b", "--beforeRegionStartLength", type=int, default=0,
        help="Distance (bp) upstream of the region start to include. Default: 0.",
    )
    sr.add_argument(
        "-a", "--afterRegionStartLength", type=int, default=0,
        help="Distance (bp) downstream of the region end to include. Default: 0.",
    )
    sr.add_argument(
        "--unscaled5prime", type=int, default=0,
        help="Number of bp at the 5' end of each region left unscaled. Default: 0.",
    )
    sr.add_argument(
        "--unscaled3prime", type=int, default=0,
        help="Number of bp at the 3' end of each region left unscaled. Default: 0.",
    )
    sr.add_argument(
        "--startLabel", default="TSS",
        help="Label used for the region-start tick in downstream plots. Default: TSS.",
    )
    sr.add_argument(
        "--endLabel", default="TES",
        help="Label used for the region-end tick in downstream plots. Default: TES.",
    )
    return p


def build_command(mode, args):
    """Translate a parsed argparse Namespace into a computeMatrix CLI command list."""
    cmd = ["computeMatrix", mode]

    cmd += ["--scoreFileName"] + args.scoreFileName
    cmd += ["--regionsFileName"] + args.regionsFileName
    cmd += ["--outFileName", args.outFileName]

    if args.outFileNameMatrix:
        cmd += ["--outFileNameMatrix", args.outFileNameMatrix]
    if args.outFileSortedRegions:
        cmd += ["--outFileSortedRegions", args.outFileSortedRegions]

    cmd += ["--binSize", str(args.binSize)]
    cmd += ["--sortRegions", args.sortRegions]
    cmd += ["--sortUsing", args.sortUsing]
    if args.sortUsingSamples:
        cmd += ["--sortUsingSamples"] + [str(x) for x in args.sortUsingSamples]
    cmd += ["--averageTypeBins", args.averageTypeBins]
    if args.missingDataAsZero:
        cmd.append("--missingDataAsZero")
    if args.skipZeros:
        cmd.append("--skipZeros")
    if args.minThreshold is not None:
        cmd += ["--minThreshold", str(args.minThreshold)]
    if args.maxThreshold is not None:
        cmd += ["--maxThreshold", str(args.maxThreshold)]
    if args.blackListFileName:
        cmd += ["--blackListFileName", args.blackListFileName]
    if args.samplesLabel:
        cmd += ["--samplesLabel"] + args.samplesLabel
    if args.smartLabels:
        cmd.append("--smartLabels")
    cmd += ["--numberOfProcessors", str(args.numberOfProcessors)]
    if args.quiet:
        cmd.append("--quiet")
    if args.verbose:
        cmd.append("--verbose")

    if args.metagene:
        cmd.append("--metagene")
    cmd += ["--transcriptID", args.transcriptID]
    cmd += ["--exonID", args.exonID]
    cmd += ["--transcript_id_designator", args.transcript_id_designator]

    cmd += ["-b", str(args.beforeRegionStartLength)]
    cmd += ["-a", str(args.afterRegionStartLength)]

    if mode == "reference-point":
        cmd += ["--referencePoint", args.referencePoint]
        if args.nanAfterEnd:
            cmd.append("--nanAfterEnd")
    else:  # scale-regions
        cmd += ["--regionBodyLength", str(args.regionBodyLength)]
        cmd += ["--unscaled5prime", str(args.unscaled5prime)]
        cmd += ["--unscaled3prime", str(args.unscaled3prime)]
        cmd += ["--startLabel", args.startLabel]
        cmd += ["--endLabel", args.endLabel]

    return cmd


def main():
    parser = argparse.ArgumentParser(
        description="Python wrapper exposing all deepTools computeMatrix options."
    )
    subparsers = parser.add_subparsers(dest="mode", required=True, help="computeMatrix sub-command")
    build_reference_point_parser(subparsers)
    build_scale_regions_parser(subparsers)

    args = parser.parse_args()
    cmd = build_command(args.mode, args)

    print("Running command:\n  " + " ".join(shlex.quote(c) for c in cmd), file=sys.stderr)
    result = subprocess.run(cmd)
    sys.exit(result.returncode)


if __name__ == "__main__":
    main()