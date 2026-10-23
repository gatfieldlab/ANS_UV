#!/usr/bin/env python3
"""
plot_profile.py

A comprehensive Python wrapper exposing (essentially) every option available
in deepTools' `plotProfile` command line tool.

It builds the appropriate `plotProfile` command line from the supplied
arguments and runs it via subprocess, giving you a scriptable, reusable,
and version-controllable way to launch plotProfile jobs instead of
hand-writing (or copy-pasting) long CLI invocations.

Requirements
------------
deepTools must be installed and `plotProfile` must be on PATH:
    pip install deeptools
    # or
    conda install -c bioconda deeptools

plotProfile draws a line/area profile (average signal over the region set)
based on the matrix produced by `computeMatrix`. It is typically used right
after computeMatrix, as an alternative or complement to plotHeatmap.

Basic usage
-----------
    python plot_profile.py \\
        -m matrix.gz \\
        -o profile.png \\
        --plotType lines \\
        --perGroup \\
        --plotTitle "Signal around TSS"

Run `python plot_profile.py --help` to see every available option
and its default.
"""

import argparse
import shlex
import subprocess
import sys


def build_parser():
    parser = argparse.ArgumentParser(
        description="Python wrapper exposing all deepTools plotProfile options."
    )

    required = parser.add_argument_group("Required arguments")
    required.add_argument(
        "-m", "--matrixFile", required=True,
        help="Matrix file produced by computeMatrix (gzipped .gz file).",
    )
    required.add_argument(
        "-out", "--outFileName", required=True,
        help="Output image file (extension determines format: png, pdf, svg, eps).",
    )

    output = parser.add_argument_group("Output options")
    output.add_argument(
        "--outFileNameData",
        help="Optional: save the values underlying the profile plot to this tab-separated file.",
    )
    output.add_argument(
        "--outFileSortedRegions",
        help="Optional: BED file with regions in the order/clusters used for the plot.",
    )

    proc = parser.add_argument_group("Region clustering")
    proc.add_argument(
        "--kmeans", type=int,
        help="Number of k-means clusters to compute on the regions (mutually exclusive with --hclust).",
    )
    proc.add_argument(
        "--hclust", type=int,
        help="Number of hierarchical clusters to compute on the regions (mutually exclusive with --kmeans).",
    )
    proc.add_argument(
        "--silhouette", action="store_true",
        help="Compute silhouette score when clustering (--kmeans/--hclust); slower for many regions.",
    )
    proc.add_argument(
        "--linkage", default="average",
        help="Linkage method used by --hclust (any method accepted by scipy, e.g. average, complete, "
             "single, ward). Default: average.",
    )

    calc = parser.add_argument_group("Profile calculation")
    calc.add_argument(
        "--averageType",
        choices=["mean", "median", "min", "max", "std", "sum"],
        default="mean",
        help="Statistic used to summarize values in each bin across regions. Default: mean.",
    )
    calc.add_argument(
        "--nanToZero", action="store_true",
        help="Treat missing (NaN) values as zero instead of excluding them.",
    )
    calc.add_argument(
        "--proj", action="store_true",
        help="For plotType 'se'/'std', plot the individual sample lines under the shaded area.",
    )

    style = parser.add_argument_group("Appearance: plot style")
    style.add_argument(
        "--plotType",
        choices=["lines", "fill", "se", "std", "overlapped_lines", "heatmap"],
        default="lines",
        help="How to render the profile: line plot, filled area, standard error/deviation band, "
             "overlapped individual lines, or a compact heatmap-style profile. Default: lines.",
    )
    style.add_argument(
        "--colors", nargs="+",
        help="List of colors, one per sample or per group (e.g. red blue '#00FF00').",
    )
    style.add_argument(
        "--plotHeight", type=float, default=7,
        help="Height in cm of each plot panel (range 3-100). Default: 7.",
    )
    style.add_argument(
        "--plotWidth", type=float, default=11,
        help="Width in cm of each plot panel (range 1-100). Default: 11.",
    )
    style.add_argument(
        "--dpi", type=int, default=200,
        help="Image resolution in dots-per-inch. Default: 200.",
    )
    style.add_argument(
        "--numPlotsPerRow", type=int, default=8,
        help="Number of panels to draw per row before wrapping to a new row. Default: 8.",
    )
    style.add_argument(
        "--perGroup", action="store_true",
        help="Plot all samples for each region group together in one panel, rather than all "
             "groups per sample.",
    )

    labels = parser.add_argument_group("Labels and titles")
    labels.add_argument("--plotTitle", default="", help="Title printed on top of the plot.")
    labels.add_argument("--xAxisLabel", default="", help="Label for the x-axis.")
    labels.add_argument("--yAxisLabel", default="", help="Label for the y-axis.")
    labels.add_argument(
        "--startLabel", default="TSS",
        help="Label for the region start tick (scale-regions matrices). Default: TSS.",
    )
    labels.add_argument(
        "--endLabel", default="TES",
        help="Label for the region end tick (scale-regions matrices). Default: TES.",
    )
    labels.add_argument(
        "--refPointLabel", default="TSS",
        help="Label for the reference point tick (reference-point matrices). Default: TSS.",
    )
    labels.add_argument(
        "--regionsLabel", nargs="+",
        help="Labels for each region set/group, in the order used in computeMatrix -R.",
    )
    labels.add_argument(
        "--samplesLabel", nargs="+",
        help="Labels for each sample/bigWig, in the order used in computeMatrix -S.",
    )
    labels.add_argument(
        "--legendLocation", default="best",
        help="Location of the legend (e.g. best, upper-right, upper-left, lower-right, "
             "lower-left, none, ...). Default: best.",
    )

    misc = parser.add_argument_group("Miscellaneous")
    misc.add_argument(
        "--yMin", type=float, nargs="+", help="Minimum value(s) for the y-axis, one per panel or one for all.",
    )
    misc.add_argument(
        "--yMax", type=float, nargs="+", help="Maximum value(s) for the y-axis, one per panel or one for all.",
    )
    misc.add_argument(
        "--plotFileFormat", choices=["png", "pdf", "svg", "eps", "plotly"],
        help="Force a specific image format regardless of the -out file extension.",
    )
    misc.add_argument("--verbose", action="store_true", help="Verbose output.")

    return parser


def build_command(args):
    cmd = ["plotProfile", "--matrixFile", args.matrixFile, "--outFileName", args.outFileName]

    if args.outFileNameData:
        cmd += ["--outFileNameData", args.outFileNameData]
    if args.outFileSortedRegions:
        cmd += ["--outFileSortedRegions", args.outFileSortedRegions]

    if args.kmeans:
        cmd += ["--kmeans", str(args.kmeans)]
    if args.hclust:
        cmd += ["--hclust", str(args.hclust)]
    if args.silhouette:
        cmd.append("--silhouette")
    cmd += ["--linkage", args.linkage]

    cmd += ["--averageType", args.averageType]
    if args.nanToZero:
        cmd.append("--nanToZero")
    if args.proj:
        cmd.append("--proj")

    cmd += ["--plotType", args.plotType]
    if args.colors:
        cmd += ["--colors"] + args.colors
    cmd += ["--plotHeight", str(args.plotHeight)]
    cmd += ["--plotWidth", str(args.plotWidth)]
    cmd += ["--dpi", str(args.dpi)]
    cmd += ["--numPlotsPerRow", str(args.numPlotsPerRow)]
    if args.perGroup:
        cmd.append("--perGroup")

    if args.plotTitle:
        cmd += ["--plotTitle", args.plotTitle]
    if args.xAxisLabel:
        cmd += ["--xAxisLabel", args.xAxisLabel]
    if args.yAxisLabel:
        cmd += ["--yAxisLabel", args.yAxisLabel]
    cmd += ["--startLabel", args.startLabel]
    cmd += ["--endLabel", args.endLabel]
    cmd += ["--refPointLabel", args.refPointLabel]
    if args.regionsLabel:
        cmd += ["--regionsLabel"] + args.regionsLabel
    if args.samplesLabel:
        cmd += ["--samplesLabel"] + args.samplesLabel
    cmd += ["--legendLocation", args.legendLocation]

    if args.yMin:
        cmd += ["--yMin"] + [str(x) for x in args.yMin]
    if args.yMax:
        cmd += ["--yMax"] + [str(x) for x in args.yMax]
    if args.plotFileFormat:
        cmd += ["--plotFileFormat", args.plotFileFormat]
    if args.verbose:
        cmd.append("--verbose")

    return cmd


def main():
    parser = build_parser()
    args = parser.parse_args()

    if args.kmeans and args.hclust:
        parser.error("--kmeans and --hclust are mutually exclusive; choose one clustering method.")

    cmd = build_command(args)
    print("Running command:\n  " + " ".join(shlex.quote(c) for c in cmd), file=sys.stderr)
    result = subprocess.run(cmd)
    sys.exit(result.returncode)


if __name__ == "__main__":
    main()