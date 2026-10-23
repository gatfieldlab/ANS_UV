#!/usr/bin/env python3
"""
plot_heatmap.py

A comprehensive Python wrapper exposing (essentially) every option available
in deepTools' `plotHeatmap` command line tool.

It builds the appropriate `plotHeatmap` command line from the supplied
arguments and runs it via subprocess, giving you a scriptable, reusable,
and version-controllable way to launch plotHeatmap jobs instead of
hand-writing (or copy-pasting) long CLI invocations.

Requirements
------------
deepTools must be installed and `plotHeatmap` must be on PATH:
    pip install deeptools
    # or
    conda install -c bioconda deeptools

plotHeatmap draws a heatmap (with an optional summary profile on top) for
scores over a set of genomic regions, based on the matrix produced by
`computeMatrix`. It is typically used right after computeMatrix.

Basic usage
-----------
    python plot_heatmap.py \\
        -m matrix.gz \\
        -o heatmap.png \\
        --colorMap RdBu \\
        --whatToShow "plot, heatmap and colorbar" \\
        --plotTitle "Signal around TSS"

Run `python plot_heatmap.py --help` to see every available option
and its default.
"""

import argparse
import shlex
import subprocess
import sys


def build_parser():
    parser = argparse.ArgumentParser(
        description="Python wrapper exposing all deepTools plotHeatmap options."
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
        "--outFileNameMatrix",
        help="Optional: save the underlying (possibly sorted/clustered) matrix values to this file.",
    )
    output.add_argument(
        "--outFileNameData",
        help="Optional: save the values underlying the summary/profile plot to this tab-separated file.",
    )
    output.add_argument(
        "--outFileSortedRegions",
        help="Optional: BED file with regions in the order/clusters used for the heatmap.",
    )

    proc = parser.add_argument_group("Region clustering / sorting")
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
    proc.add_argument(
        "--sortRegions", choices=["descend", "ascend", "no", "keep"], default="descend",
        help="Whether/how to sort regions before plotting. Default: descend.",
    )
    proc.add_argument(
        "--sortUsing",
        choices=["mean", "median", "max", "min", "sum", "region_length"],
        default="mean",
        help="Statistic used for sorting when --sortRegions is not 'no'/'keep'. Default: mean.",
    )
    proc.add_argument(
        "--sortUsingSamples", type=int, nargs="+",
        help="1-indexed list of sample numbers to base sorting on, e.g. 1 3.",
    )

    style = parser.add_argument_group("Appearance: colors")
    style.add_argument(
        "--colorMap", nargs="+", default=["RdYlBu"],
        help="Matplotlib colormap name(s), one per sample, or one for all. Default: RdYlBu.",
    )
    style.add_argument(
        "--colorList", nargs="+",
        help="Custom list of colors to build a colormap from (overrides --colorMap), "
             "e.g. white,blue or '#FFFFFF,#0000FF'.",
    )
    style.add_argument(
        "--colorNumber", type=int, default=256,
        help="Number of levels used to build a custom colormap from --colorList. Default: 256.",
    )
    style.add_argument(
        "--missingDataColor", default="black",
        help="Color used for missing data points. Default: black.",
    )
    style.add_argument(
        "--alpha", type=float, default=1.0,
        help="Transparency of the heatmap (0 to 1). Default: 1.0.",
    )
    style.add_argument(
        "--zMin", type=float, nargs="+",
        help="Minimum value(s) for the heatmap intensity scale, one per sample or one for all.",
    )
    style.add_argument(
        "--zMax", type=float, nargs="+",
        help="Maximum value(s) for the heatmap intensity scale, one per sample or one for all.",
    )
    style.add_argument(
        "--interpolationMethod", default="auto",
        help="Interpolation method for the heatmap image (e.g. auto, nearest, bilinear, bicubic). "
             "Default: auto.",
    )

    layout = parser.add_argument_group("Appearance: layout")
    layout.add_argument(
        "--dpi", type=int, default=200,
        help="Image resolution in dots-per-inch. Default: 200.",
    )
    layout.add_argument(
        "--heatmapHeight", type=float, default=28,
        help="Height in cm of the heatmap (range 3-100). Default: 28.",
    )
    layout.add_argument(
        "--heatmapWidth", type=float, default=4,
        help="Width in cm of the heatmap per sample (range 1-100). Default: 4.",
    )
    layout.add_argument(
        "--whatToShow",
        choices=["plot, heatmap and colorbar", "plot and heatmap", "heatmap only", "heatmap and colorbar"],
        default="plot, heatmap and colorbar",
        help="Which elements to draw. Default: 'plot, heatmap and colorbar'.",
    )
    layout.add_argument(
        "--boxAroundHeatmaps", choices=["yes", "no"], default="yes",
        help="Whether to draw a box around each heatmap. Default: yes.",
    )
    layout.add_argument(
        "--perGroup", action="store_true",
        help="Plot all samples for each region group together, rather than all groups per sample.",
    )
    layout.add_argument(
        "--averageTypeSummaryPlot",
        choices=["mean", "median", "min", "max", "std", "sum"],
        default="mean",
        help="Statistic shown in the summary plot above the heatmap. Default: mean.",
    )
    layout.add_argument(
        "--numPlotsPerRow", type=int,
        help="Number of heatmap panels to draw per row before wrapping to a new row.",
    )

    labels = parser.add_argument_group("Labels and titles")
    labels.add_argument("--plotTitle", default="", help="Title printed on top of the plot.")
    labels.add_argument("--xAxisLabel", default="", help="Label for the x-axis of the summary plot.")
    labels.add_argument("--yAxisLabel", default="", help="Label for the y-axis of the summary plot.")
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
        "--labelRotation", type=float, default=0.0,
        help="Rotation, in degrees, of the x-axis labels. Default: 0.",
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
        "--legendLocation",
        default="best",
        help="Location of the legend in the summary plot (e.g. best, upper-right, upper-left, "
             "lower-right, lower-left, none, ...). Default: best.",
    )

    misc = parser.add_argument_group("Miscellaneous")
    misc.add_argument(
        "--yMin", type=float, nargs="+", help="Minimum value(s) for the y-axis of the summary plot.",
    )
    misc.add_argument(
        "--yMax", type=float, nargs="+", help="Maximum value(s) for the y-axis of the summary plot.",
    )
    misc.add_argument(
        "--plotFileFormat", choices=["png", "pdf", "svg", "eps", "plotly"],
        help="Force a specific image format regardless of the -out file extension.",
    )
    misc.add_argument("--verbose", action="store_true", help="Verbose output.")

    return parser


def build_command(args):
    cmd = ["plotHeatmap", "--matrixFile", args.matrixFile, "--outFileName", args.outFileName]

    if args.outFileNameMatrix:
        cmd += ["--outFileNameMatrix", args.outFileNameMatrix]
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
    cmd += ["--sortRegions", args.sortRegions]
    cmd += ["--sortUsing", args.sortUsing]
    if args.sortUsingSamples:
        cmd += ["--sortUsingSamples"] + [str(x) for x in args.sortUsingSamples]

    cmd += ["--colorMap"] + args.colorMap
    if args.colorList:
        cmd += ["--colorList"] + args.colorList
    cmd += ["--colorNumber", str(args.colorNumber)]
    cmd += ["--missingDataColor", args.missingDataColor]
    cmd += ["--alpha", str(args.alpha)]
    if args.zMin:
        cmd += ["--zMin"] + [str(x) for x in args.zMin]
    if args.zMax:
        cmd += ["--zMax"] + [str(x) for x in args.zMax]
    cmd += ["--interpolationMethod", args.interpolationMethod]

    cmd += ["--dpi", str(args.dpi)]
    cmd += ["--heatmapHeight", str(args.heatmapHeight)]
    cmd += ["--heatmapWidth", str(args.heatmapWidth)]
    cmd += ["--whatToShow", args.whatToShow]
    cmd += ["--boxAroundHeatmaps", args.boxAroundHeatmaps]
    if args.perGroup:
        cmd.append("--perGroup")
    cmd += ["--averageTypeSummaryPlot", args.averageTypeSummaryPlot]
    if args.numPlotsPerRow:
        cmd += ["--numPlotsPerRow", str(args.numPlotsPerRow)]

    if args.plotTitle:
        cmd += ["--plotTitle", args.plotTitle]
    if args.xAxisLabel:
        cmd += ["--xAxisLabel", args.xAxisLabel]
    if args.yAxisLabel:
        cmd += ["--yAxisLabel", args.yAxisLabel]
    cmd += ["--startLabel", args.startLabel]
    cmd += ["--endLabel", args.endLabel]
    cmd += ["--refPointLabel", args.refPointLabel]
    cmd += ["--labelRotation", str(args.labelRotation)]
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