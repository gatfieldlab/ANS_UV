# SGA → bigWig → deepTools Pipeline

A set of Python scripts to go from `.sga` count files to **heatmap** and
**average profile** figures (deepTools), using **transcripts** (rather than a
genome) as references.

```
.sga ──► 01_sga_to_bigwig ──► 02_bigwig_mean ──► 03_bigwig_compare ──► 04_compute_matrix ──┬─► 05_plotHeatmap
         (CPM + 1 bp bigWig)   (replicate mean)    (log2, difference…)   (matrix)           └─► 05_bis_plotProfile
```

| Step | Script | Underlying tool | Purpose |
|-----:|--------|-----------------|---------|
| 1 | `01_sga_to_bigwig.py` | `bedGraphToBigWig` (UCSC) + pandas | CPM normalization and `.sga` → bigWig conversion (1 bp) |
| 2 | `02_bigwig_mean.py` | `bigwigAverage` (deepTools) | Mean of several bigWig files (replicates) |
| 3 | `03_bigwig_compare.py` | `bigwigCompare` (deepTools) | Compare two bigWig files (log2, ratio, subtraction…) |
| 4 | `04_compute_matrix.py` | `computeMatrix` (deepTools) | Score matrix around/over regions |
| 5 | `05_plotHeatmap.py` | `plotHeatmap` (deepTools) | Heatmap (with summary profile and colorbar) |
| 5bis | `05_bis_plotProfile.py` | `plotProfile` (deepTools) | Average profile (lines, fill, SE/STD bands…) |

Scripts 02 to 05bis are thin *wrappers*: they build the corresponding deepTools
command line and run it via `subprocess`. They print the exact command that is
run (on stderr, or with `--verbose` for 02 and 03).

---

## Requirements

- Python ≥ 3.9 (`list[str]` annotations are used)
- [deepTools](https://deeptools.readthedocs.io/) (scripts 02 to 05bis)
- `bedGraphToBigWig` from the UCSC utilities (script 01)
- `pandas` (script 01)

```bash
conda install -c bioconda deeptools ucsc-bedgraphtobigwig pandas
# or
pip install deeptools pandas   # bedGraphToBigWig must be installed separately
```

Check your installation:

```bash
bigwigAverage --version
bigwigCompare --version
computeMatrix --version
which bedGraphToBigWig
```

---

## 1. `01_sga_to_bigwig.py` — SGA → CPM → bigWig

For each `.sga` file in a folder, the script:

1. normalizes counts to **CPM** (`count × 1e6 / total counts`);
2. writes a temporary bedGraph at 1 bp resolution;
3. converts it to bigWig with `bedGraphToBigWig`.

**Expected SGA format** (TSV, no header):
`transcript  ref_type  position  strand  count`

Chromosome sizes (`chrom.sizes`) are computed from a transcript FASTA.

```bash
python 01_sga_to_bigwig.py /path/to/sga_folder --fasta transcripts.fa
python 01_sga_to_bigwig.py /path/to/sga_folder --fasta transcripts.fa --keep-cpm
```

| Option | Description |
|--------|-------------|
| `input_dir` | Folder containing the `.sga` files |
| `--fasta` | Transcript FASTA (required) |
| `--keep-cpm` | Also save the CPM-normalized `.sga` files |

**Outputs** (created next to the input folder):

- `bigwig_CPM/<name>_CPM.bw` — one bigWig per `.sga` file
- `bigwig_CPM/chrom.sizes`
- `<folder>_normalized_CPM/<name>_normCPM.sga` — only with `--keep-cpm`

Files whose total count is 0 are skipped.

---

## 2. `02_bigwig_mean.py` — Mean of bigWig files

Computes the mean of several bigWig files (typically replicates already in CPM)
using `bigwigAverage`.

```bash
python 02_bigwig_mean.py -i rep1.bw rep2.bw rep3.bw -o mean.bw
python 02_bigwig_mean.py -i bigwig_CPM/*.bw -o mean.bw --binSize 1 -p 8
```

| Option | Default | Description |
|--------|---------|-------------|
| `-i, --input` | — | Input bigWig files (required) |
| `-o, --output` | — | Output bigWig (required) |
| `--binSize` | `10` | Bin size in bases |
| `-p, --numberOfProcessors` | `4` | Number of processors |
| `--outFileFormat` | `bigwig` | `bigwig` or `bedgraph` |
| `--skipNonCoveredRegions` | off | Skip bins without coverage |
| `--blackListFileName` | — | BED file of regions to exclude |
| `--skip-check` | off | Skip the deepTools availability check |
| `-v, --verbose` | off | Print the command being run |

---

## 3. `03_bigwig_compare.py` — Compare two bigWig files

Wrapper around `bigwigCompare`.

```bash
# log2(treatment / control)
python 03_bigwig_compare.py -b1 treatment.bw -b2 control.bw -o log2.bw --operation log2

# Difference
python 03_bigwig_compare.py -b1 s1.bw -b2 s2.bw -o diff.bw --operation subtract
```

| Option | Default | Description |
|--------|---------|-------------|
| `-b1, --bigwig1` | — | First bigWig (required) |
| `-b2, --bigwig2` | — | Second bigWig (required) |
| `-o, --output` | — | Output file (required) |
| `--operation` | `log2` | `log2`, `ratio`, `subtract`, `add`, `mean`, `reciprocal_ratio`, `first`, `second` |
| `--binSize` | `10` | Bin size in bp |
| `-p, --numberOfProcessors` | `4` | Number of processors |
| `--pseudocount PC1 PC2` | `1 1` | Pseudocounts (only passed for `log2`, `ratio`, `reciprocal_ratio`) |
| `--skipNonCoveredRegions` | off | Skip regions without coverage |
| `--blackListFileName` | — | Blacklist BED file |
| `--skip-check` | off | Skip the deepTools availability check |
| `-v, --verbose` | off | Print the command being run |

---

## 4. `04_compute_matrix.py` — computeMatrix

Wrapper exposing nearly all `computeMatrix` options, with both sub-commands:

- **`reference-point`**: values relative to a single point (TSS, TES, center);
- **`scale-regions`**: all regions are rescaled to a common length.

```bash
# Around a reference point
python 04_compute_matrix.py reference-point \
    -S sample1.bw sample2.bw \
    -R regions.bed \
    -o matrix.gz \
    --referencePoint TSS -b 1000 -a 1000 --binSize 10

# Region scaling
python 04_compute_matrix.py scale-regions \
    -S sample1.bw -R genes.bed -o matrix.gz \
    --regionBodyLength 1000 -b 500 -a 500
```

**Common options**

| Option | Default | Description |
|--------|---------|-------------|
| `-S, --scoreFileName` | — | Score bigWig file(s) (required) |
| `-R, --regionsFileName` | — | BED/GTF region file(s) (required) |
| `-o, --outFileName` | — | Output `.gz` matrix (required) |
| `--outFileNameMatrix` | — | Raw values as TSV |
| `--outFileSortedRegions` | — | BED of regions in final order |
| `-bs, --binSize` | `10` | Bin size (bases) |
| `--sortRegions` | `keep` | `descend`, `ascend`, `no`, `keep` |
| `--sortUsing` | `mean` | `mean`, `median`, `max`, `min`, `sum`, `region_length` |
| `--sortUsingSamples` | — | Samples (1-indexed) used for sorting |
| `--averageTypeBins` | `mean` | Per-bin summary statistic |
| `--missingDataAsZero` | off | Treat missing data as 0 |
| `--skipZeros` | off | Skip regions that are entirely 0/NaN |
| `--minThreshold` / `--maxThreshold` | — | Exclude regions outside thresholds |
| `-bl, --blackListFileName` | — | Blacklist BED file |
| `--samplesLabel` | — | Sample names (same order as `-S`) |
| `--smartLabels` | off | Infer labels from file names |
| `-p, --numberOfProcessors` | `max` | Integer, `max` or `max/2` |
| `-q, --quiet` / `--verbose` | off | Verbosity |
| `--metagene`, `--transcriptID`, `--exonID`, `--transcript_id_designator` | — | GTF-specific options |

**`reference-point`**: `--referencePoint {TSS,TES,center}` (default `TSS`),
`-b/--beforeRegionStartLength`, `-a/--afterRegionStartLength`, `--nanAfterEnd`.

**`scale-regions`**: `-m/--regionBodyLength` (default `1000`),
`-b`, `-a`, `--unscaled5prime`, `--unscaled3prime`, `--startLabel`, `--endLabel`.

Full help: `python 04_compute_matrix.py reference-point --help` or
`python 04_compute_matrix.py scale-regions --help`.

---

## 5. `05_plotHeatmap.py` — plotHeatmap

Draws a heatmap (with a summary profile on top) from the `computeMatrix`
matrix.

```bash
python 05_plotHeatmap.py \
    -m matrix.gz \
    -out heatmap.png \
    --colorMap RdBu \
    --whatToShow "plot, heatmap and colorbar" \
    --plotTitle "Signal around TSS"
```

Main options (all listed with `--help`):

| Category | Options |
|----------|---------|
| Required | `-m/--matrixFile`, `-out/--outFileName` |
| Outputs | `--outFileNameMatrix`, `--outFileNameData`, `--outFileSortedRegions` |
| Clustering / sorting | `--kmeans`, `--hclust` (mutually exclusive), `--silhouette`, `--linkage`, `--sortRegions`, `--sortUsing`, `--sortUsingSamples` |
| Colors | `--colorMap`, `--colorList`, `--colorNumber`, `--missingDataColor`, `--alpha`, `--zMin`, `--zMax`, `--interpolationMethod` |
| Layout | `--dpi`, `--heatmapHeight`, `--heatmapWidth`, `--whatToShow`, `--boxAroundHeatmaps`, `--perGroup`, `--averageTypeSummaryPlot`, `--numPlotsPerRow` |
| Titles / labels | `--plotTitle`, `--xAxisLabel`, `--yAxisLabel`, `--startLabel`, `--endLabel`, `--refPointLabel`, `--labelRotation`, `--regionsLabel`, `--samplesLabel`, `--legendLocation` |
| Misc | `--yMin`, `--yMax`, `--plotFileFormat`, `--verbose` |

---

## 5bis. `05_bis_plotProfile.py` — plotProfile

Draws the average signal profile over the set of regions.

```bash
python 05_bis_plotProfile.py \
    -m matrix.gz \
    -out profile.png \
    --plotType lines \
    --perGroup \
    --plotTitle "Signal around TSS"
```

Main options:

| Category | Options |
|----------|---------|
| Required | `-m/--matrixFile`, `-out/--outFileName` |
| Outputs | `--outFileNameData`, `--outFileSortedRegions` |
| Clustering | `--kmeans`, `--hclust` (mutually exclusive), `--silhouette`, `--linkage` |
| Computation | `--averageType`, `--nanToZero`, `--proj` |
| Style | `--plotType {lines,fill,se,std,overlapped_lines,heatmap}`, `--colors`, `--plotHeight`, `--plotWidth`, `--dpi`, `--numPlotsPerRow`, `--perGroup` |
| Titles / labels | `--plotTitle`, `--xAxisLabel`, `--yAxisLabel`, `--startLabel`, `--endLabel`, `--refPointLabel`, `--regionsLabel`, `--samplesLabel`, `--legendLocation` |
| Misc | `--yMin`, `--yMax`, `--plotFileFormat`, `--verbose` |

---

## Full pipeline example

```bash
# 1. SGA -> CPM bigWig (1 bp)
python 01_sga_to_bigwig.py data/sga --fasta transcripts.fa
# -> data/bigwig_CPM/*.bw

# 2. Mean of replicates per condition
python 02_bigwig_mean.py -i data/bigwig_CPM/ctrl_rep*.bw -o ctrl_mean.bw --binSize 1
python 02_bigwig_mean.py -i data/bigwig_CPM/treat_rep*.bw -o treat_mean.bw --binSize 1

# 3. Treatment vs control comparison
python 03_bigwig_compare.py -b1 treat_mean.bw -b2 ctrl_mean.bw \
    -o treat_vs_ctrl_log2.bw --operation log2 --binSize 1

# 4. Matrix (regions BED defined on the transcripts from the FASTA)
python 04_compute_matrix.py scale-regions \
    -S treat_vs_ctrl_log2.bw -R transcripts.bed -o matrix.gz \
    --regionBodyLength 1000 -b 100 -a 100 --binSize 10 \
    --startLabel "5'" --endLabel "3'"

# 5. Figures
python 05_plotHeatmap.py -m matrix.gz -out heatmap.pdf --colorMap RdBu
python 05_bis_plotProfile.py -m matrix.gz -out profile.pdf --plotType se
```

---

## Things to watch out for

- **Resolution**: script 01 produces 1 bp bigWig files, but scripts 02 and 03
  use `--binSize 10` by default. To keep the original resolution, pass
  `--binSize 1`.
- **Pseudocounts**: `bigwigCompare` uses `1 1` by default. On CPM values (often
  small per position), a pseudocount of 1 can flatten the signal; adjust with
  `--pseudocount` (e.g. `0.01 0.01`).
- **Strands**: script 01 writes intervals grouped by `(transcript, strand)` into
  a single bedGraph. If an `.sga` file contains both strands for the same
  transcript, the bedGraph is no longer coordinate-sorted and
  `bedGraphToBigWig` will fail. Use one `.sga` per strand (or split strands
  upstream).
- **`-o` vs `-out`**: in `05_plotHeatmap.py` and `05_bis_plotProfile.py`, the
  output option is declared as `-out`/`--outFileName`. The docstring examples
  use `-o`, which only works through `argparse` automatic abbreviation; prefer
  `-out` or `--outFileName`.
- **deepTools compatibility**: the wrappers always forward their default values
  (`--linkage`, `--colorNumber`, etc.). If an option does not exist in your
  deepTools version, the command will fail: upgrade deepTools or remove the
  corresponding option in `build_command`.
- **Regions / FASTA**: transcript names in the regions BED/GTF must match those
  in the FASTA (and therefore in the bigWig files).
