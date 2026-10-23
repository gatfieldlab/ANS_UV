SGA-to-BigWig and BigWig averaging
This repository contains two Python scripts for converting .sga signal files into CPM-normalized BigWig tracks and calculating the mean signal across multiple BigWig files.
Scripts
1. 01_sga_to_bigwig.py
Converts each .sga file in an input directory into a CPM-normalized BigWig file at 1-bp resolution.
Workflow
1. Read transcript lengths from a FASTA file
   - Parses the transcript FASTA.
   - Writes a chrom.sizes file containing each transcript identifier and its sequence length.
   - Uses transcript lengths to keep output intervals within valid coordinate bounds.
2. Read and normalize each SGA file
   - Expects five tab-separated columns, with no header:
     transcript, ref_type, position, strand, count.
   - Calculates counts per million (CPM):
     CPM = count × 1,000,000 / sum(all counts in the SGA file)
   - Replaces the original count values with CPM values.
   - Skips files whose total count is zero.
3. Convert signal positions to bedGraph intervals
   - Converts each position to a 1-bp interval: start = position - 1, end = position.
   - Sorts positions by transcript and strand.
   - Clamps intervals to transcript boundaries and sums signal values at duplicate coordinates.
   - Writes a temporary bedGraph file.
4. Convert bedGraph to BigWig
   - Calls the UCSC bedGraphToBigWig command-line utility.
   - Removes the temporary bedGraph file after conversion.
5. Optionally save normalized SGA files
   - With --keep-cpm, saves the CPM-normalized SGA files as well as the BigWig tracks.
Usage
python 01_sga_to_bigwig.py /path/to/sga_folder --fasta /path/to/transcripts.fa
To also save normalized SGA files:
python 01_sga_to_bigwig.py /path/to/sga_folder \
  --fasta /path/to/transcripts.fa \
  --keep-cpm
Outputs
The script creates these directories next to the input directory:
- bigwig_CPM/
  - One BigWig file per input SGA file, named <sample>_CPM.bw.
  - chrom.sizes, generated from the FASTA.
- <input_folder_name>_normalized_CPM/ (only with --keep-cpm)
  - One normalized SGA file per input SGA file, named <sample>_normCPM.sga.
Requirements
- Python 3
- pandas
- UCSC bedGraphToBigWig executable available on PATH
- A transcript FASTA whose sequence identifiers match the transcript identifiers in the SGA files
Input assumptions and notes
- The script assumes that the SGA file's fifth column contains counts and that positions are 1-based.
- CPM normalization is performed separately for each SGA file using the sum of its fifth-column values.
- The SGA format includes a strand column, but the current bedGraph output uses transcript and position coordinates; strand is used to group rows during processing, not encoded in the bedGraph output.
- The FASTA sequence identifiers (the first token after >) must match the transcript names in the SGA files.
- If an SGA file has a total count of zero, it is skipped.
- bedGraphToBigWig must be installed separately; the script reports an error if it is unavailable.
2. 02_bigwig_mean.py
Calculates the mean signal across multiple input BigWig files using deepTools' bigwigAverage.
Workflow
1. Checks that bigwigAverage is available on PATH, unless --skip-check is specified.
2. Verifies that each input file exists.
3. Creates the output directory if necessary.
4. Builds and runs the bigwigAverage command with the selected inputs and options.
5. Reports an error if bigwigAverage exits with a non-zero status.
Usage
Average all .bw files in a directory:
python 02_bigwig_mean.py \
  --input /path/to/bigwig/*.bw \
  --output /path/to/output/mean.bw
Average selected replicates:
python 02_bigwig_mean.py \
  --input rep1.bw rep2.bw rep3.bw \
  --output /path/to/output/mean.bw
Use 1-bp bins and eight processors:
python 02_bigwig_mean.py \
  --input rep1.bw rep2.bw \
  --output /path/to/output/mean.bw \
  --binSize 1 \
  --numberOfProcessors 8
Write bedGraph instead of BigWig:
python 02_bigwig_mean.py \
  --input rep1.bw rep2.bw \
  --output /path/to/output/mean.bedgraph \
  --outFileFormat bedgraph
Arguments
Argument	Description	Default
-i, --input	One or more input BigWig files	Required
-o, --output	Output path	Required
--binSize	Bin size in bases used for averaging	10
--numberOfProcessors, -p	Number of processors	4
--outFileFormat	Output format: bigwig or bedgraph	bigwig
--skipNonCoveredRegions	Skip bins without coverage in all input files	Off
--blackListFileName	Optional BED file of regions to exclude	None
--skip-check	Skip the check for bigwigAverage on PATH	Off
--verbose, -v	Print the command before running it	Off


Requirements
- Python 3
- deepTools, providing bigwigAverage
- Input BigWig files that can be read by deepTools
Typical workflow
For example, if your SGA files are in sga/ and the transcript FASTA is transcripts.fa:
# Step 1: Convert SGA files to CPM-normalized BigWig tracks
python 01_sga_to_bigwig.py sga/ --fasta transcripts.fa

# Step 2: Average selected BigWig tracks
python 02_bigwig_mean.py \
  --input bigwig_CPM/sample1_CPM.bw bigwig_CPM/sample2_CPM.bw \
  --output mean_sample_CPM.bw
The second step should only combine files that use compatible reference sequences and comparable signal units. The mean is calculated by bigwigAverage using the chosen bin size; with the default --binSize 10, the output is binned at 10-bp resolution rather than retaining the original 1-bp resolution.
Troubleshooting
- bedGraphToBigWig not found in PATH: install the UCSC Genome Browser utilities and ensure bedGraphToBigWig is available in your terminal.
- bigwigAverage not found on PATH: activate the environment where deepTools is installed, or install deepTools.
- Input file not found: check the paths passed to --input, --output, and --fasta.
- BigWig files fail to open: test one input file independently before averaging all files. Confirm that the BigWig was written completely and that the transcript names and lengths match the FASTA used for conversion.
