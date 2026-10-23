#!/usr/bin/env python3
"""
Normalize .sga files to CPM then convert to bigWig at 1bp resolution.

Pipeline:
    .sga  →  CPM normalization  →  bedGraph (temp)  →  bigWig

Usage:
    python sga_to_bigwig.py /path/to/sga_folder --fasta transcripts.fa
    python sga_to_bigwig.py /path/to/sga_folder --fasta transcripts.fa --keep-cpm
"""

import os
import argparse
import subprocess
import tempfile
from collections import defaultdict
from pathlib import Path

import pandas as pd


# ── Helpers ───────────────────────────────────────────────────────────────────

def fasta_to_chrom_sizes(fasta_path: Path, output_path: Path) -> dict:
    """Parse a transcript FASTA, write a chrom.sizes file, and return sizes as dict."""
    chrom_sizes = {}
    current = None
    with open(fasta_path) as fh:
        for line in fh:
            line = line.strip()
            if line.startswith(">"):
                current = line[1:].split()[0]
                chrom_sizes[current] = 0
            elif current:
                chrom_sizes[current] += len(line)
    with open(output_path, "w") as fh:
        for transcript, size in chrom_sizes.items():
            fh.write(f"{transcript}\t{size}\n")
    print(f"  Chrom sizes → {output_path}")
    return chrom_sizes


def normalize_cpm(sga_path: Path) -> pd.DataFrame:
    """Read an SGA file and return a CPM-normalized DataFrame.

    SGA columns: transcript, ref_type, position, strand, count
    The CPM value replaces 'count' in column 5 (index 4).
    """
    df = pd.read_csv(
        sga_path, sep="\t", header=None,
        names=["transcript", "ref_type", "position", "strand", "count"],
    )
    total = df["count"].sum()
    if total == 0:
        return None
    df["count"] = df["count"] * 1e6 / total   # overwrite count with CPM
    return df


def df_to_bedgraph(df: pd.DataFrame, fh, chrom_sizes: dict) -> None:
    grouped = defaultdict(list)
    for row in df.itertuples(index=False):
        grouped[(row.transcript, row.strand)].append((row.position, row.count))

    for (transcript, strand), positions in sorted(grouped.items()):
        positions.sort()
        chrom_size = chrom_sizes.get(transcript)

        # Build clamped intervals, merging duplicates by summing signal
        merged = defaultdict(float)
        for pos, signal in positions:
            start = pos - 1
            end   = pos
            if chrom_size is not None:
                end   = min(end, chrom_size)
                start = min(start, end - 1)
            if start < 0 or start >= end:
                continue
            merged[(start, end)] += signal  # sum signal if coords collide

        for (start, end), signal in sorted(merged.items()):
            fh.write(f"{transcript}\t{start}\t{end}\t{signal}\n")

    fh.flush()


def bedgraph_to_bigwig(bedgraph_path: str, bigwig_path: Path, chrom_sizes_path: Path) -> bool:
    """Call bedGraphToBigWig via subprocess."""
    try:
        result = subprocess.run(
            ["bedGraphToBigWig", bedgraph_path, str(chrom_sizes_path), str(bigwig_path)],
            check=True, capture_output=True, text=True,
        )
        return True
    except subprocess.CalledProcessError as e:
        print(f"    ✗ bedGraphToBigWig error: {e.stderr.strip()}")
        return False
    except FileNotFoundError:
        print("    ✗ bedGraphToBigWig not found in PATH.")
        return False


# ── Main pipeline ─────────────────────────────────────────────────────────────

def process(input_dir: Path, fasta_path: Path, keep_cpm: bool = False) -> None:
    # Output directories
    bigwig_dir = input_dir.parent / "bigwig_CPM"
    bigwig_dir.mkdir(exist_ok=True)

    cpm_dir = input_dir.parent / (input_dir.name + "_normalized_CPM")
    if keep_cpm:
        cpm_dir.mkdir(exist_ok=True)

    print(f"📂 Input   : {input_dir}")
    print(f"💾 BigWig  : {bigwig_dir}")
    if keep_cpm:
        print(f"💾 CPM SGA : {cpm_dir}")

    # Chrom sizes from FASTA — returned as dict for coordinate clamping
    chrom_sizes_path = bigwig_dir / "chrom.sizes"
    chrom_sizes = fasta_to_chrom_sizes(fasta_path, chrom_sizes_path)

    # Process each SGA file
    sga_files = sorted(input_dir.glob("*.sga"))
    if not sga_files:
        print(f"⚠️  No .sga files found in {input_dir}")
        return
    print(f"\nFound {len(sga_files)} .sga file(s)\n")

    for sga_file in sga_files:
        print(f"→ {sga_file.name}")

        # Step 1: CPM normalization
        df = normalize_cpm(sga_file)
        if df is None:
            print(f"  ⚠️  Total counts = 0 → skipped")
            continue
        print(f"  ✅ CPM normalized ({len(df):,} rows)")

        # Optionally save CPM SGA
        if keep_cpm:
            cpm_path = cpm_dir / (sga_file.stem + "_normCPM.sga")
            df.to_csv(cpm_path, sep="\t", index=False, header=False)
            print(f"  💾 CPM SGA saved → {cpm_path.name}")

        # Step 2: bedGraph (temp) → bigWig
        bigwig_path = bigwig_dir / (sga_file.stem + "_CPM.bw")
        with tempfile.NamedTemporaryFile(mode="w", suffix=".bedgraph", delete=False) as tmp:
            df_to_bedgraph(df, tmp, chrom_sizes)
            tmp_path = tmp.name

        try:
            ok = bedgraph_to_bigwig(tmp_path, bigwig_path, chrom_sizes_path)
            if ok:
                print(f"  ✅ BigWig → {bigwig_path.name}")
            else:
                print(f"  ✗  BigWig generation failed for {sga_file.stem}")
        finally:
            os.remove(tmp_path)

    print("\n🎉 All done!")


# ── CLI ───────────────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(
        description="CPM-normalize .sga files and convert to bigWig.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  python sga_to_bigwig.py /path/to/sga_folder --fasta transcripts.fa
  python sga_to_bigwig.py /path/to/sga_folder --fasta transcripts.fa --keep-cpm
        """,
    )
    parser.add_argument("input_dir", help="Directory containing .sga files")
    parser.add_argument("--fasta", required=True, help="Transcript FASTA for chrom sizes")
    parser.add_argument(
        "--keep-cpm", action="store_true",
        help="Also save the CPM-normalized .sga files to disk",
    )
    args = parser.parse_args()
    process(Path(args.input_dir), Path(args.fasta), keep_cpm=args.keep_cpm)


if __name__ == "__main__":
    main()