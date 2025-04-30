#!/usr/bin/env python
# author: Jannes Spangenberg
# e-mail: jannes.spangenberg@uni-jena.de
# github: https://github.com/JannesSP
# website: https://jannessp.github.io

from argparse import ArgumentDefaultsHelpFormatter, ArgumentParser, Namespace
import os
import csv
from pathlib import Path
from collections import defaultdict
import pysam
import tqdm

def parse() -> Namespace:
    parser = ArgumentParser(
        formatter_class=ArgumentDefaultsHelpFormatter
    )
    parser.add_argument("segmentation", type=str, help="Path to the segmentation file")
    parser.add_argument("basecalls", type=str, help="Path to the original basecalled bam file")
    parser.add_argument("outfile", type=str, help="Path to the output file")
    parser.add_argument("--rna", action="store_true", help="RNA mode")
    return parser.parse_args()

def detect_format(header: list[str]) -> str:
    """Determine the format of the segmentation file based on its header."""
    header_set = set(h.strip().lower() for h in header)

    if {"readid", "motif"}.issubset(header_set):
        return "dynamont"
    elif {"read_id", "kmer_idx", "start_raw_idx", "end_raw_idx"}.issubset(header_set):
        return "f5c_resquiggle"
    elif {"model_kmer", "read_index", "start_idx", "end_idx"}.issubset(header_set):
        return "f5c_eventalign"
    elif {"aln.read_id", "seq.kmer", "dtw.start", "dtw.length"}.issubset(header_set):
        return "uncalled4"
    else:
        raise ValueError(f"Unsupported or unrecognized segmentation file format: {header}")

def load_read_map(summary_file: Path) -> dict:
    """Loads f5c_eventalign .sum file mapping read_index to read_name."""
    read_map = {}
    with summary_file.open() as f:
        for line in f:
            if line.startswith("read_index"):  # skip header
                continue
            cols = line.strip().split('\t')
            if len(cols) >= 2:
                read_map[cols[0]] = cols[1]
    return read_map

def load_basecalls_from_bam(basecalls: str) -> dict[str, str]:
    """
    Loads basecalled sequences from a BAM file produced by Dorado.

    Parameters
    ----------
    bam_path : str
        Path to the basecalled.bam file.

    Returns
    -------
    dict[str, str]
        A dictionary mapping read IDs to basecalled nucleotide sequences.
    """
    reads = {}
    with pysam.AlignmentFile(basecalls, "r" if basecalls.endswith('.sam') else "rb", check_sq=False) as bam:
        for read in bam.fetch(until_eof=True):
            reads[read.query_name] = read.query_sequence
    return reads

def convert_to_fasta(seg_file: str, basecalls: str, rna: bool, outfile: str):
    seg_path = Path(seg_file)
    out = Path(outfile).with_suffix(".fasta")

    # Load read sequences
    # print(f"Loading reads from: {basecalls}")
    reads = load_basecalls_from_bam(basecalls)

    # Detect format
    with seg_path.open() as f:
        header = f.readline().strip()
        delimiter = "\t" if '\t' in header else ','
        format_type = detect_format(header.split(delimiter))

    # print(f"Detected format: {format_type}")

    segments = defaultdict(list)

    read_map = {}
    if format_type == "f5c_eventalign":
        summary_file = seg_path.with_suffix('.sum')
        # print(f"Loading read map from: {summary_file}")
        read_map = load_read_map(summary_file)

    with seg_path.open() as f, tqdm.tqdm(desc="Translating segments to fasta", unit=" lines") as pbar:
        reader = csv.DictReader(f, delimiter=delimiter)
        for row in reader:
            match format_type:
                case "dynamont":
                    rid = row["readid"]
                    start = int(row["start"])
                    base = row["base"].replace("U", "T")
                case "f5c_resquiggle":
                    rid = row["read_id"]
                    kmer_idx = int(row["kmer_idx"])
                    seq = reads.get(rid)
                    base = seq[kmer_idx].replace("U", "T")
                    try:
                        start = int(row["start_raw_idx"])
                    except ValueError:
                        continue
                case "f5c_eventalign":
                    rid_index = row["read_index"]
                    rid = read_map.get(rid_index)
                    try:
                        start = int(row["start_idx"])
                    except ValueError:
                        continue
                    motif = row["model_kmer"].replace("U", "T")
                    base = motif[len(motif) // 2]
                case "uncalled4":
                    rid = row["aln.read_id"]
                    try:
                        start = int(row["dtw.start"])
                    except ValueError:
                        continue
                    motif = row["seq.kmer"].replace("U", "T")
                    base = motif[len(motif) // 2]
                case _:
                    raise ValueError("Unsupported format")
            segments[rid].append((start, base))
            pbar.update(1)

    # print(f"Writing output to: {out}")
    with out.open("w") as fasta:
        for rid, segs in segments.items():
            segs.sort()
            seq = "".join(m[1] for m in segs)
            if rna:
                seq = seq[::-1]
            fasta.write(f">{rid}\n{seq}\n")

def main() -> None:
    args = parse()
    segmentation = args.segmentation
    basecalls = args.basecalls
    outfile = args.outfile
    rna = args.rna

    if not os.path.exists(segmentation):
        raise FileNotFoundError(f"Segmentation file {segmentation} does not exist.")
    
    if not os.path.exists(basecalls):
        raise FileNotFoundError(f"Basecalls file {basecalls} does not exist.")
    
    if not os.path.exists(os.path.dirname(outfile)) and os.path.dirname(outfile):
        os.makedirs(os.path.dirname(outfile))
    
    convert_to_fasta(segmentation, basecalls, rna, outfile)

if __name__ == '__main__':
    main()