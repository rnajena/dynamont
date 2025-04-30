#!/usr/bin/env python
# author: Jannes Spangenberg
# e-mail: jannes.spangenberg@uni-jena.de
# github: https://github.com/JannesSP
# website: https://jannessp.github.io

from argparse import ArgumentDefaultsHelpFormatter, ArgumentParser, Namespace
import pysam
from pathlib import Path
from statistics import median
import numpy as np
import json
from tqdm import tqdm
from multiprocessing import Pool, cpu_count

def parse() -> Namespace:
    parser = ArgumentParser(
        formatter_class=ArgumentDefaultsHelpFormatter
    )
    parser.add_argument("basecalls", type=str, help="Path to the basecalled bam file")
    parser.add_argument("fasta", type=str, help="Path to the fasta file")
    return parser.parse_args()

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

def load_fasta(fasta_path: str) -> dict[str, str]:
    """
    Loads a fasta file into a dictionary mapping read IDs to sequences.

    Parameters
    ----------
    fasta_path : str
        Path to the fasta file.

    Returns
    -------
    dict[str, str]
        A dictionary mapping read IDs to sequences.
    """
    reads = {}
    with open(fasta_path) as f:
        rid, seq = None, []
        for line in f:
            if line.startswith(">"):
                if rid:
                    reads[rid] = "".join(seq)
                rid = line[1:].strip()
                seq = []
            else:
                seq.append(line.strip())
        if rid:
            reads[rid] = "".join(seq)
    return reads

def global_alignment(a: str, b: str) -> int:
    """Compute edit distance using semi-global Needleman-Wunsch."""
    n, m = len(a), len(b)
    dp = np.zeros((n+1, m+1), dtype=int)
    # Initialization: no penalty for terminal gaps
    dp[0, :] = 0
    dp[:, 0] = 0
    for i in range(1, n+1):
        for j in range(1, m+1): # banded matrix
            match = dp[i-1, j-1] + (a[i-1] != b[j-1])
            delete = dp[i-1, j] + 1
            insert = dp[i, j-1] + 1
            dp[i, j] = min(match, delete, insert)
    return dp[-1, -1]

def _compute_read_stats(args):
    """
    Worker function for a single read.
    Returns (is_present, is_identical, is_truncated, nt_changes, length)
    """
    rid, original_seq, segmented = args
    if rid not in segmented:
        return (False, False, False, None, None)

    seg_seq = segmented[rid]
    is_truncated = len(seg_seq) < len(original_seq)
    nt_changes = global_alignment(original_seq, seg_seq)
    is_identical = not nt_changes
    length = len(seg_seq)

    return (True, is_identical, is_truncated, nt_changes, length)

def compute_stats(basecalls: dict[str, str], segmented: dict[str, str], threads: int = None):
    threads = threads or cpu_count()
    stats = {
        "total": len(basecalls),
        "present": 0,
        "identical": 0,
        "truncated": 0,
        "nt_changed": 0,
        "missing": 0,
        "lengths": [],
    }

    read_args = [(rid, seq, segmented) for rid, seq in basecalls.items()]

    with Pool(threads) as pool:
        results = pool.imap_unordered(_compute_read_stats, read_args, chunksize = 128)
        with tqdm(total=len(read_args), desc="  Computing segmented read stats", unit=" reads") as pbar:
            for present, identical, truncated, nt_changed, length in results:
                pbar.update()
                if not present:
                    stats["missing"] += 1
                    continue

                stats["present"] += 1
                stats["lengths"].append(length)
                stats["identical"] += identical
                stats["truncated"] += truncated
                stats["nt_changed"] += nt_changed

    return stats

def compute_N50(lengths: list[int]) -> int:
    """
    Computes the N50 value from a list of read lengths.

    Parameters:
        lengths (list[int]): A list of integers representing read lengths.

    Returns:
        int: The N50 value, which is the length of the shortest read such that
             the sum of the lengths of all reads of this length or longer
             constitutes at least half of the total sum of all read lengths.
    """

    lengths = sorted(lengths, reverse=True)
    total = sum(lengths)
    cumsum = 0
    for l in lengths:
        cumsum += l
        if cumsum >= total / 2:
            return l
    return 0

def convert_numpy(obj):
    if isinstance(obj, dict):
        return {k: convert_numpy(v) for k, v in obj.items()}
    elif isinstance(obj, list):
        return [convert_numpy(v) for v in obj]
    elif hasattr(obj, "item"):  # NumPy scalar
        return obj.item()
    else:
        return obj

def report_stats(stats: dict, output_path):
    label = output_path.stem
    print(f"\n=== Stats for {label} ===")
    print(f"Total reads in basecalls: {stats['total']}")
    print(f"Present in segmented file: {stats['present']}")
    print(f"Thrown away (missing): {stats['missing']}")
    print(f"Identical reads: {stats['identical']}")
    print(f"Truncated reads: {stats['truncated']}")
    print(f"Total NT changes: {stats['nt_changed']}")
    
    if stats["lengths"]:
        n50 = compute_N50(stats['lengths'])
        mean_len = sum(stats['lengths']) / len(stats['lengths'])
        median_len = median(stats['lengths'])
        max_len = max(stats['lengths'])
        min_len = min(stats['lengths'])
        
        print(f"N50: {n50}")
        print(f"Mean length: {mean_len:.1f}")
        print(f"Median length: {median_len:.1f}")
        print(f"Max length: {max_len}")
        print(f"Min length: {min_len}")
        
        stats.update({
            "N50": n50,
            "mean_length": mean_len,
            "median_length": median_len,
            "max_length": max_len,
            "min_length": min_len,
        })
    else:
        print("No reads to compute length stats.")
    
    # Save to file
    output_file = output_path.with_suffix(".json")
    with open(output_file, "w") as f:
        json.dump(convert_numpy(stats), f, indent=2)

    print(f"\nStats written to {output_file}")


def main():
    args = parse()
    basecalls = load_basecalls_from_bam(args.basecalls)
    segmented = load_fasta(args.fasta)
    stats = compute_stats(basecalls, segmented)
    print("Done computing stats, now writing to file...")
    report_stats(stats, Path(args.fasta))

if __name__ == '__main__':
    main()