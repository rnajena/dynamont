#!/usr/bin/env python
# author: Jannes Spangenberg
# e-mail: jannes.spangenberg@uni-jena.de
# github: https://github.com/JannesSP
# website: https://jannessp.github.io

from argparse import ArgumentDefaultsHelpFormatter, ArgumentParser, Namespace
import os

def parse() -> Namespace:
    parser = ArgumentParser(
        formatter_class=ArgumentDefaultsHelpFormatter
    )
    parser.add_argument("flye", type=str, help="Path to the flye/assembly_info.txt output")
    parser.add_argument("svim", type=str, help="Path to the svim output")
    parser.add_argument("outfile", type=str, help="Path to the output file")
    return parser.parse_args()

def read_flye(file: str):
    """
    Reads the Flye assembly_info.txt file and calculates:
    - Total length of all contigs
    - N50 fragment length
    - Mean coverage

    Args:
        file (str): Path to the Flye assembly_info.txt file.

    Returns:
        dict: A dictionary containing the total length, N50, and mean coverage.
    """

    if not os.path.exists(file):
        return {
            "total_length": 0,
            "n50": 0,
            "mean_coverage": 0
        }

    contig_lengths = []
    coverages = []

    with open(file, 'r') as f:
        for line in f:
            if line.startswith("#"):  # Skip header
                continue
            parts = line.strip().split("\t")
            length = int(parts[1])  # Contig length
            coverage = float(parts[2])  # Contig coverage

            contig_lengths.append(length)
            coverages.append(coverage)

    # Total length of all contigs
    total_length = sum(contig_lengths)

    # Calculate N50
    sorted_lengths = sorted(contig_lengths, reverse=True)
    cumulative_length = 0
    for length in sorted_lengths:
        cumulative_length += length
        if cumulative_length >= total_length / 2:
            n50 = length
            break

    # Mean coverage
    mean_coverage = sum(coverages) / len(coverages) if coverages else 0

    return {
        "total_length": total_length,
        "n50": n50,
        "mean_coverage": mean_coverage
    }

def count_structural_variants(vcf_file: str) -> int:
    """
    Counts the number of structural variants (SVs) in a VCF file produced by SVIM.

    Args:
        vcf_file (str): Path to the VCF file.

    Returns:
        int: The number of structural variants in the VCF file.
    """
    sv_count = 0

    if not os.path.exists(vcf_file):
        return sv_count
    
    with open(vcf_file, 'r') as f:
        for line in f:
            if line.startswith("#"):  # Skip header lines
                continue
            # Parse the INFO field to check for structural variants
            parts = line.strip().split("\t")
            info_field = parts[7]  # INFO field is the 8th column in VCF
            if "SVTYPE" in info_field:  # Check if the variant is structural
                sv_count += 1

    return sv_count

def main():
    args = parse()
    flye = args.flye
    flye_results = read_flye(flye)

    svim = args.svim
    sv_count = count_structural_variants(svim)

    with open(args.outfile, 'w') as f:
        f.write(f"Total length: {flye_results['total_length']}\n")
        f.write(f"N50: {flye_results['n50']}\n")
        f.write(f"Mean coverage: {flye_results['mean_coverage']}\n")
        f.write(f"Number of structural variants: {sv_count}\n")

if __name__ == '__main__':
    main()