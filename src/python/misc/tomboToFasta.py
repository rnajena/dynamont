#!/usr/bin/env python
# author: Jannes Spangenberg
# e-mail: jannes.spangenberg@uni-jena.de
# github: https://github.com/JannesSP
# website: https://jannessp.github.io

from argparse import ArgumentDefaultsHelpFormatter, ArgumentParser, Namespace
import os
from pathlib import Path
from collections import defaultdict
from tqdm import tqdm
import h5py

def parse() -> Namespace:
    parser = ArgumentParser(
        formatter_class=ArgumentDefaultsHelpFormatter
    )
    parser.add_argument("segmentation", type=str, help="Path to tombo single fast5 directory")
    parser.add_argument("outfile", type=str, help="Path to the output file")
    parser.add_argument("--rna", action="store_true", help="RNA mode")
    return parser.parse_args()

def collect_fast5_files(seg_dir: str) -> list[Path]:
    """
    Recursively collect all .fast5 files from the given directory.

    Parameters
    ----------
    seg_dir : str
        Path to the directory to search for .fast5 files.

    Returns
    -------
    list[Path]
        A list of Path objects representing the .fast5 files.
    """
    seg_path = Path(seg_dir)
    if not seg_path.is_dir():
        raise NotADirectoryError(f"{seg_dir} is not a valid directory.")
    
    # Collect all .fast5 files recursively
    fast5_files = [file.resolve() for file in seg_path.rglob("*.fast5")]
    return fast5_files

def convert_to_fasta(seg_dir: str, rna: bool, outfile: str):
    fast5s = collect_fast5_files(seg_dir)
    out = Path(outfile).with_suffix(".fasta")

    segments = defaultdict(list)

    for file in tqdm(fast5s, desc="Processing .fast5 files", unit=" file"):
        readid = Path(file).stem

        with h5py.File(file, "r") as f5:
            try:
                # Extract the Events dataset
                # print(f5.keys())
                events = f5['Analyses/RawGenomeCorrected_000/BaseCalled_template/Events'][:]
                # print(f"ReadID: {readid}, Events: {len(events)} entries")
            except KeyError:
                # print(f"Warning: Dataset not found in {file}")
                continue

        # Process the events or store them as needed
        for event in events:
            start = event["start"]
            base = event["base"].decode("utf-8").replace("U", "T")
            segments[readid].append((start, base))

    # Write the output FASTA file
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
    outfile = args.outfile
    rna = args.rna

    if not os.path.exists(segmentation):
        raise FileNotFoundError(f"Segmentation file {segmentation} does not exist.")
    
    if not os.path.exists(os.path.dirname(outfile)) and os.path.dirname(outfile):
        os.makedirs(os.path.dirname(outfile))
    
    convert_to_fasta(segmentation, rna, outfile)

if __name__ == '__main__':
    main()