#!/usr/bin/env python
# author: Jannes Spangenberg
# e-mail: jannes.spangenberg@uni-jena.de
# github: https://github.com/JannesSP
# website: https://jannessp.github.io

from argparse import ArgumentDefaultsHelpFormatter, ArgumentParser, Namespace
import pysam

def parse() -> Namespace:
    parser = ArgumentParser(
        formatter_class=ArgumentDefaultsHelpFormatter
    )
    parser.add_argument('samfiles', nargs='+', help='List of basecalled dorado files that contain the mv:B:c tag.')
    parser.add_argument('-o', '--outfile', default='moves.tsv', help='output file')
    return parser.parse_args()

def extractMoves(bamFile : str, outfile : str, writeBatchSize=1000):
    """Extracts the positions where mv:B:c tag has 1s and maps them to base positions for each read."""
    
    # Open the SAM file using pysam
    with open(outfile, 'a') as out:
        # Open the SAM or BAM file
        with pysam.AlignmentFile(bamFile, "r" if bamFile.endswith('.sam') else "rb", check_sq=False) as bamFile:
            batch = []
            for read in bamFile.fetch(until_eof=True):
                # Check if the mv tag exists in the read
                if read.has_tag("mv"):
                    mv = read.get_tag("mv")
                    ts = read.get_tag("ts")
                    sl = read.get_tag("ns")
                    seq = read.query_sequence
                    rid = read.query_name
                    sid = read.get_tag("pi") if read.has_tag("pi") else rid

                    # mv:B:c is stored as a list of integers, extract the first as the downsample value
                    scale = mv[0]

                    # The rest of the mv tag contains the binary transitions
                    # Find the indices where there is a '1' (transition points)
                    # Map the transition positions to base positions using the scale
                    transitions = ([ts+i*scale for i, value in enumerate(mv[1:]) if value == 1] + [sl])[::-1]

                    # assert len(transitions) - 1 == len(seq)

                    # Iterate over each transition and extract the required data
                    for pos in range(len(seq)-1, -1, -1):
                        # start = sl - 1 - transitions[pos + 1]
                        # end = sl - 1 - transitions[pos]
                        start = transitions[pos + 1]
                        end = transitions[pos]
                        base = seq[pos]
                        # Extract 5-mer motif around the transition position (from pos-2 to pos+3)
                        motif_start = max(pos - 2, 0)
                        motif_end = min(pos + 3, len(seq))
                        motif = seq[motif_start:motif_end]
                        batch.append('\t'.join([rid, sid, str(pos), base, motif, str(start), str(end)]) + '\n')

                        # Write in batches to minimize file I/O
                        if len(batch) >= writeBatchSize:
                            out.writelines(batch)
                            batch.clear()
                    
                    # break

            # Write any remaining rows
            if batch:
                out.writelines(batch)

def main() -> None:
    args = parse()
    files = args.samfiles

    # Write the header if it's a new file
    outfile = args.outfile
    with open(outfile, 'w') as out:
        out.write('\t'.join(['readid', 'signalid', 'position', 'base', 'motif', 'start', 'end']) + '\n')

    for fi, file in enumerate(files):
        print(f'Extracting segmentation from file {fi}/{len(files)}', end='\r')
        extractMoves(file, outfile)

if __name__ == '__main__':
    main()