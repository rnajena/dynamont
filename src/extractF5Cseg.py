#!/usr/bin/env python
# author: Jannes Spangenberg
# e-mail: jannes.spangenberg@uni-jena.de
# github: https://github.com/JannesSP
# website: https://jannessp.github.io

from argparse import ArgumentDefaultsHelpFormatter, ArgumentParser, Namespace
import pandas as pd
import pickle
import numpy as np

def parse() -> Namespace:
    parser = ArgumentParser(
        formatter_class=ArgumentDefaultsHelpFormatter
    )
    parser.add_argument("f5c_result")
    parser.add_argument("outfile")
    switch = parser.add_mutually_exclusive_group(required=True)
    switch.add_argument("--eventalign", action="store_true")
    switch.add_argument("--resquiggle", action="store_true")
    parser.add_argument("--f5c_summary", default=None)
    return parser.parse_args()

def parseSummary(file : str) -> dict:
    df = pd.read_csv(file, sep='\t')
    return pd.Series(df.read_name.values, index=df.read_index).to_dict()

def parseEventalign(eventalign : str, readMap : dict, outfile : str) -> None:
    segmentation = {read : [] for read in readMap.values()} # {readid : [segmentborders]}
    with open(eventalign, 'r') as r:
        r.readline() # skip header
        for lidx, line in enumerate(r):
            if (lidx+1)%100000==0:
                print(f'Line {lidx+1}', end='\r')
            _, _, _, readidx, _, _, _, _, _, _, _, _, _, start, end = line.strip().split('\t')
            segmentation[readMap[int(readidx)]].extend([int(start), int(end)])
        print(f'Line {lidx}')

    for read in segmentation:
        segmentation[read] = sorted(list(set(segmentation[read])))

    with open(outfile, 'wb') as handle:
        pickle.dump(segmentation, handle, pickle.HIGHEST_PROTOCOL)

def parseResquiggle(resquiggle : str, outfile : str) -> None:
    segmentation = {} # {readid : [segmentborders]}
    with open(resquiggle, 'r') as r:
        r.readline() # skip header
        for lidx, line in enumerate(r):
            if (lidx+1)%100000==0:
                print(f'Line {lidx+1}', end='\r')
            readid, kmeridx, start, end = line.strip().split('\t')
            if start == "." or end == ".":
                continue
            if readid not in segmentation:
                segmentation[readid] = []
            segmentation[readid].append([int(start), int(kmeridx)]) #, [int(end), int(kmeridx)]])
        print(f'Line {lidx}')

    for read in segmentation:
        segmentation[read] = np.array(list(set(segmentation[read])))

    with open(outfile, 'wb') as handle:
        pickle.dump(segmentation, handle, pickle.HIGHEST_PROTOCOL)

def main() -> None:
    args = parse()
    print(args)
    if args.eventalign:
        assert args.f5c_summary is not None, f"--f5c_summary needed"
        parseEventalign(args.f5c_result, parseSummary(args.f5c_summary), args.outfile)
    elif args.resquiggle:
        parseResquiggle(args.f5c_result, args.outfile)

if __name__ == '__main__':
    main()