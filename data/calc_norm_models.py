#!/usr/bin/env python
# author: Jannes Spangenberg
# e-mail: jannes.spangenberg@uni-jena.de
# github: https://github.com/JannesSP
# website: https://jannessp.github.io

from argparse import ArgumentDefaultsHelpFormatter, ArgumentParser, Namespace
import numpy as np
from read5 import read
from itertools import product
from scipy.stats import median_abs_deviation as mad

def parse() -> Namespace:
    parser = ArgumentParser(
        formatter_class=ArgumentDefaultsHelpFormatter
    )
    parser.add_argument('raw', type=str, help='Raw ONT training data')
    parser.add_argument('dynamont', type=str, help='Dynamont segmentation')
    parser.add_argument('outfile', type=str, help='Outfile with models in 3\' -> 5\' orientation')
    return parser.parse_args()

def readDynamont(file : str) -> dict:
    res = {}
    print(f"Reading {file}")
    with open(file, 'r') as r:
        r.readline()
        for line in r:
            readid,start,end,basepos,base,motif,state = line.strip().split(',')
            # skip short motifs
            if len(motif) < 5:
                continue
            if not readid in res:
                res[readid] = []
            res[readid].append([int(start), int(end), motif.replace('U', 'T')])
    return res

def extractSegments(normSignal : np.ndarray, segments : list, models : dict) -> dict:
    for segment in segments:
        start, end, motif = segment
        # skip short segments
        if end-start < 20:
            continue
        models[motif].extend(normSignal[start+5:end-5])
    return models

def writeModels(models : dict, outfile : str) -> None:
    print(f"Writing {len(models)} models to {outfile}")
    with open(outfile + "_mean", 'w') as w:
        w.write('kmer\tlevel_mean\tlevel_stdv\n')
        for motif in [''.join(motif) for motif in list(product('ACGT', repeat=5))]:
            w.write(f'{motif}\t{np.mean(models[motif[::-1]]):.5f}\t{np.std(models[motif[::-1]]):.5f}\n')

    with open(outfile + "_median", 'w') as w:
        w.write('kmer\tlevel_mean\tlevel_stdv\n')
        for motif in [''.join(motif) for motif in list(product('ACGT', repeat=5))]:
            w.write(f'{motif}\t{np.median(models[motif[::-1]]):.5f}\t{mad(models[motif[::-1]]):.5f}\n')


def main() -> None:
    args = parse()
    r5 = read(args.raw)
    dynamont = readDynamont(args.dynamont)

    models = {''.join(motif) : [] for motif in list(product('ACGT', repeat=5))}
    print("Extracting segments")
    for rnum, readid in enumerate(dynamont):
        if not readid in r5.getReads():
            continue
        if (rnum+1)%1000==0:
            print(f"Extracting from read {rnum+1}", end='\r')
        if (rnum+1)%10000==0:
            writeModels(models, args.outfile + f"_{rnum+1}")
        models = extractSegments(r5.getZNormSignal(readid, "mean"), dynamont[readid], models)
    writeModels(models, args.outfile)

if __name__ == '__main__':
    main()