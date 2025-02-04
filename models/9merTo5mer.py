#!/usr/bin/env python
# author: Jannes Spangenberg
# e-mail: jannes.spangenberg@uni-jena.de
# github: https://github.com/JannesSP
# website: https://jannessp.github.io

"""
Using 9mer Gaussian models might lead to overfitting.
Using this script, we reduce the 9mer models to 5mer models by only using the middle 5 bases.
Position-wise conversion:
012345678 -> 23456
"""

from argparse import ArgumentDefaultsHelpFormatter, ArgumentParser, Namespace
import pandas as pd
import numpy as np

def parse() -> Namespace:
    parser = ArgumentParser(
        formatter_class=ArgumentDefaultsHelpFormatter
    )
    parser.add_argument("inmodel", metavar="PATH", help="Path to 9mer model")
    parser.add_argument("outmodel", metavar="PATH", help="Path to write 5mer model")
    return parser.parse_args()

def readKmerModels(filepath : str) -> dict:
    '''
    Returns
    -------
    models : dict
        format: {kmer : (mean, stdev)}
    '''
    df = pd.read_csv(filepath, sep='\t')
    models = pd.Series(zip(df.level_mean.values, df.level_stdv.values), index=df.kmer).to_dict()
    return models

def reduce(inmodel : dict) -> dict:
    outmodel = {}
    for kmer in inmodel:

        rKmer = kmer[2:7]
        if rKmer not in outmodel:
            outmodel[rKmer] = [[], []]

        outmodel[rKmer][0].append(inmodel[kmer][0])
        outmodel[rKmer][1].append(inmodel[kmer][1])

    for kmer in outmodel:
        outmodel[kmer] = [np.mean(outmodel[kmer][0], dtype=np.float128), np.mean(outmodel[kmer][1], dtype=np.float128)]

    return outmodel

def writeModel(model : dict, outpath : str) -> None:
    with open(outpath, 'w') as w:
        w.write('kmer\tlevel_mean\tlevel_stdv\n')
        for kmer in model:
            mean = model[kmer][0]
            stdev = model[kmer][1]
            w.write(f'{kmer}\t{mean}\t{stdev}\n')

def main() -> None:
    args = parse()
    print("Reading model from " + args.inmodel)
    model = readKmerModels(args.inmodel)
    print("Reducing model")
    outmodel = reduce(model)
    print("Writing model to " + args.outmodel)
    writeModel(outmodel, args.outmodel)

if __name__ == '__main__':
    main()