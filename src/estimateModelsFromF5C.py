#!/usr/bin/env python
# author: Jannes Spangenberg
# e-mail: jannes.spangenberg@uni-jena.de
# github: https://github.com/JannesSP
# website: https://jannessp.github.io

from argparse import ArgumentDefaultsHelpFormatter, ArgumentParser, Namespace
from fileio import getFiles, loadFastx, readPolyAStartEnd
from itertools import product
from read5 import read
import pickle
from OnlineMeanVar import LocShift
from os.path import join

PAMODELS = {''.join(motif) : LocShift(30) for motif in product('ACGU', repeat= 5)}
NORMPAMODELS = {''.join(motif) : LocShift(30) for motif in product('ACGU', repeat= 5)}
POLYASTDMODELS = {''.join(motif) : LocShift(30) for motif in product('ACGU', repeat= 5)}

def parse() -> Namespace:
    parser = ArgumentParser(
        formatter_class=ArgumentDefaultsHelpFormatter
    )
    parser.add_argument('--raw', type=str, required=True, help='Path to raw ONT training data')
    parser.add_argument('--fastx', type=str, required=True, help='Basecalls of ONT training data')
    parser.add_argument('--segmentationPickle', type=str, default=None, help='f5c resquiggle segmentation pickle file')
    parser.add_argument('--polya', type=str, required=True, help='Poly A table from nanopolish polya containing the transcript starts')
    parser.add_argument('--out', type=str, required=True, help='Outpath to write files')
    return parser.parse_args()

def main() -> None:
    args = parse()

    rawFiles = getFiles(args.raw, True)
    print(f'ONT Files: {len(rawFiles)}')
    basecalls = loadFastx(args.fastx)
    print(f'Segmenting {len(basecalls)} reads')
    polya = readPolyAStartEnd(args.polya)

    start(rawFiles, basecalls, args.segmentationPickle, polya, args.out)
    print('Done')

def start(files : list, basecalls : dict, segPickle : str, polyA : dict, outdir : str) -> None:
    with open(segPickle, 'rb') as handle:
        borderMap = pickle.load(handle)

    for fi, file in enumerate(files):
        r5 = read(file)

        for ri, readid in enumerate(r5.getReads()):
            if (ri + 1) % 100 == 0:
                print(f'Analysing read {ri + 1} of file {fi}', end='\r')

            if readid not in basecalls:
                continue
            if readid not in borderMap:
                continue
            if readid not in polyA:
                continue

            borders = borderMap[readid]
            basecall = basecalls[readid]

            polyAstart = polyA[readid][0]
            polyAend = polyA[readid][1]

            pASignal = r5.getpASignal(readid)
            normSignal = r5.getZNormSignal(readid)

            if polyAend - polyAstart > 30:
                polyAStdSignal = r5.getPolyAStandardizedSignal(readid, polyAstart, polyAend)
            else:
                polyAStdSignal = None

            for segment in borders:
                try:
                    # IMPORTANT: models are in 5' -> 3' direction
                    motif = basecall[int(segment[2]) - 2 : int(segment[2]) + 3]
                    # change it to 3' -> 5' direction, similar to the sequencing direction
                    motif = motif[::-1].replace('T', 'U')
                except IndexError:
                    continue
                start = int(segment[0])
                end = int(segment[1])
                
                PAMODELS[motif].append(pASignal[start : end])
                NORMPAMODELS[motif].append(normSignal[start : end])
                if polyAStdSignal is not None:
                    POLYASTDMODELS[motif].append(polyAStdSignal[start : end])

    print(f'Done')
    print(f'Writing models to {outdir} ...')

    writeModel(PAMODELS, outdir, 'pAmodels.model')
    writeModel(NORMPAMODELS, outdir, 'normpAmodels.model')
    writeModel(POLYASTDMODELS, outdir, 'polyAstdmodels.model')

def writeModel(models : dict, outdir : str, name : str) -> None:
    with open(join(outdir, name), 'w') as model:
        model.write('kmer\tlevel_mean\tlevel_stdv\n')

        for motif in models:
            mean, stdev = models[motif].meanStdev()
            model.write(f'{motif}\t{mean:.6f}\t{stdev:.6f}\n')

if __name__ == '__main__':
    main()