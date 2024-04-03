#!/usr/bin/env python
# author: Jannes Spangenberg
# e-mail: jannes.spangenberg@uni-jena.de
# github: https://github.com/JannesSP
# website: https://jannessp.github.io

import sys
sys.path.append('/home/yi98suv/projects/dynamont/src')

from argparse import ArgumentDefaultsHelpFormatter, ArgumentParser, Namespace
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from os.path import exists, join
from os import makedirs
from dynamont.FileIO import getFiles, loadFastx, readPolyAStartEnd
from read5 import read
import pickle
from hampel import hampel

KMERMODELS = pd.read_csv('/home/yi98suv/projects/dynamont/data/template_median69pA_extended.model', sep='\t', index_col = "kmer")

def parse() -> Namespace:
    parser = ArgumentParser(
        formatter_class=ArgumentDefaultsHelpFormatter
    )
    parser.add_argument('--raw', type=str, required=True, help='Path to raw ONT training data')
    parser.add_argument('--fastx', type=str, required=True, help='Basecalls of ONT training data')
    parser.add_argument('--polya', type=str, required=True, help='Poly A table from nanopolish polya containing the transcript starts')
    parser.add_argument('--out', type=str, required=True, help='Outpath to write files')
    parser.add_argument('--readid', type=str, required=True, help='Read to plot')
    parser.add_argument('--seg_pickle', type=str, required=True, help='f5c segmentation pickle file')
    return parser.parse_args()

def plotBorders(signal : np.ndarray, hampel_signal : np.ndarray, polyAend : int, read : str, readid : str, outpath : str, borders : np.ndarray = None):
    '''
    Input
    -----
    segments : np.ndarray
        in 3' -> 5' orientation
    read : str
        in 3' -> 5' orientation
    '''

    getBase = {True : lambda x : read[::-1][x], False: lambda x : x}

    fig, ax = plt.subplots(figsize=(110,8), dpi=300)
    fig.suptitle(f'{readid} segmentation in 3\' -> 5\' orientation')
    x=np.arange(-polyAend, len(signal) - polyAend, 1)
    ax.plot(x, signal, alpha=0.4, color='black', label='signal', linewidth=1)
    ax.plot(x, hampel_signal, color='purple', label='hampel(signal)', linewidth=0.5)
    ax.set_ylim((min(signal)-15, max(signal)))
    ax.set_ylabel('Signal pico Ampere')
    ax.set_xticks(np.arange(0, len(signal), 2000))

    resquiggle = np.issubdtype(type(borders[0][2]), np.integer)

    if borders is not None:
        label = 'resquiggle' if resquiggle else 'eventalign'
        ax.vlines(borders[:, 0].astype(int), ymin=min(signal)-5, ymax=max(signal), colors='darkgreen', linestyles='--', label=f'f5c {label} segmentation', linewidth=0.7)

        for border in borders:
            # border indices are for read 5' 3' direction
            base = getBase[resquiggle](border[2]) # read[::-1][border[2]] if resquiggle else border[2]
            plt.text(int(border[0]) + np.abs(int(border[1]) - int(border[0]))/2 - 3, min(signal)-5, base, fontdict={'size' : 8, 'color':'darkgreen'})
            
    # plot kmer model range
    for i, segment in enumerate(borders):

        motif = ''
        for j in range(i-2, i+3, 1):
            if j < 0:
                motif += 'A'
            elif j >= len(borders):
                motif += 'N'
            else:
                base = borders[j][2]
                motif += getBase[resquiggle](base)
        motif = motif.replace('U', 'T')

        x = segment[0]
        width = segment[1] - segment[0]
        mean, stdev = KMERMODELS.loc[motif][['level_mean', 'level_stdv']]
        height = 2*stdev
        ax.add_patch(Rectangle((x, mean-stdev), width, height, alpha=0.4, facecolor="yellow"))
        # plt.text(x + width/2 - 3, min(signal)-10, getBase[resquiggle](segment[2]), fontdict={'size' : 8, 'color':'red'})

    plt.legend()
    plt.grid(True, 'both', 'both')
    plt.savefig(join(outpath, readid + '.svg'))
    plt.savefig(join(outpath, readid + '.pdf'))

def start(files, basecalls, targetID, polyA, outdir, seg_pickle) -> tuple:
    for file in files:
        r5 = read(file)
        if targetID in r5.getReads():
            # if seg_pickle:
            with open(seg_pickle, 'rb') as handle:
                borderMap = pickle.load(handle)
                if targetID in borderMap:
                    borders = borderMap[targetID]
                    try:
                        borders[:, 0] = borders[:, 0].astype(int) - polyA[targetID][1]
                        borders[:, 1] = borders[:, 1].astype(int) - polyA[targetID][1]
                    except:
                        print(borders)
                        exit(1)
                    # print(borders)
                else:
                    borders = None
                    print(f'WARNING: no border found in {handle}')
            # signal = r5.getpASignal(targetID)[polyA[targetID]:]
            if polyA[targetID][1] - polyA[targetID][0] > 30:
                signal = r5.getPolyAStandardizedSignal(targetID, polyAstart=polyA[targetID][0], polyAend=polyA[targetID][1])
            else:
                signal = r5.getpASignal(targetID)
            hampel_std_signal = hampel(signal, 20, 2.).filtered_data
            plotBorders(signal, hampel_std_signal, polyA[targetID][1], basecalls[targetID][::-1], targetID, outdir, borders)

def main() -> None:
    args = parse()
    if not exists(args.out):
        makedirs(args.out)
    polya = readPolyAStartEnd(args.polya)
    rawFiles = getFiles(args.raw, True)
    print(f'ONT Files: {len(rawFiles)}')
    basecalls = loadFastx(args.fastx)
    # print("5' -> 3'", len(basecalls[args.readid]), basecalls[args.readid].replace("U", "T"))
    print(f'Segmenting {len(basecalls)} reads')
    start(rawFiles, basecalls, args.readid, polya, args.out, args.seg_pickle)

if __name__ == '__main__':
    main()