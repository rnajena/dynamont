#!/usr/bin/env python
# author: Jannes Spangenberg
# e-mail: jannes.spangenberg@uni-jena.de
# github: https://github.com/JannesSP
# website: https://jannessp.github.io

from argparse import ArgumentDefaultsHelpFormatter, ArgumentParser, Namespace
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from os.path import exists, join, dirname
from os import makedirs, name
from fileio import getFiles, loadFastx, openCPPScriptParams, stopFeeding, calcZ, readPolyAStartEnd, feedSegmentation
from read5 import read
from hampel import hampel

KMERMODELS = pd.read_csv('/home/yi98suv/projects/dynamont/data/template_median69pA_extended.model', sep='\t', index_col = "kmer")
# KMERMODELS = pd.read_csv('/home/yi98suv/projects/dynamont/data/template_median69pA_reduced.model', sep='\t', index_col = "kmer")

# PARAMS = {
#         "e1":1.,
#         "m2":.04,
#         "d1":.0,
#         "e2":.96,
#         "e3":.0,
#         "i1":.0,
#         "m3":.1,
#         "i2":.0,
#         "m4":.1,
#         "d2":.0,
# }

# E. coli
# PARAMS = {'e1': 1.0, 'm2': 0.062403286666666655, 'd1': 0.0, 'e2': 0.9373050666666667, 'e3': 0.0002915558866666667, 'i1': 0.0, 'm3': 0, 'i2': 0, 'm4': 0, 'd2': 0}
# PARAMS = {'e1': 1.0, 'm2': 0.024169237500000006, 'd1': 4.263524041666667e-06, 'e2': 0.9746092083333332, 'e3': 0.0005571478833333333, 'i1': 0.0006601715291666667, 'm3': 0.8175128750000001, 'i2': 0.182487125, 'm4': 1.0, 'd2': 1.648725e-100}
PARAMS = {'e1': 1.0, 'm2': 0.021304436000000003, 'd1': 0.00052215122, 'e2': 0.9321389999999997, 'e3': 0.044451505, 'i1': 0.00158291295, 'm3': 0.7157543000000001, 'i2': 0.2842457, 'm4': 0.48615565, 'd2': 0.5138443286000001}

# simulated
# PARAMS = {'e1': 1.0, 'm2': 0.025623768750000008, 'd1': 1.950319e-29, 'e2': 0.9743754999999998, 'e3': 1.0193476874999999e-69, 'i1': 8.142517119374999e-07, 'm3': 1.0, 'i2': 6.209504555e-15, 'm4': 1.0, 'd2': 6.981955249999999e-25}

def parse() -> Namespace:
    parser = ArgumentParser(
        formatter_class=ArgumentDefaultsHelpFormatter
    )
    parser.add_argument('--raw', type=str, required=True, help='Path to raw ONT training data')
    parser.add_argument('--fastx', type=str, required=True, help='Basecalls of ONT training data')
    parser.add_argument('--polya', type=str, required=True, help='Poly A table from nanopolish polya containing the transcript starts')
    parser.add_argument('--out', type=str, required=True, help='Outpath to write files')
    parser.add_argument('--readid', type=str, required=True, help='Read to plot')
    parser.add_argument('--resquigglePickle', type=str, default=None, help='f5c resquiggle segmentation pickle file')
    parser.add_argument('--eventalignPickle', type=str, default=None, help='f5c eventalign segmentation pickle file')
    return parser.parse_args()

def plotBorders(signal : np.ndarray, hampel_signal : np.ndarray, polyAend : int, read : str, segments : np.ndarray, readid : str, outpath : str, resquiggleBorders : np.ndarray, eventalignBorders : np.ndarray):
    '''
    Input
    -----
    segments : np.ndarray
        in 3' -> 5' orientation
    read : str
        in 3' -> 5' orientation
    '''
    nPlots = 1
    if resquiggleBorders is not None:
        nPlots += 1
    if eventalignBorders is not None:
        nPlots += 1

    fig, ax = plt.subplots(nrows = nPlots, sharex=True, sharey=True, figsize=(110,12), dpi=300)
    fig.suptitle(f'{readid} segmentation in 3\' -> 5\' orientation')
    x=np.arange(-polyAend, len(signal) - polyAend, 1)
    for axis in range(nPlots):
        ax[axis].plot(x, signal, alpha=0.4, color='black', label='signal', linewidth=1)
        ax[axis].plot(x, hampel_signal, color='purple', label='hampel(signal)', linewidth=0.5)
        ax[axis].set_ylim((min(signal)-15, max(signal)))
        ax[axis].set_ylabel('Signal pico Ampere')
        ax[axis].set_xticks(np.arange(0, len(signal), 2000))
        ax[axis].grid(True, 'both', 'both')

    plotNumber = 0
    # F5C EVENTALIGN
    if eventalignBorders is not None:
        ax[plotNumber].vlines(eventalignBorders[:, 0].astype(int), ymin=min(signal)-5, ymax=max(signal), colors='darkgreen', linestyles='--', label=f'f5c eventalign segmentation', linewidth=0.7)
        for border in eventalignBorders:
            # border indices are for read 5' 3' direction
            base = border[2]
            ax[plotNumber].text(int(border[0]) + np.abs(int(border[1]) - int(border[0]))/2 - 3, min(signal)-5, base, fontdict={'size' : 7, 'color':'darkgreen'})
        for i, segment in enumerate(eventalignBorders):
            motif = ''
            for j in range(i-2, i+3, 1):
                if j < 0:
                    motif += 'A'
                elif j >= len(eventalignBorders):
                    motif += 'N'
                else:
                    base = eventalignBorders[j][2]
                    motif += base
            motif = motif.replace('U', 'T')
            x = int(segment[0])
            width = int(segment[1]) - int(segment[0])
            mean, stdev = KMERMODELS.loc[motif][['level_mean', 'level_stdv']]
            height = 2*stdev
            ax[plotNumber].add_patch(Rectangle((x, mean-stdev), width, height, alpha=0.4, facecolor="yellow"))
        ax[plotNumber].legend()
        plotNumber += 1

    # F5C RESQUIGGLE
    if resquiggleBorders is not None:
        ax[plotNumber].vlines(resquiggleBorders[:, 0].astype(int), ymin=min(signal)-5, ymax=max(signal), colors='darkgreen', linestyles='--', label=f'f5c resquiggle segmentation', linewidth=0.7)
        for border in resquiggleBorders:
            # border indices are for read 5' 3' direction
            base = read[::-1][border[2]]
            ax[plotNumber].text(int(border[0]) + np.abs(int(border[1]) - int(border[0]))/2 - 3, min(signal)-5, base, fontdict={'size' : 7, 'color':'darkgreen'})
        for i, segment in enumerate(resquiggleBorders):
            motif = ''
            for j in range(i-2, i+3, 1):
                if j < 0:
                    motif += 'A'
                elif j >= len(resquiggleBorders):
                    motif += 'N'
                else:
                    base = resquiggleBorders[j][2]
                    motif += read[::-1][base]
            motif = motif.replace('U', 'T')
            x = segment[0]
            width = segment[1] - segment[0]
            mean, stdev = KMERMODELS.loc[motif][['level_mean', 'level_stdv']]
            height = 2*stdev
            ax[plotNumber].add_patch(Rectangle((x, mean-stdev), width, height, alpha=0.4, facecolor="yellow"))
        ax[plotNumber].legend()
        plotNumber += 1

    # OUR SEGMENTATION
    ax[plotNumber].vlines(segments[:, 0], ymin=min(signal)-10, ymax=max(signal), colors='red', linestyles='--', label='our segmentation', linewidth=0.7, alpha=0.6)
    # plot kmer model range
    for segment in segments:
        
        # extract motif from read
        pos = segment[2]
        # add As to 3' end
        if pos-2 < 0:
            motif = ('A' * abs(pos-2)) + read[max(0, pos-2) : pos+3]
        # add Ns to 5' end
        elif pos+3 > len(read):
            motif = read[pos-2 : pos+3] + 'N'*(pos+3-len(read))
        else:
            motif = read[pos-2 : pos+3]
        # plot motif in 3' - 5' direction
        motif = motif.replace('U', 'T')
        x = segment[0]
        width = segment[1] - segment[0]

        if segment[3] == 'M':
            # draw kmer range as rectangle
            mean, stdev = KMERMODELS.loc[motif][['level_mean', 'level_stdv']]
            height = 2*stdev
            ax[plotNumber].add_patch(Rectangle((x, mean-stdev), width, height, alpha=0.4, facecolor="yellow"))
            # write motif
            ax[plotNumber].text(x + width/2 - 3, min(signal)-10, read[pos], fontdict={'size' : 7, 'color':'red'})
        elif segment[3] == 'D':
            height = max(signal) - min(signal)
            ax[plotNumber].add_patch(Rectangle((x, min(signal)), width, height, alpha=0.3, facecolor="grey"))
            ax[plotNumber].text(x + width/2 - 3, min(signal)-14, 'Del', rotation=90, fontdict={'size' : 6, 'color' : 'grey'})
        elif segment[3] == 'I':
            ax[plotNumber].text(x - 3, min(signal)-14, 'Ins', rotation=90, fontdict={'size' : 6, 'color' : 'grey'})
        ax[plotNumber].legend()

    # plt.legend()
    plt.savefig(join(outpath, readid + '.svg'))
    plt.savefig(join(outpath, readid + '.pdf'))

def segmentRead(standardizedSignal : np.ndarray, polyAend : int, read : str, readid : str, outdir : str, resquiggleBorders : np.ndarray, eventalignBorders : np.ndarray):
    '''
    Takes read in 3' -> 5' direction
    '''

    # Create pipe between python script and cp segmentation script
    CPP_SCRIPT = join(dirname(__file__), 'segmentation')
    if name == 'nt': # check for windows
        CPP_SCRIPT+='.exe'
    pipe = openCPPScriptParams(CPP_SCRIPT, PARAMS)

    hampel_std_signal = hampel(standardizedSignal, 20, 2.).filtered_data
    # hampel_raw_signal = hampel(rawSignal, 20, 2.).filtered_data
    # segments, probs = feedSegmentation(hampel_std_signal[polyAstart:], read, pipe)
    segments, probs = feedSegmentation(hampel_std_signal, read, pipe)
    # print(segments)
    segments[:, 0] = segments[:, 0] - polyAend
    segments[:, 1] = segments[:, 1] - polyAend
    stopFeeding(pipe)

    with open(join(outdir, readid + '.txt'), 'w') as w:
        w.write('\n'.join(list(map(str, segments))))

    # print(len(segments), segments)
    # print(len(read), read)

    plotBorders(standardizedSignal, hampel_std_signal, polyAend, read, segments, readid, outdir, resquiggleBorders, eventalignBorders)
    print(calcZ(standardizedSignal, read, PARAMS, CPP_SCRIPT))

def start(files, basecalls, targetID, polyA, out, resquigglePickle, eventalignPickle) -> tuple:
    for file in files:
        r5 = read(file)
        if targetID in r5.getReads():
            
            import pickle
            if resquigglePickle:
                with open(resquigglePickle, 'rb') as handle:
                    borderMap = pickle.load(handle)
                    if targetID in borderMap:
                        resquiggleBorders = borderMap[targetID]
                        try:
                            resquiggleBorders[:, 0] = resquiggleBorders[:, 0].astype(int) - polyA[targetID][1]
                            resquiggleBorders[:, 1] = resquiggleBorders[:, 1].astype(int) - polyA[targetID][1]
                        except:
                            print(resquiggleBorders)
                            exit(1)
                        # print(borders)
                    else:
                        resquiggleBorders = None
                        print(f'WARNING: no border found in {handle}')

            if eventalignPickle:
                with open(eventalignPickle, 'rb') as handle:
                    borderMap = pickle.load(handle)
                    if targetID in borderMap:
                        eventalignBorders = borderMap[targetID]
                        try:
                            eventalignBorders[:, 0] = eventalignBorders[:, 0].astype(int) - polyA[targetID][1]
                            eventalignBorders[:, 1] = eventalignBorders[:, 1].astype(int) - polyA[targetID][1]
                        except:
                            print(eventalignBorders)
                            exit(1)
                        # print(borders)
                    else:
                        eventalignBorders = None
                        print(f'WARNING: no border found in {handle}')

            # change read from 5'-3' to 3'-5'
            # signal = r5.getpASignal(targetID)[polyA[targetID]:]
            print('polyA:', polyA[targetID])
            if polyA[targetID][1] - polyA[targetID][0] > 30:
                signal = r5.getPolyAStandardizedSignal(targetID, polyAstart=polyA[targetID][0], polyAend=polyA[targetID][1])
            else:
                signal = r5.getpASignal(targetID)
            segmentRead(signal, polyA[targetID][1], basecalls[targetID][::-1], targetID, out, resquiggleBorders, eventalignBorders)

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
    start(rawFiles, basecalls, args.readid, polya, args.out, args.resquigglePickle, args.eventalignPickle)

if __name__ == '__main__':
    main()