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
from os.path import exists, join, dirname
from os import makedirs, name
from dynamont.FileIO import feedSegmentation, SegmentationError
import read5
from hampel import hampel
import pickle
import pysam

STANDARDONTMODEL = pd.read_csv("/home/yi98suv/projects/dynamont/data/template_median69pA_extended.model", sep='\t', index_col = "kmer")

def parse() -> Namespace:
    parser = ArgumentParser(
        formatter_class=ArgumentDefaultsHelpFormatter
    )
    parser.add_argument('-r', '--raw',   type=str, required=True, metavar="POD5|FAST5|SLOW5", help='Raw ONT training data')
    parser.add_argument('-b', '--basecalls', type=str, required=True, metavar="BAM", help='Basecalls of ONT training data as .bam file')
    parser.add_argument('-o', '--outdir',   type=str, required=True, metavar="PATH", help='Outpath to write files')
    parser.add_argument('--readid', type=str, required=True, help='Read to plot')
    parser.add_argument('--pore',  type=str, required=True, choices=["rna_r9", "dna_r9", "rna_rp4", "dna_r10_260bps", "dna_r10_400bps"], help='Pore generation used to sequence the data')
    # optional
    parser.add_argument('--resquigglePickle', type=str, default=None, help='f5c resquiggle segmentation pickle file')
    parser.add_argument('--eventalignPickle', type=str, default=None, help='f5c eventalign segmentation pickle file')
    parser.add_argument('--mode', choices=['basic', 'banded', 'resquiggle'], required=True)
    parser.add_argument('--model_path', type=str, default=join(dirname(__file__), '..', '..', 'data', 'norm_models', 'rna_r9.4_180mv_70bps.model'), help='Kmer model file')
    parser.add_argument('--minSegLen', type=int, default=1, help='Minmal allowed segment length')
    parser.add_argument('--probability', action="store_true", help="Output the segment border probability per position.")
    
    return parser.parse_args()

def plotBorders(signal : np.ndarray, normSignal : np.ndarray, polyAend : int, read : str, segments : np.ndarray, probs : np.ndarray, readid : str, outpath : str, resquiggleBorders : np.ndarray, eventalignBorders : np.ndarray, kmermodels : pd.DataFrame, probability : bool):
    '''
    Input
    -----
    segments : np.ndarray
        in 5' -> 3' orientation, but readpos is on 5' -> 3' orientation
        [[start : int, end : int, readpos : int, state : str] ...]
        already shifted by polyAend
    resquiggleBorders : np.ndarray
        already shifted by polyAend
    eventalignBorders : np.ndarray
        already shifted by polyAend
    read : str
        in 5' -> 3' orientation
    '''
    basecolors = {
        'A':'#377eb8',
        'a':'#377eb8',
        'C':'#ff7f00',
        'c':'#ff7f00',
        'G':'#4daf4a',
        'g':'#4daf4a',
        'T':'#f781bf',
        't':'#f781bf',
        'U':'#f781bf',
        'u':'#f781bf',
        'N':'#ffffff',
        'n':'#ffffff',
        }

    nPlots = 1
    if resquiggleBorders is not None:
        nPlots += 1
    if eventalignBorders is not None:
        nPlots += 1

    fig, ax = plt.subplots(nrows = nPlots, sharex=True, figsize=(130,10), dpi=400)
    if nPlots == 1:
        ax = [ax]
    fig.suptitle(f'{readid} segmentation in 3\' -> 5\' orientation')
    x=np.arange(-polyAend, len(signal) - polyAend, 1)
    for axis in range(nPlots):
        if axis < nPlots-1:
            ax[axis].plot(x, signal, color='black', label='PolyA Standardized Signal', linewidth=0.8)
            ax[axis].set_ylim((0, 150))
        else:
            ax[axis].plot(x, normSignal, color='black', label='Z Normalised Signal', linewidth=0.8)
            # ax[axis].plot(normSignal, color='black', label='Z Normalised Signal', linewidth=0.8)
            ax[axis].set_ylim((-6, 3))
            if probability:
                ax2 = ax[axis].twinx()
                ax2.plot(x, probs, linewidth=1, label="log(Border Probability)", alpha=0.8)
                # ax2.plot(probs, linewidth=1, label="log(Border Probability)", alpha=0.8)
                ax2.set_ylim((-8, 28))
                ax2.set_xticks(np.arange(-8, 29, 4))
                ax2.set_ylabel('log(Border Probability)')
        ax[axis].set_ylabel('Signal pico Ampere')
        ax[axis].set_xticks(np.arange(-((polyAend//1000)*1000), len(signal), 1000))
        ax[axis].grid(True, 'both', 'y')

    plotNumber = 0
    # F5C EVENTALIGN
    if eventalignBorders is not None:
        # ax[plotNumber].vlines(eventalignBorders[:, 0].astype(int), ymin=20, ymax=150, colors='darkgreen', linestyles='--', label=f'f5c eventalign segmentation', linewidth=0.7)

        for i, segment in enumerate(eventalignBorders):
            motif = ''
            for j in range(i-2, i+3, 1):
                if j < 0:
                    motif += 'N'
                elif j >= len(eventalignBorders):
                    motif += 'N'
                else:
                    base = eventalignBorders[j][2]
                    motif += base
            # motif is stored in 3' - 5' direction
            motif = motif.replace('U', 'T')

            # ax[plotNumber].text(int(border[0]) + (int(border[1]) - int(border[0]))/2 - 6, 30, base, fontdict={'size' : 7, 'color':'black'})
            ax[plotNumber].text(int(segment[0]) + (int(segment[1]) - int(segment[0]))/2 - 6, 10, motif, fontdict={'size' : 6, 'color':'black'}, rotation=90)
            ax[plotNumber].vlines([int(segment[0])], ymin=0, ymax=150, colors=basecolors[read[segment[2]]], linestyles='--', linewidth=0.7)
            ax[plotNumber].add_patch(Rectangle((int(segment[0]), 0), int(segment[1]) - int(segment[0]), 150, alpha=0.4, facecolor=basecolors[read[segment[2]]]))

            x = int(segment[0])
            width = int(segment[1]) - int(segment[0])
            mean, stdev = STANDARDONTMODEL.loc[motif][['level_mean', 'level_stdv']]
            height = 3.92*stdev
            ax[plotNumber].add_patch(Rectangle((x, mean-1.96*stdev), width, height, alpha=0.4, facecolor="grey"))
            ax[plotNumber].hlines(y=mean, xmin=int(segment[0]), xmax=int(segment[1]), color='grey', linestyle='--', alpha=0.8)
            ax[plotNumber].hlines(y=mean+1.96*stdev, xmin=int(segment[0]), xmax=int(segment[1]), color='grey', linewidth=1, linestyle=':', alpha=0.8)
            ax[plotNumber].hlines(y=mean-1.96*stdev, xmin=int(segment[0]), xmax=int(segment[1]), color='grey', linewidth=1, linestyle=':', alpha=0.8)
        
        ax[plotNumber].vlines([int(segment[0])], ymin=0, ymax=150, colors=basecolors[read[segment[2]]], linestyles='--', linewidth=0.7, label = "f5c eventalign segmentation")
        ax[plotNumber].hlines(y=mean, xmin=int(segment[0]), xmax=int(segment[1]), color='grey', alpha=0.8, linestyle='--', label="ONT Model Mean")
        ax[plotNumber].hlines(y=mean+1.96*stdev, xmin=int(segment[0]), xmax=int(segment[1]), color='grey', linewidth=1, linestyle=':', alpha=0.8, label="95% conf. Interval")
        ax[plotNumber].legend()
        plotNumber += 1

    # F5C RESQUIGGLE
    if resquiggleBorders is not None:
        # ax[plotNumber].vlines(resquiggleBorders[:, 0].astype(int), ymin=20, ymax=150, colors='darkgreen', linestyles='--', label=f'f5c resquiggle segmentation', linewidth=0.7)

        for i, segment in enumerate(resquiggleBorders):
            motif = ''
            for j in range(i-2, i+3, 1):
                if j < 0:
                    motif += 'N'
                elif j >= len(resquiggleBorders):
                    motif += 'N'
                else:
                    base = resquiggleBorders[j][2]
                    motif += read[base]
            # motif is stored in 3' - 5' direction
            motif = motif.replace('U', 'T')

            # ax[plotNumber].text(int(border[0]) + (int(border[1]) - int(border[0]))/2 - 6, 30, base, fontdict={'size' : 7, 'color':'black'})
            ax[plotNumber].text(int(segment[0]) + (int(segment[1]) - int(segment[0]))/2 - 6, 10, motif[::-1], fontdict={'size' : 6, 'color':'black'}, rotation=90)
            ax[plotNumber].vlines([int(segment[0])], ymin=0, ymax=150, colors=basecolors[read[segment[2]]], linestyles='--', linewidth=0.7)
            ax[plotNumber].add_patch(Rectangle((int(segment[0]), 0), int(segment[1]) - int(segment[0]), 150, alpha=0.4, facecolor=basecolors[read[segment[2]]]))

            x = segment[0]
            width = segment[1] - segment[0]
            mean, stdev = STANDARDONTMODEL.loc[motif][['level_mean', 'level_stdv']]
            height = 3.92*stdev
            ax[plotNumber].add_patch(Rectangle((x, mean-1.96*stdev), width, height, alpha=0.4, facecolor="grey"))
            ax[plotNumber].hlines(y=mean, xmin=int(segment[0]), xmax=int(segment[1]), color='grey', linestyle='--', alpha=0.8)
            ax[plotNumber].hlines(y=mean+1.96*stdev, xmin=int(segment[0]), xmax=int(segment[1]), color='grey', linewidth=1, linestyle=':', alpha=0.8)
            ax[plotNumber].hlines(y=mean-1.96*stdev, xmin=int(segment[0]), xmax=int(segment[1]), color='grey', linewidth=1, linestyle=':', alpha=0.8)

        ax[plotNumber].vlines([int(segment[0])], ymin=0, ymax=150, colors=basecolors[read[segment[2]]], linestyles='--', linewidth=0.7, label = "f5c resquiggle segmentation")
        ax[plotNumber].hlines(y=mean, xmin=int(segment[0]), xmax=int(segment[1]), color='grey', alpha=0.8, linestyle='--', label="ONT Model Mean")
        ax[plotNumber].hlines(y=mean+1.96*stdev, xmin=int(segment[0]), xmax=int(segment[1]), color='grey', linewidth=1, linestyle=':', alpha=0.8, label="95% conf. Interval")
        ax[plotNumber].legend()
        plotNumber += 1

    # OUR SEGMENTATION
    for segment in segments:
        # extract motif from read
        pos = segment[2]
        motif = segment[6]
        if motif == 'NA':
            motif = segment[4]
            # add Ns to 5' end
            if pos-2 < 0:
                continue
                # motif = (' ' * abs(pos-2)) + motif
            # add As to 3' end
            elif pos+3 > len(read):
                continue
                # motif = motif + (' '*(pos+3-len(read)))
        base = motif[len(motif)//2]

        # motif in 5' - 3' direction
        motif = motif.replace('U', 'T')
        x = segment[0]
        width = segment[1] - segment[0]

        ax[plotNumber].vlines([int(segment[0])], ymin=-6, ymax=3, colors=basecolors[base], linestyles='--', linewidth=0.7)
        ax[plotNumber].add_patch(Rectangle((int(segment[0]), -6), width, 9, alpha=0.4, facecolor=basecolors[base]))

        match segment[5]:
            # for 3D: move in t and k, but not n
            # for 2D: move in t and n
            case 'M' | 'P':
                mean, stdev = kmermodels.loc[motif][['level_mean', 'level_stdv']]
                # print(mean, stdev)
                height = 3.92*stdev
                ax[plotNumber].add_patch(Rectangle((x, mean-1.96*stdev), width, height, alpha=0.4, facecolor="grey"))
                # write motif
                # ax[plotNumber].text(x + width/2 - 6, -3.5, read[pos], fontdict={'size' : 7, 'color':'black'})
                ax[plotNumber].text(x + width/2 - 6, -3.5, motif, fontdict={'size' : 6, 'color':'black'}, rotation=90)
                # draw kmer range as rectangle
                ax[plotNumber].hlines(y=mean, xmin=int(segment[0]), xmax=int(segment[1]), color='grey', linestyle='--', alpha=0.8)
                ax[plotNumber].hlines(y=mean+1.96*stdev, xmin=int(segment[0]), xmax=int(segment[1]), color='grey', linewidth=1, linestyle=':', alpha=0.8)
                ax[plotNumber].hlines(y=mean-1.96*stdev, xmin=int(segment[0]), xmax=int(segment[1]), color='grey', linewidth=1, linestyle=':', alpha=0.8)
            # for 2D: move in t, but not n and k
            case 'D':
                ax[plotNumber].add_patch(Rectangle((x, min(normSignal)), width, 9, alpha=0.3, facecolor="grey"))
                ax[plotNumber].text(x + width/2 - 3, -3, f'Del {motif}', rotation=90, fontdict={'size' : 6, 'color' : 'grey'})
            # for 3D: move in n, but not t and k
            case 'I':
                ax[plotNumber].text(x - 3, -3, f'Ins {motif}', rotation=90, fontdict={'size' : 6, 'color' : 'grey'})
            # move in t and n, but not k
            # case 'S':
            #     pass

    ax[plotNumber].vlines([int(segment[0])], ymin=20, ymax=150, colors=basecolors[base], linestyles='--', linewidth=0.7, label = "Dynamont Segmentation")
    ax[plotNumber].hlines(y=mean, xmin=int(segment[0]), xmax=int(segment[1]), color='grey', alpha=0.8, linestyle='--', label="Model Mean")
    ax[plotNumber].hlines(y=mean+1.96*stdev, xmin=int(segment[0]), xmax=int(segment[1]), color='grey', linewidth=1, linestyle=':', alpha=0.8, label="95% conf. Interval")
    ax[plotNumber].plot([0, 0.1], [-6.1, -6], c='blue', label='log(Border Probability)')
    ax[plotNumber].legend()

    # plt.legend()
    plt.savefig(join(outpath, readid + '.svg'), dpi=500)
    plt.savefig(join(outpath, readid + '.pdf'), dpi=500)

def segmentRead(signal : np.ndarray, normSignal : np.ndarray, ts : int, read : str, readid : str, outdir : str, resquiggleBorders : np.ndarray, eventalignBorders : np.ndarray, mode : str, modelpath : str, minSegLen : int, probability : bool, pore : str):
    '''
    Takes read in 3' -> 5' direction
    '''
    print(f"Segmenting {readid}")

    kmermodels = pd.read_csv(modelpath, sep='\t', index_col = "kmer")

    if name == 'nt': # check for windows
        CPP_SCRIPT+='.exe'
    if mode == 'basic' or mode == 'banded':
        mode = 'dynamont_NT' if mode == 'basic' else 'dynamont_NT_banded'
        PARAMS = {
            'e1': 1.0,
            'm1': 0.03,
            'e2': 0.97
            }
    elif mode == 'resquiggle':
        mode = 'dynamont_NTK'
        PARAMS = {
            'a1': 0.015537006200000003,
            'a2': 0.34805735,
            'p1': 0.008918512069999999,
            'p2': 0.06678974900000001,
            'p3': 0.047764121937800004,
            's1': 0.0031768724926999996,
            's2': 0.013912044499999998,
            's3': 0.49361836,
            'e1': 1.0,
            'e2': 0.9968231100000001,
            'e3': 0.9910815099999999,
            'e4': 0.9020646600000001,
            'i1': 0.0016965344,
            'i2': 0.11056013505290001
            }
    else:
        print(f'Mode {mode} not implemented')
        exit(1)
    
    CPP_SCRIPT = join(dirname(__file__), f'{mode}')

    PARAMS['m'] = modelpath
    PARAMS['c'] = minSegLen
    PARAMS['p'] = probability
    PARAMS['r'] = pore

    segments, probs = feedSegmentation(normSignal[ts:], read, CPP_SCRIPT, PARAMS)

    # check for resquiggle how many new segments were inserted
    print('Read bases: ', len(read), 'Segments: ', len(segments))

    if not len(segments):
        raise SegmentationError(readid)

    segments[:, 0] = segments[:, 0] - ts
    segments[:, 1] = segments[:, 1] - ts

    with open(join(outdir, readid + '.txt'), 'w') as w:
        w.write('\n'.join(list(map(str, segments))))

    plotBorders(signal, normSignal, ts, read[::-1], segments, probs, readid, outdir, resquiggleBorders, eventalignBorders, kmermodels, probability)

def start(dataPath : str, basecalls : str, targetID : str, outdir : str, resquigglePickle : str, eventalignPickle : str, mode : str, modelpath : str, minSegLen : int, probability : bool, pore : str) -> tuple:

    with pysam.AlignmentFile(basecalls, "r" if basecalls.endswith('.sam') else "rb", check_sq=False) as samfile:
        for basecalled_read in samfile.fetch(until_eof=True):

            # init read
            readid = basecalled_read.query_name
            if readid != targetID:
                continue
            seq = basecalled_read.query_sequence
            qual = basecalled_read.get_tag("qs")
            ts = basecalled_read.get_tag("ts")
            rawFile = join(dataPath, basecalled_read.get_tag("fn"))
            r5 = read5.read(rawFile)
            print("Avg read quality: ", qual, "Transcript start in signal: ", ts)
            
            resquiggleBorders = None
            if resquigglePickle:
                print("Prepare f5c resquiggle")
                with open(resquigglePickle, 'rb') as handle:
                    borderMap = pickle.load(handle)
                    if targetID in borderMap:
                        resquiggleBorders = borderMap[targetID]
                        try:
                            resquiggleBorders[:, 0] = resquiggleBorders[:, 0].astype(int) - ts
                            resquiggleBorders[:, 1] = resquiggleBorders[:, 1].astype(int) - ts
                        except KeyError:
                            print(targetID, "not in polyA")
                            # print("ERROR resquiggle", resquiggleBorders)
                            exit(1)
                        # print(borders)
                    else:
                        resquiggleBorders = None
                        print(f'WARNING: no border found in {handle}')

            eventalignBorders = None
            if eventalignPickle:
                print("Prepare f5c eventalign")
                with open(eventalignPickle, 'rb') as handle:
                    borderMap = pickle.load(handle)
                    if targetID in borderMap:
                        eventalignBorders = borderMap[targetID]
                        try:
                            eventalignBorders[:, 0] = eventalignBorders[:, 0].astype(int) - ts
                            eventalignBorders[:, 1] = eventalignBorders[:, 1].astype(int) - ts
                        except KeyError:
                            print(targetID, "not in polyA")
                            # print("ERROR eventalign", eventalignBorders)
                            exit(1)
                        # print(borders)
                    else:
                        eventalignBorders = None
                        print(f'WARNING: no border found in {handle}')
            
            signal = r5.getPolyAStandardizedSignal(targetID, ts-105, ts-5)
            normSignal = r5.getZNormSignal(targetID, "mean")
            normSignal = hampel(normSignal, 6, 5.).filtered_data # small window and high variance allowed: just to filter outliers that result from sensor errors, rest of the original signal should be kept

            # change read from 5'-3' to 3'-5'
            segmentRead(signal, normSignal, ts, seq[::-1], targetID, outdir, resquiggleBorders, eventalignBorders, mode, modelpath, minSegLen, probability, pore)

def main() -> None:
    args = parse()
    if not exists(args.outdir):
        makedirs(args.outdir)
    start(args.raw, args.basecalls, args.readid, args.outdir, args.resquigglePickle, args.eventalignPickle, args.mode, args.model_path, args.minSegLen - 1, args.probability, args.pore)

if __name__ == '__main__':
    main()