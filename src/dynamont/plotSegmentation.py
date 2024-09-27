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
from dynamont.FileIO import getFiles, loadFastx, readPolyAStartEnd, feedSegmentation, SegmentationError
from read5 import read
from hampel import hampel
import pickle
# from matplotlib import ticker

STANDARDONTMODEL = pd.read_csv("/home/yi98suv/projects/dynamont/data/template_median69pA_extended.model", sep='\t', index_col = "kmer")

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
    parser.add_argument('--mode', choices=['basic', 'basic_sparsed', 'indel', '3d', '3d_sparsed', '3d_reduced'], required=True)
    parser.add_argument('--model', type=str, default=join(dirname(__file__), '..', '..', 'data', 'norm_models', 'rna_r9.4_180mv_70bps_extended_stdev0_5.model'), help='Kmer model file')
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
        in 3' -> 5' orientation
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
                ax2.set_ylim((-4, 14))
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
    # ax[plotNumber].vlines(segments[:, 0], ymin=-4, ymax=4, colors='red', linestyles='--', label='our segmentation', linewidth=0.6)
    # plot kmer model range
    # print(read)
    for segment in segments:
        # print(segment)
        # exit(1)
        # print(segment)
        # extract motif from read
        pos = segment[2]
        motif = segment[6]
        if motif == 'NA':
            motif = segment[4]
            # add Ns to 5' end
            if pos-2 < 0:
                motif = ('N' * abs(pos-2)) + motif
            # add As to 3' end
            elif pos+3 > len(read):
                motif = motif + ('N'*(pos+3-len(read)))
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

    # plt.figure(figsize=(10, 10), dpi=500)
    # plt.xlim((19000, 21000))
    # plt.savefig(join(outpath, readid + '_poster.svg'), dpi=500)
    # plt.savefig(join(outpath, readid + '_poster.pdf'), dpi=500)

    # plt.figure(figsize=(12,8))
    # plt.gcf().set_size_inches(12, 8)
    # plt.xlim(18000, 18400)
    # plt.ylim(45, 110)
    # plt.gca().xaxis.set_major_locator(ticker.MultipleLocator(100))
    # plt.savefig(join(outpath, readid + '_ausschnitt.svg'), dpi=300)


def segmentRead(signal : np.ndarray, normSignal : np.ndarray, polyAstart : int, polyAend : int, read : str, readid : str, outdir : str, resquiggleBorders : np.ndarray, eventalignBorders : np.ndarray, mode : str, modelpath : str, minSegLen : int, probability : bool):
    '''
    Takes read in 3' -> 5' direction
    '''
    print(f"Segmenting {readid}")

    kmermodels = pd.read_csv(modelpath, sep='\t', index_col = "kmer")

    if name == 'nt': # check for windows
        CPP_SCRIPT+='.exe'
    if mode == 'indel':
        PARAMS = {
            'e1': 1.0,
            'm1': 0.047923630044895006,
            'd1': 0.0795918281369677,
            'e2': 0.8719994721033932,
            'i1': 0.00048506999577316453,
            'm2': 0.006990024477599084,
            'i2': 0.019397852195642322,
            'm3': 0.9806021470164152,
            'd2': 0.9930099746422482
            }
    elif mode == 'basic':
        mode = 'basic_sparsed'
        PARAMS = {
            'e1': 1.0,
            'm1': 0.03189915859979101,
            'e2': 0.9681008434763126
            }
    elif mode == '3d':
        mode = '3d_sparsed'
        # a1, a2, p1, p2, p3, s1, s2, s3, e1, e2, e3, e4, i1, i2
        # PARAMS = {
        #     'a1': 0.029304646362105145,
        #     'a2': 0.4589929809984935,
        #     'p1': 0.44631725138060524,
        #     'p2': 0.22329377700463882,
        #     'p3': 0.14091038339921583,
        #     's1': 0.010868603047271925,
        #     's2': 0.0019430768236690966,
        #     's3': 0.24016500608414348,
        #     'e1': 1.0,
        #     'e2': 0.9891315629003756,
        #     'e3': 0.5536827484371548,
        #     'e4': 0.7449262248177602,
        #     'i1': 0.0005322041406781533,
        #     'i2': 0.15993149851340224
        #     }
        PARAMS = {'a1': 0.030570760232128354, 'a2': 0.6101330423858251, 'p1': 0.0383208584222504, 'p2': 0.08685413259505585, 'p3': 0.015545502978626788, 's1': 0.0010061921517067177, 's2': 0.001534865889565432, 's3': 0.13149698483529376, 'e1': 1.0, 'e2': 0.9989938127816652, 'e3': 0.9616791335752949, 'e4': 0.8801771554834845, 'i1': 0.0008630902783022252, 'i2': 0.24282445342753348}
    elif mode == '3d_reduced':
        mode = '3d_sparsed_reduced'
        # a1, a2, p1, p2, p3, s1, s2, s3, e1, e2, e3, e4, i1, i2
        PARAMS = {'a1': 0.030570760232128354, 'a2': 0.6101330423858251, 'p1': 0.0383208584222504, 'p2': 0.08685413259505585, 'p3': 0.015545502978626788, 's1': 0.0010061921517067177, 's2': 0.001534865889565432, 's3': 0.13149698483529376, 'e1': 1.0, 'e2': 0.9989938127816652, 'e3': 0.9616791335752949, 'e4': 0.8801771554834845, 'i1': 0.0008630902783022252, 'i2': 0.24282445342753348}
        # PARAMS = {
        #     'a1': 0.030582340361124723,
        #     'a2': 0.5519926024794314,
        #     'p1': 0.0760246390600083,
        #     'p2': 0.05755050242472108,
        #     'p3': 0.013846319526016653,
        #     's1': 0.0011123625979078346,
        #     's2': 0.0005556412414201246,
        #     's3': 0.0901126314806047,
        #     'e1': 1.0,
        #     'e2': 0.9988876515465803,
        #     'e3': 0.9239753253210568,
        #     'e4': 0.9105210852714761,
        #     'i1': 0.0007904494265846322,
        #     'i2': 0.3440484496408686
        # }
    
    CPP_SCRIPT = join(dirname(__file__), '..', 'dynamont', f'segmentation_{mode}')

    PARAMS['m'] = modelpath
    PARAMS['c'] = minSegLen
    PARAMS['p'] = probability

    # filter outliers
    # hampel_std_signal = hampel(standardizedSignal, 20, 2.).filtered_data
    # hampel_raw_signal = hampel(rawSignal, 20, 2.).filtered_data
    # segments = feedSegmentation(hampel_std_signal[polyAstart:], read, CPP_SCRIPT, PARAMS)

    # feedSignal = normSignal[polyAend:]
    segments, probs = feedSegmentation(normSignal, read, CPP_SCRIPT, PARAMS)

    if not len(segments):
        # print(segments)
        # print(str(list(normSignal[-300:])).replace(" ", "")[1:-1], end=' ')
        # print(read[-30:])
        raise SegmentationError(readid)

    segments[:, 0] = segments[:, 0] - polyAend
    segments[:, 1] = segments[:, 1] - polyAend

    with open(join(outdir, readid + '.txt'), 'w') as w:
        w.write('\n'.join(list(map(str, segments))))

    plotBorders(signal, normSignal, polyAend, read[::-1], segments, probs, readid, outdir, resquiggleBorders, eventalignBorders, kmermodels, probability)
    # print(calcZ(normSignal, read, PARAMS, CPP_SCRIPT))

def start(files, basecalls, targetID, polyA, out, resquigglePickle, eventalignPickle, mode, modelpath, minSegLen, probability) -> tuple:
    for file in files:
        r5 = read(file)
        if targetID in r5.getReads():
            
            resquiggleBorders = None
            if resquigglePickle:
                print("Prepare f5c resquiggle")
                with open(resquigglePickle, 'rb') as handle:
                    borderMap = pickle.load(handle)
                    if targetID in borderMap:
                        resquiggleBorders = borderMap[targetID]
                        try:
                            resquiggleBorders[:, 0] = resquiggleBorders[:, 0].astype(int) - polyA[targetID][1]
                            resquiggleBorders[:, 1] = resquiggleBorders[:, 1].astype(int) - polyA[targetID][1]
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
                            eventalignBorders[:, 0] = eventalignBorders[:, 0].astype(int) - polyA[targetID][1]
                            eventalignBorders[:, 1] = eventalignBorders[:, 1].astype(int) - polyA[targetID][1]
                        except KeyError:
                            print(targetID, "not in polyA")
                            # print("ERROR eventalign", eventalignBorders)
                            exit(1)
                        # print(borders)
                    else:
                        eventalignBorders = None
                        print(f'WARNING: no border found in {handle}')
            
            # fill batch
            polyAstart = polyA[targetID][0] if targetID in polyA else 0
            polyAend = polyA[targetID][1] if targetID in polyA else 0
            if polyAend - polyAstart > 30:
                signal = r5.getPolyAStandardizedSignal(targetID, polyA[targetID][0]+5, polyA[targetID][1]-5)
            else:
                signal = r5.getpASignal(targetID)
            # if targetID in polyA and polyA[targetID][1] - polyA[targetID][0] > 30:
            normSignal = r5.getZNormSignal(targetID, "mean") #[polyA[targetID][1]:]
            normSignal = hampel(normSignal, 6, 5.).filtered_data # small window and high variance allowed: just to filter outliers that result from sensor errors, rest of the original signal should be kept

            # change read from 5'-3' to 3'-5'
            segmentRead(signal, normSignal, polyAstart, polyAend, basecalls[targetID][::-1], targetID, out, resquiggleBorders, eventalignBorders, mode, modelpath, minSegLen, probability)

def main() -> None:
    args = parse()
    # print(args)
    if not exists(args.out):
        makedirs(args.out)
    polya = readPolyAStartEnd(args.polya)
    rawFiles = getFiles(args.raw, True)
    print(f'ONT Files: {len(rawFiles)}')
    basecalls = loadFastx(args.fastx)
    # print("5' -> 3'", len(basecalls[args.readid]), basecalls[args.readid].replace("U", "T"))
    # print(f'Segmenting {len(basecalls)} reads')
    start(rawFiles, basecalls, args.readid, polya, args.out, args.resquigglePickle, args.eventalignPickle, args.mode, args.model, args.minSegLen - 1, args.probability)

if __name__ == '__main__':
    main()