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
from FileIO import feedSegmentation, SegmentationError, hampelFilter
import read5
import pysam

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
    parser.add_argument('--model_path', type=str, required=True, help='Kmer model file')
    parser.add_argument('--mode', required=True, choices=['basic', 'banded', 'resquiggle'])
    parser.add_argument('--probability', action="store_true", help="Output the segment border probability per position.")
    
    return parser.parse_args()

def plotBorders(normSignal : np.ndarray, ts : int, read : str, segments : np.ndarray, probs : np.ndarray, readid : str, outpath : str, kmerModels : pd.DataFrame, prob : bool):
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
    baesColors = {
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

    plt.figure(figsize=(130, 10))
    plt.title(f'{readid} segmentation in 3\' -> 5\' orientation')
    plt.plot(normSignal, color='black', label='Z Normalised Signal', linewidth=0.8)
    plt.ylim((-6, 3))
    plt.ylabel('Signal pico Ampere')
    plt.xticks(np.arange(0, len(normSignal), 1000))
    plt.grid(True, 'both', 'y')

    # OUR SEGMENTATION
    for segment in segments:
        # start, end, pos, read[pos], read[max(0, pos-2):min(len(read), pos+3)], state, prob, polish
        # extract motif from read
        pos = segment[2]
        # NTK version
        motif = segment[7]
        # for NT version
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
        width = segment[1] - segment[0]

        plt.vlines([int(segment[0])], ymin=-6, ymax=3, colors=baesColors[base], linestyles='--', linewidth=0.7)
        ax = plt.gca()
        ax.add_patch(Rectangle((int(segment[0]), -6), width, 9, alpha=0.4, facecolor=baesColors[base]))

        match segment[5]:
            # for 3D: move in t and k, but not n
            # for 2D: move in t and n
            case 'M' | 'P':
                mean, stdev = kmerModels.loc[motif][['level_mean', 'level_stdv']]
                height = 3.92*stdev # 1.96 * 2 * stdev
                ax.add_patch(Rectangle((segment[0], mean-1.96*stdev), width, height, alpha=0.4, facecolor="grey"))
                # write motif
                plt.text(segment[0] + width/2 - 6, -3.5, motif, fontdict={'size' : 6, 'color':'black'}, rotation=90)
                # draw kmer range as rectangle
                plt.hlines(y=mean, xmin=int(segment[0]), xmax=int(segment[1]), color='grey', linestyle='--', alpha=0.8)
                plt.hlines(y=mean+1.96*stdev, xmin=int(segment[0]), xmax=int(segment[1]), color='grey', linewidth=1, linestyle=':', alpha=0.8)
                plt.hlines(y=mean-1.96*stdev, xmin=int(segment[0]), xmax=int(segment[1]), color='grey', linewidth=1, linestyle=':', alpha=0.8)
            case 'I':
                plt.text(segment[0] - 3, -3, f'Ins {motif}', rotation=90, fontdict={'size' : 6, 'color' : 'grey'})
            # move in t and n, but not k
            # case 'S':
            #     pass

    plt.vlines([int(segment[0])], ymin=20, ymax=150, colors=baesColors[base], linestyles='--', linewidth=0.7, label = "Dynamont Segmentation")
    plt.hlines(y=mean, xmin=int(segment[0]), xmax=int(segment[1]), color='grey', alpha=0.8, linestyle='--', label="Model Mean")
    plt.hlines(y=mean+1.96*stdev, xmin=int(segment[0]), xmax=int(segment[1]), color='grey', linewidth=1, linestyle=':', alpha=0.8, label="95% conf. Interval")
    plt.plot([0, 0.1], [-6.1, -6], c='blue', label='log(Border Probability)')
    plt.legend()

    if prob:
        plt.twinx()
        # print(x, polyAend, probs)
        plt.plot(np.arange(ts, ts+len(probs)), probs, linewidth=1, label="log(Border Probability)", alpha=0.8)
        plt.ylim((-8, 28))
        plt.yticks(np.arange(-8, 29, 4))
        plt.ylabel('log(Border Probability)')

    # plt.legend()
    plt.savefig(join(outpath, readid + '.svg'), dpi=500)
    plt.savefig(join(outpath, readid + '.pdf'), dpi=500)

def segmentRead(normSignal : np.ndarray, ts : int, read : str, readid : str, outdir : str, mode : str, modelPath : str, probability : bool, pore : str):
    '''
    Takes read in 3' -> 5' direction
    '''
    print(f"Segmenting {readid}")

    kmerModels = pd.read_csv(modelPath, sep='\t', index_col = "kmer")

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
            'a1': 0.012252440188168037,
            'a2': 0.246584724985145,
            'p1': 0.04477093133243305,
            'p2': 0.007687811003133089,
            'p3': 0.4469623669791557,
            's1': 0.05321209670114726,
            's2': 0.0007555035568187239,
            's3': 0.21999557711272136,
            'e1': 1.0,
            'e2': 0.9467879033992115,
            'e3': 0.9552290685034269,
            'e4': 0.9792321612614708,
            'i1': 7.208408117990252e-05,
            'i2': 0.08645733058947891
            }
    else:
        print(f'Mode {mode} not implemented')
        exit(1)
    
    CPP_SCRIPT = join(dirname(__file__), f'{mode}')

    PARAMS['m'] = modelPath
    PARAMS['p'] = probability
    PARAMS['r'] = pore

    segments, borderProbs = feedSegmentation(normSignal[ts:], read, CPP_SCRIPT, ts, PARAMS)

    # check for resquiggle how many new segments were inserted
    print('Read bases:', len(read), ' Segments:', len(segments), ' Border probs:', len(borderProbs), ' Signal:', len(normSignal[ts:]))

    if not len(segments):
        raise SegmentationError(readid)

    with open(join(outdir, readid + '.txt'), 'w') as w:
        w.write('\n'.join(list(map(str, segments))))

    if "rna" in pore: # change orientation back from 3'-5' to 5'-3'
        read = read[::-1]

    plotBorders(normSignal, ts, read, segments, borderProbs, readid, outdir, kmerModels, probability)

def start(dataPath : str, basecalls : str, targetID : str, outdir : str, mode : str, modelpath : str, probability : bool, pore : str) -> tuple:

    with pysam.AlignmentFile(basecalls, "r" if basecalls.endswith('.sam') else "rb", check_sq=False) as samfile:
        for basecalledRead in samfile.fetch(until_eof=True):

            # init read, sometimes a read got split by the basecaller and got a new id
            readid = basecalledRead.get_tag("pi") if basecalledRead.has_tag("pi") else basecalledRead.query_name
            
            # print(readid)
            if readid != targetID:
                continue
            seq = basecalledRead.query_sequence
            # qs = basecalled_read.get_tag("qs") # quality score
            sp = basecalledRead.get_tag("sp") if basecalledRead.has_tag("sp") else 0 # if split read get start offset of the signal
            ts = basecalledRead.get_tag("ts") # ts:i: 	the number of samples trimmed from the start of the signal
            ns = basecalledRead.get_tag("ns") # ns:i: 	the number of samples in the signal prior to trimming

            rawFile = join(dataPath, basecalledRead.get_tag("fn"))
            r5 = read5.read(rawFile)

            # print(basecalled_read.get_tag('sm'), basecalled_read.get_tag('sd'))
            # print('median complete', r5.getShift(readid, 'median'), r5.getScale(readid, 'median'))
            # print('median sp:sp+ns', np.median(r5.getpASignal(readid)[sp:sp+ns]), mad(r5.getpASignal(readid)[sp:sp+ns]))
            # print('median sp+ts:sp+ns', np.median(r5.getpASignal(readid)[sp+ts:sp+ns]), mad(r5.getpASignal(readid)[sp+ts:sp+ns]))
            # print('mean', r5.getShift(readid, 'mean'), r5.getScale(readid, 'mean'))
            # print('mean sp:sp+ns', np.mean(r5.getpASignal(readid)[sp:sp+ns]), np.std(r5.getpASignal(readid)[sp:sp+ns]))
            # print('mean sp+ts:sp+ns', np.mean(r5.getpASignal(readid)[sp+ts:sp+ns]), np.std(r5.getpASignal(readid)[sp+ts:sp+ns]))

            # normSignal = r5.getZNormSignal(readid, "median")[sp:sp+ns].astype(np.float32)
            shift = basecalledRead.get_tag("sm")
            scale = basecalledRead.get_tag("sd")
            signal = r5.getpASignal(readid)[sp+ts:sp+ns].astype(np.float32)
            normSignal = (signal - shift) / scale
            normSignal = hampelFilter(normSignal, 6, 5.) # small window and high variance allowed: just to filter outliers that result from sensor errors, rest of the original signal should be kept

            # change read from 5'-3' to 3'-5'
            if "rna" in pore:
                seq = seq[::-1]
                
            segmentRead(normSignal, ts, seq, basecalledRead.query_name, outdir, mode, modelpath, probability, pore)

def main() -> None:
    args = parse()
    if not exists(args.outdir):
        makedirs(args.outdir)
    start(args.raw, args.basecalls, args.readid, args.outdir, args.mode, args.model_path, args.probability, args.pore)

if __name__ == '__main__':
    main()