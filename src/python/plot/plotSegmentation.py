#!/usr/bin/env python
# author: Jannes Spangenberg
# e-mail: jannes.spangenberg@uni-jena.de
# github: https://github.com/JannesSP
# website: https://jannessp.github.io

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import read5_ont
import pysam
import seaborn as sns
from argparse import ArgumentDefaultsHelpFormatter, ArgumentParser, Namespace
from matplotlib.patches import Rectangle
from os.path import exists, join
from os import makedirs, name
from src.python.segmentation.FileIO import feedSegmentation, SegmentationError, hampelFilter

def parse() -> Namespace:
    parser = ArgumentParser(
        formatter_class=ArgumentDefaultsHelpFormatter
    )
    parser.add_argument('-r', '--raw',   type=str, required=True, metavar="POD5|FAST5", help='Raw ONT training data')
    parser.add_argument('-b', '--basecalls', type=str, required=True, metavar="BAM", help='Basecalls of ONT training data as .bam file')
    parser.add_argument('-o', '--outdir',   type=str, required=True, metavar="PATH", help='Outpath to write files')
    parser.add_argument('--readid', type=str, required=True, help='Read to plot')
    parser.add_argument('--pore',  type=str, required=True, choices=["rna_r9", "dna_r9", "rna_rp4", "dna_r10_260bps", "dna_r10_400bps"], help='Pore generation used to sequence the data')
    # optional
    parser.add_argument('--f5cResquiggle', type=str, default=None, help='f5c resquiggle segmentation file')
    parser.add_argument('--model_path', type=str, required=True, help='Kmer model file')
    parser.add_argument('--mode', required=True, choices=['basic', 'resquiggle'])
    parser.add_argument('--probability', action="store_true", help="Output the segment border probability per position.")
    parser.add_argument("--changepoints", type=str, default=None, metavar="HDF5", help="HDF5 file with ground truth change points")
    return parser.parse_args()

def plotBorders(normSignal : np.ndarray, start : int, end : int, read : str, segments : np.ndarray, probs : np.ndarray, readid : str, outpath : str, kmerModels : pd.DataFrame, prob : bool, resquiggleBorders : np.ndarray, pore : str, changepoints : np.ndarray):
    '''
    Input
    -----
    segments : np.ndarray
        in 5' -> 3' orientation
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
    
    sns.set_theme()

    motifLength = 5 if "r9" in pore else 9

    lb, ub = -5, 3
    nPlots = 1
    if resquiggleBorders is not None:
        nPlots += 1

    fig, ax = plt.subplots(nrows = nPlots, figsize=(120,10), dpi=300)
    fig.suptitle(f'{readid} segmentation in 3\' -> 5\' orientation')
    if nPlots == 1:
        ax = [ax]
    
    for axis in range(nPlots):
        ax[axis].plot(normSignal, color='black', label='Normalized Signal', linewidth=0.8)
        ax[axis].set_ylim((lb, ub))
        ax[axis].set_ylabel('Normalized pA Signal')
        ax[axis].set_xticks(np.arange(0, len(normSignal), 1000))
        # ax[axis].grid(True, 'both', 'y')

    plotNumber = 0
    # F5C RESQUIGGLE
    if resquiggleBorders is not None:
        for segment in resquiggleBorders:
            if segment[2] + 5 > len(read):
                continue
            motif = ''
            for j in range(segment[2], segment[2] + motifLength):
                if j < 0:
                    motif += 'N'
                elif j >= len(read):
                    motif += 'N'
                else:
                    motif += read[j]
            motif = motif.replace('U', 'T')
            base = motif[len(motif)//2]

            ax[plotNumber].text(int(segment[0]) + (int(segment[1]) - int(segment[0]))/2 - 6, -3.5, motif, fontdict={'size' : 6, 'color':'black'}, rotation=90)
            ax[plotNumber].vlines([int(segment[0])], ymin=lb, ymax=ub, colors=basecolors[base], linestyles='--', linewidth=0.7)
            ax[plotNumber].add_patch(Rectangle((int(segment[0]), lb), int(segment[1]) - int(segment[0]), ub-lb, alpha=0.4, edgecolor=basecolors[base], facecolor=basecolors[base]))

        ax[plotNumber].vlines([int(segment[0])], ymin=lb, ymax=ub, colors=basecolors[base], linestyles='--', linewidth=0.7, label = "f5c resquiggle segmentation")
        ax[plotNumber].legend('upper right')
        
        if changepoints is not None:
            ax[plotNumber].set_xticks(changepoints, minor=True)
            ax[plotNumber].tick_params(axis='x', which='minor', direction='out', length=5, color='red', labelbottom=False, bottom=True, top=True)

        plotNumber += 1

    # OUR SEGMENTATION
    for segment in segments:
        # print(segment)
        # break
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

        ax[plotNumber].vlines([int(segment[0])], ymin=lb, ymax=ub, colors=basecolors[base], linestyles='--', linewidth=0.7)
        ax[plotNumber].add_patch(Rectangle((int(segment[0]), lb), width, ub-lb, alpha=0.4, facecolor=basecolors[base], edgecolor=basecolors[base]))

        match segment[5]:
            # for 3D: move in t and k, but not n
            # for 2D: move in t and n
            case 'M' | 'P':
                mean, stdev = kmerModels.loc[motif][['level_mean', 'level_stdv']]
                height = 3.92*stdev # 1.96 * 2 * stdev
                ax[plotNumber].add_patch(Rectangle((segment[0], mean-1.96*stdev), width, height, alpha=0.3, facecolor="grey", edgecolor="grey"))
                # write motif
                ax[plotNumber].text(segment[0] + width/2 - 6, -3.5, motif, fontdict={'size' : 6, 'color':'black'}, rotation=90)
                # draw kmer range as rectangle
                ax[plotNumber].hlines(y=mean, xmin=int(segment[0]), xmax=int(segment[1]), color='grey', linestyle='--', alpha=0.6)
                ax[plotNumber].hlines(y=mean+1.96*stdev, xmin=int(segment[0]), xmax=int(segment[1]), color='grey', linewidth=1, linestyle=':', alpha=0.6)
                ax[plotNumber].hlines(y=mean-1.96*stdev, xmin=int(segment[0]), xmax=int(segment[1]), color='grey', linewidth=1, linestyle=':', alpha=0.6)
            case 'I':
                ax[plotNumber].text(segment[0] - 3, -3, f'Ins {motif}', rotation=90, fontdict={'size' : 6, 'color' : 'grey'})
            # move in t and n, but not k
            # case 'S':
            #     pass

    ax[plotNumber].vlines([int(segment[0])], ymin=lb, ymax=ub, colors=basecolors[base], linestyles='--', linewidth=0.7, label = "Dynamont Segmentation")
    ax[plotNumber].hlines(y=mean, xmin=int(segment[0]), xmax=int(segment[1]), color='grey', alpha=0.8, linestyle='--', label="Model Mean")
    ax[plotNumber].hlines(y=mean+1.96*stdev, xmin=int(segment[0]), xmax=int(segment[1]), color='grey', linewidth=1, linestyle=':', alpha=0.8, label="95% conf. Interval")
    ax[plotNumber].plot([0, 0.1], [-6.1, -6], c='black', label='log(Border Probability)')
    ax[plotNumber].legend(loc='upper right')
    
    if changepoints is not None:
        ax[plotNumber].set_xticks(changepoints, minor=True)
        ax[plotNumber].tick_params(axis='x', which='minor', direction='out', length=5, color='red', labelbottom=False, bottom=True, top=True)

    if prob:
        twinax = ax[plotNumber].twinx()
        twinax.plot(np.arange(start, start+len(probs)), probs, linewidth=0.7, color='black', label="log(Border Probability)")
        twinax.set_ylim((-6, 42))
        twinax.set_yticks(np.arange(-6, 43, 6))
        twinax.set_ylabel('log(Border Probability)')

    plt.hlines([start, end], lb, ub, colors='black')
    plt.grid(False)

    plt.savefig(join(outpath, readid + '.svg'), dpi=500)
    plt.savefig(join(outpath, readid + '.pdf'), dpi=500)

def segmentRead(normSignal : np.ndarray, start : int, end : int, read : str, readid : str, outdir : str, mode : str, modelPath : str, probability : bool, pore : str, resquiggleBorders : np.ndarray, changepoints : np.ndarray):
    '''
    Takes read in 3' -> 5' direction
    '''
    print(f"Segmenting {readid}")

    kmerModels = pd.read_csv(modelPath, sep='\t', index_col = "kmer")
    PARAMS = {}

    if mode == 'basic':
        CPP_SCRIPT = 'dynamont-NT-banded'
    elif mode == 'resquiggle':
        CPP_SCRIPT = 'dynamont-NTC'
    else:
        print(f'Mode {mode} not implemented')
        exit(1)
        
    if name == 'nt': # check for windows
        CPP_SCRIPT+='.exe'

    PARAMS['m'] = modelPath
    PARAMS['p'] = probability
    PARAMS['r'] = pore

    # if "r9" in pore:
    kmerSize = 5
    # else:
        # kmerSize = 9

    segments, borderProbs = feedSegmentation(normSignal[start:end], read, CPP_SCRIPT, start, kmerSize, PARAMS) # , heatmap

    # sns.set_theme()
    # plt.figure(dpi=200)
    # im = plt.imshow(heatmap, aspect='auto', cmap='viridis')
    # cbar = plt.colorbar(im)  # Add a color bar
    # cbar.set_label('Probability')  # Set the color bar label
    # # sns.heatmap(heatmap, cmap='viridis', annot=False, cbar_kws={'label': 'Probability'}, linewidths=0)
    # plt.title('HMM Probability Matrix with MAP Path')
    # plt.ylabel("Signal Data Points")
    # plt.xlabel("Nucleotide Position")
    # plt.tight_layout()
    # plt.grid(False)
    # plt.savefig("heatmap.svg", dpi=300, bbox_inches='tight', pad_inches=0)
    # plt.savefig("heatmap.pdf", dpi=300, bbox_inches='tight', pad_inches=0)
    # plt.close()

    # check for resquiggle how many new segments were inserted
    print('Read bases:', len(read), ' Segments:', len(segments), ' Border probs:', len(borderProbs), ' Signal:', end-start)

    if not len(segments):
        raise SegmentationError(readid)

    with open(join(outdir, readid + '.txt'), 'w') as w:
        w.write('\n'.join(list(map(str, segments))))

    if "rna" in pore: # change orientation back from 3'-5' to 5'-3'
        read = read[::-1]

    plotBorders(normSignal, start, end, read, segments, borderProbs, readid, outdir, kmerModels, probability, resquiggleBorders, pore, changepoints)

def start(dataPath : str, basecalls : str, targetID : str, outdir : str, mode : str, modelpath : str, probability : bool, pore : str, f5cReadMap : dict, changepoints : np.ndarray) -> tuple:

    resquiggleBorders = None
    if f5cReadMap and targetID in f5cReadMap:
        resquiggleBorders = f5cReadMap[targetID]

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
            r5 = read5_ont.read(rawFile)

            # normSignal = r5.getZNormSignal(readid, "median")[sp:sp+ns].astype(np.float32)
            shift = basecalledRead.get_tag("sm")
            scale = basecalledRead.get_tag("sd")
            if pore in ["dna_r9", "rna_r9"]:
                # for r9 pores, shift and scale are stored for pA signal in bam
                signal = r5.getpASignal(readid)
            else:
                # for new pores, shift and scale is directly applied to stored integer signal (DACs)
                # this way the conversion from DACs to pA is skipped
                signal = r5.getSignal(readid)
            signal = (signal - shift) / scale
            hampelFilter(signal, 6, 5.) # small window and high variance allowed: just to filter outliers that result from sensor errors, rest of the original signal should be kept

            # change read from 5'-3' to 3'-5'
            if "rna" in pore:
                seq = seq[::-1]
                if not seq.startswith("AAAAAAAAA"):
                    seq = "AAAAAAAAA" + seq

            segmentRead(signal, sp+ts, sp+ns, seq, basecalledRead.query_name, outdir, mode, modelpath, probability, pore, resquiggleBorders, changepoints)

def readF5CResquiggle(file: str) -> dict:
    """
    Parses a TSV file to extract read IDs along with their corresponding start and end positions.

    Parameters
    ----------
    file : str
        Path to the TSV file containing resquiggle segmentation output.

    Returns
    -------
    dict
        A dictionary mapping each read ID to a set of start and end positions.
    """
    print("Reading f5c resquiggle output from " + file)
    readMap = {}
    with open(file, 'r') as f:
        next(f) # skip header
        # contig    position    reference_kmer  read_index  strand  event_index event_level_mean    event_stdv  event_length    model_kmer  model_mean  model_stdv  standardized_level  start_idx   end_idx
        for line in f:
            readid, kmer_idx, start, end = line.strip().split('\t')
            if start == '.' or end == '.':
                continue
            if readid not in readMap:
                readMap[readid] = []
            readMap[readid].append([int(start), int(end), int(kmer_idx)])
    return readMap

def readChangepoints(file : str, targetID : str) -> np.ndarray:
    import h5py
    print("Reading changepoints from " + file)
    with h5py.File(file, 'r') as h5:
        for readid in h5.keys():
            if readid != targetID:
                continue

            changepoints = h5[readid + "/waveletEdge"][:]
            res = changepoints

    return res

def main() -> None:
    args = parse()
    if not exists(args.outdir):
        makedirs(args.outdir)

    f5cReadMap = readF5CResquiggle(args.f5cResquiggle) if args.f5cResquiggle is not None else None
    changepoints = readChangepoints(args.changepoints, args.readid) if args.changepoints is not None else None

    start(args.raw, args.basecalls, args.readid, args.outdir, args.mode, args.model_path, args.probability, args.pore, f5cReadMap, changepoints)

if __name__ == '__main__':
    main()