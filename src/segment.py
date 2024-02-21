#!/usr/bin/env python
# author: Jannes Spangenberg
# e-mail: jannes.spangenberg@uni-jena.de
# github: https://github.com/JannesSP
# website: https://jannessp.github.io

from argparse import ArgumentDefaultsHelpFormatter, ArgumentParser, Namespace
from os import makedirs, name
from os.path import dirname, exists, join
from sys import maxsize
from itertools import repeat
import multiprocessing as mp

import matplotlib.pyplot as plt
import numpy as np
from read5 import read
from fileio import feedSegmentation, stopFeeding, openCPPScript, getFiles, loadFastx

def parse() -> Namespace:
    parser = ArgumentParser(
        formatter_class=ArgumentDefaultsHelpFormatter
    )
    parser.add_argument('ONTpath', type=str, help='Path to one file or folder containing ONT raw data, either {fast5, slow5, blow5 or pod5}')
    parser.add_argument('FastQ', type=str, help='Multi fastq files containing basecalls for provided ONT data')
    # parser.add_argument('nanopolish_polyA', type=str, help='nanopolish polya output table')
    parser.add_argument('outfile_basename', type=str, help='Output directory')
    parser.add_argument('-r', '--recursive', action='store_true', help='Look for ONT data recursively in folder')
    parser.add_argument('-p', '--processes', type=int, default=1, help='Number of processes for multiprocessing')
    parser.add_argument('--rna', action='store_true', help='Input signal is RNA data')
    
    # TODO parameters for testing, might remove later
    # parser.add_argument('-ri', '--readid', default=None, help='Plot only the segmentation for a single read with the provided readid.')
    # parser.add_argument('segmentation', default=None, help='Segmentation from nanopolish/f5c with readnames as readids')
    return parser.parse_args()

def segmentRead(signal : np.ndarray, read : str, readid : str, out : str):
    '''
    Takes read in 3' -> 5' direction
    '''

    # TODO for later mutliprocessing move the pipe out of this function
    # Create pipe between python script and cp segmentation script
    CPP_SCRIPT = join(dirname(__file__), 'segment_affine_deletion')
    if name == 'nt': # check for windows
        CPP_SCRIPT+='.exe'
    pipe = openCPPScript(CPP_SCRIPT)

    np.set_printoptions(formatter={'float': lambda x: "{0:0.1f}".format(x)}, threshold=maxsize)
    # filter signal with hampel filter
    # signal = hampel(signal, window_size=30).filtered_data
    segments, probs = feedSegmentation(signal, read, pipe)
    # borders += transcript_start

    stopFeeding(pipe)

    res = np.stack((list(repeat(readid, len(segments))), segments[:, 2], list(read) + ['S'], segments[:,3], segments[:,0], segments[:,1]), axis = 1)
    # print("Result: ", res[:10], len(res))
    res = res[np.apply_along_axis(lambda row: int(row[1]), arr=res, axis=1).argsort()]
    
    with open(out + '.txt', 'a+') as w:
        w.write(formatOut(res))

        plotSegmentationProbs(signal, probs, res, out + '.svg')
    # plotHampelSegmentationProbs(hampel_signal, signal, probs, borders[:-1], read, outpath)

    # TODO for testing, might remove later
    # plotSegmentationProbs(hampel_signal, signal, hampel_probs, hampel_borders, original_probs, original_borders, read, segmentation, outpath)
    # plotBorders(hampel_signal, signal, hampel_borders, borders, read, segmentation, outpath)

def formatOut(res : np.ndarray) -> str:
    return '\n'.join(','.join(map(str, e)) for e in res) + '\n'

def plotSegmentationProbs(signal : np.ndarray, probs : np.ndarray, segments : np.ndarray, outpath : str):
    
    matches = segments[segments[:, 3] == 'M']

    fig, ax1 = plt.subplots(figsize=(120,8), dpi=300)
    fig.suptitle('Signal segmentation in 3\' -> 5\' orientation')
    ax1.plot(signal, color='blue', label='original signal', linewidth=0.5)
    ax1.set_ylim((0, max(signal)+10))
    for i in range(len(matches)):
        ax1.annotate(matches[i, 2], (matches[i, 4], -10), fontsize=7, annotation_clip=False)
    ax1.set_ylabel('Signal pico Ampere')
    ax1.set_xticks(np.arange(0, len(signal), 2000))

    ax2 = ax1.twinx()
    ymin = -20
    ymax = 40
    ax2.plot(np.arange(0, len(probs)), probs, label='border probabilities', color='grey', alpha=0.6, linewidth=1)
    ax2.fill_between(np.arange(0, len(probs)), ymin, probs, color='grey', alpha=0.4)
    ax2.vlines(matches[:,-2], ymin=ymin, ymax=ymax, colors='red', linestyles='--', label='our segmentation', linewidth=1, alpha=0.6)
    ax2.set_ylim((ymin, ymax))
    ax2.grid(True, 'major', 'y', linestyle=':', alpha=0.6)
    ax2.set_xlabel('Position in the ONT Signal')
    ax2.set_ylabel('Summed Transition-Probability per Signal-Datapoint')
    h1, l1 = ax1.get_legend_handles_labels()
    h2, l2 = ax2.get_legend_handles_labels()
    plt.legend([*h1, *h2], [*l1, *l2])
    plt.tight_layout()
    plt.savefig(outpath)

def plotHampelSegmentationProbs(hampel_signal : np.ndarray, signal : np.ndarray, probs : np.ndarray, borders : np.ndarray, read : str, outpath : str):

    fig, ax1 = plt.subplots(figsize=(120,8), dpi=300)
    fig.suptitle('Signal segmentation in 3\' -> 5\' orientation')
    ax1.plot(signal, color='blue', label='original signal', linewidth=2)
    ax1.plot(hampel_signal, color='red', label='hampel signal', linewidth=0.5)
    ax1.set_ylim((0, max(signal)+10))
    for i, e in enumerate(borders):
        ax1.annotate(read[i], (e, -10), fontsize=7, annotation_clip=False)
    ax1.set_ylabel('Signal pico Ampere')
    ax1.set_xticks(np.arange(0, len(signal), 2000))

    ax2 = ax1.twinx()
    ymin = -20
    ymax = 40
    ax2.plot(np.arange(0, len(probs)), probs, label='hampel probabilities', color='grey', alpha=0.6, linewidth=1)
    ax2.fill_between(np.arange(0, len(probs)), ymin, probs, color='grey', alpha=0.4)
    ax2.vlines(borders, ymin=ymin, ymax=ymax, colors='red', linestyles='--', label='our segmentation', linewidth=1, alpha=0.6)
    ax2.set_ylim((ymin, ymax))
    ax2.grid(True, 'major', 'y', linestyle=':', alpha=0.6)
    ax2.set_xlabel('Position in the ONT Signal')
    ax2.set_ylabel('Summed Transition-Probability per Signal-Datapoint')
    h1, l1 = ax1.get_legend_handles_labels()
    h2, l2 = ax2.get_legend_handles_labels()
    plt.legend([*h1, *h2], [*l1, *l2])
    plt.tight_layout()
    plt.savefig(join(outpath, 'segmentation_probs.svg'))

# def plotSegmentationProbs(hampel_signal : np.ndarray, signal : np.ndarray, hampel_probs : np.ndarray, hampel_borders : np.ndarray, original_probs : np.ndarray, original_borders : np.ndarray, read : str, f5c_segs : list, outpath : str):

#     fig, ax1 = plt.subplots(figsize=(120,8), dpi=300)
#     fig.suptitle('Signal segmentation in 3\' -> 5\' orientation')
#     ax1.plot(signal, color='blue', label='original signal', linewidth=2)
#     ax1.plot(hampel_signal, color='red', label='hampel signal', linewidth=0.5)
#     ax1.set_ylim((0, max(signal)+10))
#     ax1.vlines([e[0] for e in f5c_segs], ymin=0, ymax=max(signal)+10, colors='grey', linestyles='--', label='f5c segmentation', linewidth=1, alpha=0.6)
#     for e in f5c_segs:
#         ax1.annotate(e[0], (e[0] + 1, max(signal) + 3), rotation=90, fontsize=7)
#         ax1.annotate(e[1], (e[0] + 1, max(signal)), fontsize=7)
#     for i, e in enumerate(hampel_borders):
#         ax1.annotate(read[i], (e, -10), fontsize=7, annotation_clip=False)
#     ax1.set_ylabel('Signal pico Ampere')
#     ax1.set_xticks(np.arange(0, len(signal), 2000))
#     # ax1.grid(True, 'both', 'x')

#     ax2 = ax1.twinx()
#     ymin = -20
#     ymax = 40
#     ax2.plot(np.arange(0, len(original_probs)), original_probs, label='original probabilities', color='grey', linewidth=1)
#     ax2.fill_between(np.arange(0, len(original_probs)), ymin, original_probs, color='grey', alpha=0.6)
#     ax2.plot(np.arange(0, len(hampel_probs)), hampel_probs, label='hampel probabilities', color='red', alpha=0.6, linewidth=1)
#     ax2.vlines(original_borders, ymin=ymin, ymax=ymax, colors='blue', linestyles='--', label='our segmentation', linewidth=1, alpha=0.6)
#     ax2.vlines(hampel_borders, ymin=ymin, ymax=ymax, colors='red', linestyles='--', label='our segmentation', linewidth=1, alpha=0.6)
#     ax2.set_ylim((ymin, ymax))
#     # ax2.set_xlim((0, len(signal)))
#     ax2.grid(True, 'major', 'y', linestyle=':', alpha=0.6)
#     ax2.set_xlabel('Position in the ONT Signal')
#     ax2.set_ylabel('Summed Transition-Probability per Signal-Datapoint')
#     # ax2.set_xticks(np.arange(0, len(segmentprobs), l))
#     h1, l1 = ax1.get_legend_handles_labels()
#     h2, l2 = ax2.get_legend_handles_labels()
#     plt.legend([*h1, *h2], [*l1, *l2])
#     plt.tight_layout()
#     plt.savefig(join(outpath, 'segmentation_probs.svg'))

def plotBorders(hampel_signal : np.ndarray, signal : np.ndarray, hampel_borders : np.ndarray, borders : np.ndarray, read : str, segmentation : list, outpath : str):

    fig, ax1 = plt.subplots(figsize=(120,8), dpi=300)
    fig.suptitle('Signal segmentation in 3\' -> 5\' orientation')
    ax1.plot(signal, color='blue', label='original signal', linewidth=2)
    ax1.plot(hampel_signal, color='red', label='hampel_signal', linewidth=0.5)
    ax1.set_ylim((0, max(signal)+10))
    ax1.vlines([e[0] for e in segmentation], ymin=0, ymax=max(signal)+10, colors='grey', linestyles='--', label='f5c segmentation', linewidth=1, alpha=0.6)
    for e in segmentation:
        ax1.annotate(e[0], (e[0] + 1, max(signal) + 3), rotation=90, fontsize=7)
        ax1.annotate(e[1], (e[0] + 1, max(signal)), fontsize=7)
    ax1.set_ylabel('Signal pico Ampere')
    ax1.set_xticks(np.arange(0, len(signal), 2000))
    ax1.vlines(borders, ymin=0, ymax=max(signal)+10, colors='blue', linestyles='--', label='segmentation of original', linewidth=1, alpha=0.6)
    ax1.vlines(hampel_borders, ymin=0, ymax=max(signal)+10, colors='red', linestyles='--', label='segmentation of hampel', linewidth=1, alpha=0.6)
    h1, l1 = ax1.get_legend_handles_labels()
    plt.legend(h1, l1)
    plt.tight_layout()
    plt.savefig(join(outpath, 'segmentation_probs.svg'))

def loadSegmentation(file : str) -> list:
    '''
    Read segmentation from nanopolish/f5c eventalign csv
    '''
    res = []
    with open(file, 'r') as seg:
        seg.readline() # header
        for line in seg:
            line = line.strip().split()
            res.append([int(line[-2]), line[9][2]])
    return res

def main() -> None:
    args = parse()
    inp = args.ONTpath
    out = args.outfile_basename
    fastq = args.FastQ
    isRna = args.rna
    # nanoPolyA = readNanoPolyA(args.nanopolish_polyA).set_index('readname')
    # seg = loadSegmentation(args.segmentation)

    assert exists(inp), f'Input path {inp} does not exist!'
    assert exists(fastq), f'Fastq path {fastq} does not exist!'

    if not exists(dirname(out)):
        makedirs(dirname(out))
    # TODO change
    outfile = out + '_result.csv'
    outsum = out + '_summary.csv'
    o = open(outsum, 'w')
    o.write('readindex,readid\n')
    with open(outfile, 'w') as w:
        w.write('readindex,position,base,segmentstart,segmentend\n')

    files = getFiles(inp, args.recursive)
    print(f'ONT Files: {len(files)}')
    basecalls = loadFastx(fastq)
    print(f'Segmenting {len(basecalls)} reads')

    # TODO
    # pool = mp.Pool(args.processes)
    # polyAIndex = nanoPolyA.index.to_list()
    i = -1
    for file in files:
        r5 = read(file)
        for readid in r5:
            if (readid not in basecalls): # or (readid not in polyAIndex)
                continue
            # transcript_start = int(nanoPolyA.loc[[readid]]['transcript_start'][0].item())
            # if transcript_start == -1:
            #     continue
            i+=1
            o.write(f'{i},{readid}\n')
            print(readid)
            # invert read from 5'->3' to 3'->5' to match the signal orientation
            # pool.apply_async(segmentRead, (r5.getpASignal(readid), basecalls[readid][::-1], i, transcript_start, outfile))

            segmentRead(r5.getpASignal(readid), basecalls[readid][::-1], i, out)
            exit(1)

            # pool.apply_async(segmentRead, (r5.getpASignal(readid), basecalls[readid][::-1], i, out))

            # change basecall orientation to 3' -> 5' because models are in 3' -> 5' direction
            # segmentRead(r5.getpASignal(readid), basecalls[readid][::-1], i, transcript_start, outfile)
    
    o.close()
    # pool.close()
    # pool.join()

if __name__ == '__main__':
    main()