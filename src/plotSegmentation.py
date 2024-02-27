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
import re
from os.path import exists, join, dirname
from os import makedirs, name
from fileio import getFiles, loadFastx, openCPPScriptParams, stopFeeding, calcZ
from read5 import read
from subprocess import Popen
import h5py as h5

KMERMODELS = pd.read_csv('/home/yi98suv/projects/ont_segmentation/data/template_median69pA_extended.model', sep='\t', index_col = "kmer")

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
PARAMS = {'e1': 1.0, 'm2': 0.062403286666666655, 'd1': 0.0, 'e2': 0.9373050666666667, 'e3': 0.0002915558866666667, 'i1': 0.0, 'm3': 0, 'i2': 0, 'm4': 0, 'd2': 0}

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
    parser.add_argument('--seg_pickle', type=str, default=None, help='f5c segmentation pickle file')
    return parser.parse_args()

def readPolyA(file : str) -> pd.DataFrame:
    df = pd.read_csv(file, usecols=['readname', 'transcript_start'], sep='\t')
    df = df.astype({'readname' : str, 'transcript_start' : int})
    # df.set_index('readname', inplace=True)
    return pd.Series(df.transcript_start.values, index=df.readname).to_dict()

def plotBorders(signal : np.ndarray, read : str, segments : np.ndarray, readid : str, outpath : str, borders : np.ndarray = None):
    '''
    Input
    -----
    segments : np.ndarray
        in 3' -> 5' orientation
    read : str
        in 3' -> 5' orientation
    '''
    fig, ax = plt.subplots(figsize=(120,8), dpi=500)
    fig.suptitle(f'{readid} segmentation in 3\' -> 5\' orientation')
    ax.plot(signal, color='blue', label='original signal', linewidth=1)
    ax.set_ylim((0, max(signal)+10))
    ax.set_ylabel('Signal pico Ampere')
    ax.set_xticks(np.arange(0, len(signal), 2000))
    if borders is not None:
        ax.vlines(borders, ymin=0, ymax=max(signal)+10, colors='green', linestyles='-', label='true segmentation', linewidth=1, alpha=0.6)
        # for j in range(len(borders) - 1):
        #     plt.text(borders[j] + (borders[j+1]-borders[j])/2 - 3, 45, read[j], fontdict={'size' : 6, 'color':'green'})
    ax.vlines(segments[:, 0], ymin=0, ymax=max(signal)+10, colors='red', linestyles='--', label='segmentation', linewidth=1, alpha=0.6)
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
            mean, stdev = KMERMODELS.loc[motif][['level_mean', 'sd_mean']]
            height = 6*stdev
            ax.add_patch(Rectangle((x, mean-3*stdev), width, height, alpha=0.3, facecolor="yellow", edgecolor='red'))
            # write motif
            plt.text(x + width/2 - 3, 40, read[pos], fontdict={'size' : 6, 'color':'red'})
        elif segment[3] == 'D':
            height = max(signal) - min(signal)
            ax.add_patch(Rectangle((x, min(signal)), width, height, alpha=0.3, facecolor="grey", edgecolor='black'))
            plt.text(x + width/2 - 3, 35, 'Del', rotation=90, fontdict={'size' : 6, 'color' : 'grey'})
        elif segment[3] == 'I':
            plt.text(x - 3, 35, 'Ins', rotation=90, fontdict={'size' : 6, 'color' : 'grey'})

    plt.legend()
    plt.grid(True, 'both', 'both')
    plt.savefig(join(outpath, readid + '.svg'))
    plt.savefig(join(outpath, readid + '.pdf'))

# https://stackoverflow.com/questions/32570029/input-to-c-executable-python-subprocess
def feedSegmentation(signal : np.ndarray, read : str, pipe : Popen) -> np.ndarray:
    '''
    Parse & feed signal & read to the C++ segmentation script.

    Parameters
    ----------
    signal : np.ndarray
    read : str
    stream
        Open stdin stream of the C++ segmentation algorithm

    Returns
    -------
    segmentation : np.ndarray
        in 3'-5' orientation
    '''
    # prepare cookie for segmentation
    cookie = f"{str(signal.tolist()).replace(' ', '').replace('[', '').replace(']', '')} {read}\n"
    # transfer data to bytes - needed in Python 3
    cookie = bytes(cookie, 'UTF-8')
    # feed cookie to segmentation
    pipe.stdin.write(cookie)
    pipe.stdin.flush()
    output = pipe.stdout.readline().strip().decode('UTF-8')
    output = re.findall('\D[\d,]+', output)
    # print(output)
    probs = pipe.stdout.readline().strip()[:-1].decode('UTF-8')
    print("Output length", len(output), len(probs))
    try:
        probs = np.array(probs.split(','), dtype=float)[:-1] # TODO border means this position is included in the segment? then exclude last element
    except:
        print(pipe.stderr.readlines())
        # print(cookie)
        exit(1)
    segments = []
    for i in range(len(output)):
        state = output[i][0]
        basepos, start = map(int, output[i][1:].split(','))
        try:
            _, end = map(int, output[i+1][1:].split(','))
        except IndexError:
            end = len(signal)
        segments.append([start, end, basepos, state])
    segments = np.array(segments, dtype=object)
    assert segments.size, f"{signal}\n==\n{read}"
    return segments, probs

def segmentRead(signal : np.ndarray, read : str, readid : str, outdir : str, borders : np.ndarray = None):
    '''
    Takes read in 3' -> 5' direction
    '''

    # Create pipe between python script and cp segmentation script
    CPP_SCRIPT = join(dirname(__file__), 'segmentation')
    if name == 'nt': # check for windows
        CPP_SCRIPT+='.exe'
    pipe = openCPPScriptParams(CPP_SCRIPT, PARAMS)

    segments, probs = feedSegmentation(signal, read, pipe)
    stopFeeding(pipe)

    with open(join(outdir, readid + '.txt'), 'w') as w:
        w.write('\n'.join(list(map(str, segments))))

    # print(len(segments), segments)
    # print(len(read), read)

    plotBorders(signal, read, segments, readid, outdir, borders)
    print(calcZ(signal, read, PARAMS, CPP_SCRIPT))

def start(files, basecalls, targetID, polyA, out, seg_pickle) -> tuple:
    for file in files:
        r5 = read(file)
        if targetID in r5.getReads():
            if seg_pickle:
                import pickle
                with open(seg_pickle, 'rb') as handle:
                    borderMap = pickle.load(handle)
                    borders = np.array(borderMap[targetID]) - polyA[targetID]
            else:
                try:
                    # borders from simulated data
                    borders = h5.File(file)[f'read_{targetID}/Raw/Borders']
                except:
                    # no borders
                    borders = None
                # change read from 5'-3' to 3'-5'
            segmentRead(r5.getpASignal(targetID)[polyA[targetID]:], basecalls[targetID][::-1], targetID, out, borders)

def main() -> None:
    args = parse()
    if not exists(args.out):
        makedirs(args.out)
    polya = readPolyA(args.polya)
    files = getFiles(args.raw, True)
    print(f'ONT Files: {len(files)}')
    basecalls = loadFastx(args.fastx)
    # print("5' -> 3'", len(basecalls[args.readid]), basecalls[args.readid].replace("U", "T"))
    print(f'Segmenting {len(basecalls)} reads')
    start(files, basecalls, args.readid, polya, args.out, args.seg_pickle)

if __name__ == '__main__':
    main()