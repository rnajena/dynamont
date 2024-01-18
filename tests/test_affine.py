#!/usr/bin/env python
# author: Jannes Spangenberg
# e-mail: jannes.spangenberg@uni-jena.de
# github: https://github.com/JannesSP
# website: https://jannessp.github.io

import gzip
from argparse import ArgumentDefaultsHelpFormatter, ArgumentParser, Namespace
from os import makedirs, name
from os.path import dirname, exists, join
from pathlib import Path
from subprocess import PIPE, Popen
from sys import maxsize
from itertools import repeat
import multiprocessing as mp
# from scipy.signal import argrelextrema
import re

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from Bio import SeqIO
# from hampel import hampel
from read5 import read

TERM_STRING = "$"

def openCPPScript(cpp_script : str) -> Popen:
    '''
    Open cpp script with Popen.

    Parameters
    ----------
    cpp_script : str
        Path of script

    Returns
    -------
    subprocess : Popen
    '''
    return Popen([cpp_script], shell=True, stdout=PIPE, stdin=PIPE, stderr=PIPE)

# https://stackoverflow.com/questions/32570029/input-to-c-executable-python-subprocess
def feedSegmentation(signal : list, read : str, pipe : Popen) -> np.ndarray:
    '''
    Parse & feed signal & read to the C++ segmentation script.

    Parameters
    ----------
    signal : list
    read : str
    stream
        Open stdin stream of the C++ segmentation algorithm

    Returns
    -------
    segmentation : np.ndarray
    '''
    # prepare cookie for segmentation
    cookie = f"{str(signal).replace(' ', '').replace('[', '').replace(']', '')} {read}\n"
    # transfer data to bytes - needed in Python 3
    cookie = bytes(cookie, 'UTF-8')
    # feed cookie to segmentation
    pipe.stdin.write(cookie)
    # print(cookie)
    pipe.stdin.flush()
    output = pipe.stdout.readline().strip().decode('UTF-8')
    output = re.findall('\D[\d,]+', output)
    # print(output)
    probs = pipe.stdout.readline().strip()[:-1].decode('UTF-8')
    print(len(output), len(probs))
    try:
        probs = np.array(probs.split(','), dtype=float)[:-1] # TODO border means this position is included in the segment? then exclude last element
    except:
        print(pipe.stderr.readlines())
        print(cookie)
        exit(1)
    segments = []
    for i in range(len(output)):
        state = output[i][0]
        basepos, start = map(int, output[i][1:].split(','))
        try:
            _, end = map(int, output[i+1][1:].split(','))
        except IndexError:
            end = len(signal)
        segments.append([start, end, len(read) - basepos, state])
    segments = np.array(segments, dtype=object)
    # print("Segments: ", segments[:10], len(segments))
    return segments, probs

def stopFeeding(pipe : Popen) -> None:
    pipe.stdin.write(bytes(f'{TERM_STRING} {TERM_STRING}\n', 'UTF-8'))
    pipe.stdin.close()
    pipe.stdout.close()

def getFiles(filepath : str, rec : bool) -> list:
    '''
    Returns
    -------
    files : list
        a list of input files
    '''

    if filepath.endswith(".fast5") or filepath.endswith(".slow5") or filepath.endswith(".blow5") or filepath.endswith(".pod5"):
        return [filepath]

    func = {True : Path(filepath).rglob, False : Path(filepath).glob}[rec]
    files = []

    for path in func('*.fast5'):
        files.append(str(path))

    for path in func('*.slow5'):
        files.append(str(path))

    for path in func('*.blow5'):
        files.append(str(path))

    for path in func('*.pod5'):
        files.append(str(path))

    return files

def loadFastqs(path : str) -> dict:
    '''
    Returns
    -------
    readDict : dict
        {readid : read}
    '''
    readDict = {}
    if path.endswith(".gz"):
        path = gzip.open(path, "rt")
    for record in SeqIO.parse(path, 'fastq'):
        readDict[record.id] = str(record.seq)
    return readDict

def readNanoPolyA(file : str) -> pd.DataFrame:
    return pd.read_csv(file, sep='\t')

def segmentRead(signal : np.ndarray, read : str, readid : str, out : str):
    '''
    Takes read in 3' -> 5' direction
    '''

    # Create pipe between python script and cp segmentation script
    # CPP_SCRIPT = join(dirname(__file__), 'segment')
    CPP_SCRIPT = join(dirname(__file__), '..', 'src', 'segment_affine')
    # CPP_SCRIPT = join(dirname(__file__), '..', 'src', 'segment_affine_deletion')
    if name == 'nt': # check for windows
        CPP_SCRIPT+='.exe'
    pipe = openCPPScript(CPP_SCRIPT)

    np.set_printoptions(formatter={'float': lambda x: "{0:0.1f}".format(x)}, threshold=maxsize)
    # filter signal with hampel filter
    # signal = hampel(signal, window_size=30).filtered_data
    segments, probs = feedSegmentation(signal, read, pipe)
    # borders += transcript_start
    print(segments)
    print(probs)

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
    res = []
    with open(file, 'r') as seg:
        seg.readline() # header
        for line in seg:
            line = line.strip().split()
            res.append([int(line[-2]), line[9][2]])
    return res

def main() -> None:
    out = dirname(__file__)
    outfile = join(out, 'result.csv')
    outsum = join(out, 'summary.csv')
    o = open(outsum, 'w')
    o.write('readindex,readid\n')
    signal=[108.90141,108.90141,108.90141,108.90141,108.90141,108.90141,108.90141,108.90141,108.90141,108.90141,107.754234,107.754234,107.754234,107.754234,107.754234,107.754234,107.754234,107.754234,107.754234,107.754234,87.18801,87.18801,87.18801,87.18801,87.18801,87.18801,87.18801,87.18801,87.18801,87.18801,73.256,73.256,73.256,73.256,73.256,73.256,73.256,73.256,73.256,73.256,73.256,73.256,73.256,101.34711,101.34711,101.34711,101.34711,101.34711,101.34711,101.34711,101.34711,101.34711,101.34711,101.34711,101.34711,101.34711,101.34711,101.34711]
    read="TGCCAAAAA"
    with open(outfile, 'w') as w:
        w.write('readindex,position,base,segmentstart,segmentend\n')

        segmentRead(signal, read, 0, out)
    
    o.close()

if __name__ == '__main__':
    main()