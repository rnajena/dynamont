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
from scipy.signal import argrelextrema

import matplotlib.pyplot as plt
import numpy as np
from Bio import SeqIO
from hampel import hampel
from read5 import read

TERM_STRING = "$"

def parse() -> Namespace:
    parser = ArgumentParser(
        formatter_class=ArgumentDefaultsHelpFormatter
    )
    parser.add_argument('ONTpath', type=str, help='File to one file or folder containing ONT raw data, either {fast5, slow5, blow5 or pod5}')
    parser.add_argument('FastQ', type=str, help='Multi fastq files containing basecalls for provided ONT data')
    parser.add_argument('outpath', type=str, help='Output directory')
    parser.add_argument('-r', '--recursive', action='store_true', help='Look for ONT data recursively in folder')
    parser.add_argument('-p', '--processes', type=int, default=1, help='Number of processes for multiprocessing')
    parser.add_argument('--rna', action='store_true', help='Input signal is RNA data')
    
    # TODO parameters for testing, might remove later
    parser.add_argument('-ri', '--readid', default=None, help='Plot only the segmentation for a single read with the provided readid.')
    parser.add_argument('segmentation', default=None, help='Segmentation from nanopolish/f5c with readnames as readids')    
    return parser.parse_args()

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
    return Popen([cpp_script], shell=True, stdout=PIPE, stdin=PIPE)

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
    '''
    # prepare cookie for segmentation
    signal = str(signal.tolist()).replace(' ', '').replace('[', '').replace(']', '')
    cookie = f'{signal} {read}\n'
    # transfer data to bytes - needed in Python 3
    cookie = bytes(cookie, 'UTF-8')
    # feed cookie to segmentation
    pipe.stdin.write(cookie)
    pipe.stdin.flush()
    # prepare segmentation and return
    borders = pipe.stdout.readline().strip()[:-1].decode('UTF-8')
    probs = pipe.stdout.readline().strip()[:-1].decode('UTF-8')
    borders = np.array(borders.split(','), dtype=float)
    probs = np.array(probs.split(','), dtype=float)[:-1] # TODO border means this position is included in the segment? then exclude last element
    return borders, probs

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

def segmentRead(signal : np.ndarray, read : str, segmentation : list, outfile : str, outpath : str):
    '''
    Takes read in 5' -> 3' direction
    '''

    # TODO cut polyA
    signal = signal[4663:]
    segmentation = [(s[0] - 4663, s[1]) for s in segmentation]

    print('Readlength', len(read))
    # TODO for later mutliprocessing move the pipe out of this function
    # Create pipe between python script and cp segmentation script
    CPP_SCRIPT = join(dirname(__file__), 'segment')
    if name == 'nt': # check for windows
        CPP_SCRIPT+='.exe'
    pipe = openCPPScript(CPP_SCRIPT)

    np.set_printoptions(formatter={'float': lambda x: "{0:0.1f}".format(x)}, threshold=maxsize)
    # print(np.array_repr(signal).replace('\n', '').replace(' ', ''), read)
    hampel_signal = hampel(signal, window_size=30).filtered_data
    # print(hampel_signal)
    hampel_borders, hampel_probs = feedSegmentation(hampel_signal, read, pipe)
    original_borders, original_probs = feedSegmentation(signal, read, pipe)
    # print(max(hampel_probs))
    # print(len(signal), len(hampel_probs), len(original_probs))
    stopFeeding(pipe)

    with open(outfile, 'a+') as w:
        w.write(f"Segmentpobabilies (log):\n{np.array_str(original_borders)}\nRead:\n{read}\n")

    plotSegmentationProbs(hampel_signal, signal, hampel_probs, hampel_borders, original_probs, original_borders, read, segmentation, outpath)
    # plotBorders(hampel_signal, signal, hampel_borders, borders, read, segmentation, outpath)

def plotSegmentationProbs(hampel_signal : np.ndarray, signal : np.ndarray, hampel_probs : np.ndarray, hampel_borders : np.ndarray, original_probs : np.ndarray, original_borders : np.ndarray, read : str, f5c_segs : list, outpath : str):

    fig, ax1 = plt.subplots(figsize=(120,8), dpi=300)
    fig.suptitle('Signal segmentation in 3\' -> 5\' orientation')
    ax1.plot(signal, color='blue', label='original signal', linewidth=2)
    ax1.plot(hampel_signal, color='red', label='hampel signal', linewidth=0.5)
    ax1.set_ylim((0, max(signal)+10))
    ax1.vlines([e[0] for e in f5c_segs], ymin=0, ymax=max(signal)+10, colors='grey', linestyles='--', label='f5c segmentation', linewidth=1, alpha=0.6)
    for e in f5c_segs:
        ax1.annotate(e[0], (e[0] + 1, max(signal) + 3), rotation=90, fontsize=7)
        ax1.annotate(e[1], (e[0] + 1, max(signal)), fontsize=7)
    ax1.set_ylabel('Signal pico Ampere')
    ax1.set_xticks(np.arange(0, len(signal), 2000))
    # ax1.grid(True, 'both', 'x')

    ax2 = ax1.twinx()
    ymin = -20
    ymax = 40
    ax2.plot(np.arange(0, len(original_probs)), original_probs, label='original probabilities', color='grey', linewidth=1)
    ax2.fill_between(np.arange(0, len(original_probs)), ymin, original_probs, color='grey', alpha=0.6)
    ax2.plot(np.arange(0, len(hampel_probs)), hampel_probs, label='hampel probabilities', color='red', alpha=0.6, linewidth=1)
    ax2.vlines(original_borders, ymin=ymin, ymax=ymax, colors='blue', linestyles='--', label='our segmentation', linewidth=1, alpha=0.6)
    ax2.vlines(hampel_borders, ymin=ymin, ymax=ymax, colors='red', linestyles='--', label='our segmentation', linewidth=1, alpha=0.6)
    ax2.set_ylim((ymin, ymax))
    # ax2.set_xlim((0, len(signal)))
    ax2.grid(True, 'major', 'y', linestyle=':', alpha=0.6)
    ax2.set_xlabel('Position in the ONT Signal')
    ax2.set_ylabel('Summed Transition-Probability per Signal-Datapoint')
    # ax2.set_xticks(np.arange(0, len(segmentprobs), l))
    h1, l1 = ax1.get_legend_handles_labels()
    h2, l2 = ax2.get_legend_handles_labels()
    plt.legend([*h1, *h2], [*l1, *l2])
    plt.tight_layout()
    plt.savefig(join(outpath, 'segmentation_probs.svg'))

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
    args = parse()
    inp = args.ONTpath
    out = args.outpath
    fastq = args.FastQ
    isRna = args.rna

    seg = loadSegmentation(args.segmentation)

    assert exists(inp), f'Input path {inp} does not exist!'
    assert exists(fastq), f'Output path {fastq} does not exist!'

    if not exists(out):
        makedirs(out)
    # TODO change
    outfile = join(out, 'segmentation_result.txt')

    files = getFiles(inp, args.recursive)
    print(f'ONT Files: {len(files)}')
    basecalls = loadFastqs(fastq)
    print(f'Segmenting {len(basecalls)} reads')

    # TODO
    # processes = args.processes
    for file in files:
        r5 = read(file)
        if args.readid is not None:
            if args.readid not in r5.getReads():
                continue
            # change basecall orientation to 3' -> 5'
            segmentRead(r5.getpASignal(args.readid), basecalls[args.readid][::-1], seg, outfile, out)
        else:
            for readid in r5:
                segmentRead(r5.getpASignal(readid), basecalls[readid][::-1], seg, outfile, out)

    # #########################
    # Test
    # ret = feedSegmentation(np.array([107.2,108.0,108.9,111.2,105.7,104.3,107.1,105.7]), 'CAAAAA', pipe)
    # print(ret)
    # l = 5

    # signal = np.array(
    #     [96.86647]*l + # AAGCT
    #     [84.4253]*l + # AGCTA
    #     [95.15087]*l + # GCTAG
    #     [117.71723]*l + # CTAGC
    #     [87.507645]*l + # TAGCA
    #     [79.56085]*l + # AGCAT
    #     [80.2712]*l + # GCATT
    #     [82.19202]*l + # CATTG
    #     [118.77031]*l + # ATTGA
    #     [118.66758]*l + # TTGAT
    #     [99.587555]*l + # TGATC
    #     [83.94988]*l + # GATCC
    #     [68.43068]*l + # ATCCG
    #     [111.2103]*l + # TCCGA
    #     [103.08109]*l + # CCGAG
    #     [128.6229]*l + # CGAGA
    #     [106.73295]*l + # GAGAC
    #     [100.980286]*l + # AGACT
    #     [76.73934]*l + # GACTA
    #     [100.237816]*l) # ACTAA
    # read = 'GCTAGCACTGATCCGAGACT' # TODO spaeter vlt: sollte 5' -> 3' sein

    # for signal, read in records:
    #     print(feedSegmentation(signal, read, pipe))

    # read is 5' -> 3'
    # 96.8,96.8,96.8,96.8,96.8,96.8,96.8,96.8,96.8,96.8,84.4,84.4,84.4,84.4,84.4,84.4,84.4,84.4,84.4,84.4,95.1,95.1,95.1,95.1,95.1,95.1,95.1,95.1,95.1,95.1,117.7,117.7,117.7,117.7,117.7,117.7,117.7,117.7,117.7,117.7,87.5,87.5,87.5,87.5,87.5,87.5,87.5,87.5,87.5,87.5,79.5,79.5,79.5,79.5,79.5,79.5,79.5,79.5,79.5,79.5,80.2,80.2,80.2,80.2,80.2,80.2,80.2,80.2,80.2,80.2,82.1,82.1,82.1,82.1,82.1,82.1,82.1,82.1,82.1,82.1,118.7,118.7,118.7,118.7,118.7,118.7,118.7,118.7,118.7,118.7,118.6,118.6,118.6,118.6,118.6,118.6,118.6,118.6,118.6,118.6,99.5,99.5,99.5,99.5,99.5,99.5,99.5,99.5,99.5,99.5,83.9,83.9,83.9,83.9,83.9,83.9,83.9,83.9,83.9,83.9,68.4,68.4,68.4,68.4,68.4,68.4,68.4,68.4,68.4,68.4,111.2,111.2,111.2,111.2,111.2,111.2,111.2,111.2,111.2,111.2,103.0,103.0,103.0,103.0,103.0,103.0,103.0,103.0,103.0,103.0,128.6,128.6,128.6,128.6,128.6,128.6,128.6,128.6,128.6,128.6,106.7,106.7,106.7,106.7,106.7,106.7,106.7,106.7,106.7,106.7,100.9,100.9,100.9,100.9,100.9,100.9,100.9,100.9,100.9,100.9,76.7,76.7,76.7,76.7,76.7,76.7,76.7,76.7,76.7,76.7,100.2,100.2,100.2,100.2,100.2,100.2,100.2,100.2,100.2,100.2 TCAGAGCCTAGTTACGATCG

    # read is 3' -> 5'
    # 96.8,96.8,96.8,96.8,96.8,96.8,96.8,96.8,96.8,96.8,84.4,84.4,84.4,84.4,84.4,84.4,84.4,84.4,84.4,84.4,95.1,95.1,95.1,95.1,95.1,95.1,95.1,95.1,95.1,95.1,117.7,117.7,117.7,117.7,117.7,117.7,117.7,117.7,117.7,117.7,87.5,87.5,87.5,87.5,87.5,87.5,87.5,87.5,87.5,87.5,79.5,79.5,79.5,79.5,79.5,79.5,79.5,79.5,79.5,79.5,80.2,80.2,80.2,80.2,80.2,80.2,80.2,80.2,80.2,80.2,82.1,82.1,82.1,82.1,82.1,82.1,82.1,82.1,82.1,82.1,118.7,118.7,118.7,118.7,118.7,118.7,118.7,118.7,118.7,118.7,118.6,118.6,118.6,118.6,118.6,118.6,118.6,118.6,118.6,118.6,99.5,99.5,99.5,99.5,99.5,99.5,99.5,99.5,99.5,99.5,83.9,83.9,83.9,83.9,83.9,83.9,83.9,83.9,83.9,83.9,68.4,68.4,68.4,68.4,68.4,68.4,68.4,68.4,68.4,68.4,111.2,111.2,111.2,111.2,111.2,111.2,111.2,111.2,111.2,111.2,103.0,103.0,103.0,103.0,103.0,103.0,103.0,103.0,103.0,103.0,128.6,128.6,128.6,128.6,128.6,128.6,128.6,128.6,128.6,128.6,106.7,106.7,106.7,106.7,106.7,106.7,106.7,106.7,106.7,106.7,100.9,100.9,100.9,100.9,100.9,100.9,100.9,100.9,100.9,100.9,76.7,76.7,76.7,76.7,76.7,76.7,76.7,76.7,76.7,76.7,100.2,100.2,100.2,100.2,100.2,100.2,100.2,100.2,100.2,100.2 GCTAGCATTGATCCGAGACT


if __name__ == '__main__':
    main()