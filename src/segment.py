#!/usr/bin/env python
# author: Jannes Spangenberg
# e-mail: jannes.spangenberg@uni-jena.de
# github: https://github.com/JannesSP
# website: https://jannessp.github.io

from argparse import ArgumentDefaultsHelpFormatter, ArgumentParser, Namespace
from os.path import join, dirname
from os import name
import numpy as np
from subprocess import Popen, PIPE
from pathlib import Path
import matplotlib.pyplot as plt

TERM_STRING = "$"

def parse() -> Namespace:
    parser = ArgumentParser(
        formatter_class=ArgumentDefaultsHelpFormatter
    )
    parser.add_argument('ONTpath', type=str, help='Folder containing ONT raw data, either {fast5, slow5, blow5 or pod5}')
    parser.add_argument('-r', '--recursive', action='store_true', type=bool, help='Look for ONT data recursively')
    parser.add_argument('outpath', type=str, help='Output directory')
    parser.add_argument('-p', '--processes', type=int, default=1, help='Number of processes for multiprocessing')
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
    signal = str(signal.tolist()).replace(' ', '')[1:-1]
    cookie = f'{signal} {read}\n'
    # transfer data to bytes - needed in Python 3
    cookie = bytes(cookie, 'UTF-8')
    # feed cookie to segmentation
    pipe.stdin.write(cookie)
    pipe.stdin.flush()
    # prepare segmentation and return
    crumbs = pipe.stdout.readline().strip().decode('UTF-8')[:-1]
    print(crumbs)
    crumbs = np.array(crumbs.split(','), dtype=float)
    return crumbs

def stopFeeding(pipe : Popen) -> None:
    pipe.stdin.write(bytes(TERM_STRING, 'UTF-8'))
    pipe.stdin.close()
    pipe.stdout.close()

def getFiles(filepath : str, rec : bool) -> list:

    func = {True : Path(filepath).rglob, False : Path(filepath).glob}[rec]
    files = []

    for path in func('*.fast5'):
        files.append(path.name)

    # for path in func('*.slow5'):
    #     print(path.name)

    # for path in func('*.blow5'):
    #     print(path.name)

    # for path in func('*.pod5'):
    #     print(path.name)

    return files

def main() -> None:
    # args = parse()
    # inp = args.ONTpath
    # out = args.outpath

    # assert exists(inp), f'Input path {inp} does not exist!'
    # assert exists(out), f'Output path {out} does not exist!'

    # # TODO read raw data from ont files
    # recursively = args.recursive
    # processes = args.processes

    CPP_SCRIPT = join(dirname(__file__), 'segment')
    if name == 'nt': # check for windows
        CPP_SCRIPT+='.exe'
    pipe = openCPPScript(CPP_SCRIPT)




    # #########################
    # Test
    # ret = feedSegmentation(np.array([107.2,108.0,108.9,111.2,105.7,104.3,107.1,105.7]), 'CAAAAA', pipe)
    # print(ret)
    l = 30

    signal = np.array(
        [96.8]*l + # AAGCT
        [84.4]*l + # AGCTA
        [95.1]*l + # GCTAG
        [117.7]*l + # CTAGC
        [87.5]*l + # TAGCA
        [79.5]*l + # AGCAT
        [80.2]*l + # GCATT
        [82.1]*l + # CATTG
        [118.7]*l + # ATTGA
        [118.6]*l + # TTGAT
        [99.5]*l + # TGATC
        [83.9]*l + # GATCC
        [68.4]*l + # ATCCG
        [111.2]*l + # TCCGA
        [103.0]*l + # CCGAG
        [128.6]*l + # CGAGA
        [106.7]*l + # GAGAC
        [100.9]*l + # AGACT
        [76.7]*l + # GACTA
        [100.2]*l) # ACTAA
    read = 'GCTAGCATTGATCCGAGACT'[::-1] # sollte 5' -> 3' sein

    # for signal, read in records:
    #     print(feedSegmentation(signal, read, pipe))

    print(signal)
    print(read)
    segmentprobs = feedSegmentation(signal, read, pipe)[1:-1]

    fig, ax1 = plt.subplots()
    ax2 = ax1.twinx()
    ax1.plot(np.arange(0, len(segmentprobs)), segmentprobs, 'o', label='probabilities', color='grey', )
    ax1.fill_between(np.arange(0, len(segmentprobs)), min(segmentprobs), segmentprobs, color='grey', alpha=0.6)
    ax1.set_ylim((min(segmentprobs), 150))
    ax1.set_xlim((0, len(signal)))
    ax1.set_xlabel('Interspace Between Signal Points')
    ax1.set_ylabel('Transition Probability')
    ax1.set_xticks(np.arange(0, len(segmentprobs), l))
    ax1.grid(True, 'both', 'x')
    ax2.plot(signal, color='red', label='signal')
    ax2.set_ylim((0, max(signal)+10))
    ax2.set_ylabel('Signal pico Ampere')
    fig.legend()
    plt.tight_layout()
    plt.show()

    stopFeeding(pipe)

if __name__ == '__main__':
    main()