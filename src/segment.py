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
    signal = str(signal.tolist()).replace(' ', '').replace('[', '').replace(']', '')
    cookie = f'{signal} {read}\n'
    # transfer data to bytes - needed in Python 3
    cookie = bytes(cookie, 'UTF-8')
    # feed cookie to segmentation
    pipe.stdin.write(cookie)
    pipe.stdin.flush()
    # prepare segmentation and return
    crumbs = pipe.stdout.readline().strip().decode('UTF-8')[:-1]
    # print(crumbs)
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
    l = 5

    signal = np.array(
        [96.86647]*l + # AAGCT
        [84.4253]*l + # AGCTA
        [95.15087]*l + # GCTAG
        [117.71723]*l + # CTAGC
        [87.507645]*l + # TAGCA
        [79.56085]*l + # AGCAT
        [80.2712]*l + # GCATT
        [82.19202]*l + # CATTG
        [118.77031]*l + # ATTGA
        [118.66758]*l + # TTGAT
        [99.587555]*l + # TGATC
        [83.94988]*l + # GATCC
        [68.43068]*l + # ATCCG
        [111.2103]*l + # TCCGA
        [103.08109]*l + # CCGAG
        [128.6229]*l + # CGAGA
        [106.73295]*l + # GAGAC
        [100.980286]*l + # AGACT
        [76.73934]*l + # GACTA
        [100.237816]*l) # ACTAA
    read = 'GCTAGCACTGATCCGAGACT' # TODO spaeter vlt: sollte 5' -> 3' sein

    # for signal, read in records:
    #     print(feedSegmentation(signal, read, pipe))

    # read is 5' -> 3'
    # 96.8,96.8,96.8,96.8,96.8,96.8,96.8,96.8,96.8,96.8,84.4,84.4,84.4,84.4,84.4,84.4,84.4,84.4,84.4,84.4,95.1,95.1,95.1,95.1,95.1,95.1,95.1,95.1,95.1,95.1,117.7,117.7,117.7,117.7,117.7,117.7,117.7,117.7,117.7,117.7,87.5,87.5,87.5,87.5,87.5,87.5,87.5,87.5,87.5,87.5,79.5,79.5,79.5,79.5,79.5,79.5,79.5,79.5,79.5,79.5,80.2,80.2,80.2,80.2,80.2,80.2,80.2,80.2,80.2,80.2,82.1,82.1,82.1,82.1,82.1,82.1,82.1,82.1,82.1,82.1,118.7,118.7,118.7,118.7,118.7,118.7,118.7,118.7,118.7,118.7,118.6,118.6,118.6,118.6,118.6,118.6,118.6,118.6,118.6,118.6,99.5,99.5,99.5,99.5,99.5,99.5,99.5,99.5,99.5,99.5,83.9,83.9,83.9,83.9,83.9,83.9,83.9,83.9,83.9,83.9,68.4,68.4,68.4,68.4,68.4,68.4,68.4,68.4,68.4,68.4,111.2,111.2,111.2,111.2,111.2,111.2,111.2,111.2,111.2,111.2,103.0,103.0,103.0,103.0,103.0,103.0,103.0,103.0,103.0,103.0,128.6,128.6,128.6,128.6,128.6,128.6,128.6,128.6,128.6,128.6,106.7,106.7,106.7,106.7,106.7,106.7,106.7,106.7,106.7,106.7,100.9,100.9,100.9,100.9,100.9,100.9,100.9,100.9,100.9,100.9,76.7,76.7,76.7,76.7,76.7,76.7,76.7,76.7,76.7,76.7,100.2,100.2,100.2,100.2,100.2,100.2,100.2,100.2,100.2,100.2 TCAGAGCCTAGTTACGATCG

    # read is 3' -> 5'
    # 96.8,96.8,96.8,96.8,96.8,96.8,96.8,96.8,96.8,96.8,84.4,84.4,84.4,84.4,84.4,84.4,84.4,84.4,84.4,84.4,95.1,95.1,95.1,95.1,95.1,95.1,95.1,95.1,95.1,95.1,117.7,117.7,117.7,117.7,117.7,117.7,117.7,117.7,117.7,117.7,87.5,87.5,87.5,87.5,87.5,87.5,87.5,87.5,87.5,87.5,79.5,79.5,79.5,79.5,79.5,79.5,79.5,79.5,79.5,79.5,80.2,80.2,80.2,80.2,80.2,80.2,80.2,80.2,80.2,80.2,82.1,82.1,82.1,82.1,82.1,82.1,82.1,82.1,82.1,82.1,118.7,118.7,118.7,118.7,118.7,118.7,118.7,118.7,118.7,118.7,118.6,118.6,118.6,118.6,118.6,118.6,118.6,118.6,118.6,118.6,99.5,99.5,99.5,99.5,99.5,99.5,99.5,99.5,99.5,99.5,83.9,83.9,83.9,83.9,83.9,83.9,83.9,83.9,83.9,83.9,68.4,68.4,68.4,68.4,68.4,68.4,68.4,68.4,68.4,68.4,111.2,111.2,111.2,111.2,111.2,111.2,111.2,111.2,111.2,111.2,103.0,103.0,103.0,103.0,103.0,103.0,103.0,103.0,103.0,103.0,128.6,128.6,128.6,128.6,128.6,128.6,128.6,128.6,128.6,128.6,106.7,106.7,106.7,106.7,106.7,106.7,106.7,106.7,106.7,106.7,100.9,100.9,100.9,100.9,100.9,100.9,100.9,100.9,100.9,100.9,76.7,76.7,76.7,76.7,76.7,76.7,76.7,76.7,76.7,76.7,100.2,100.2,100.2,100.2,100.2,100.2,100.2,100.2,100.2,100.2 GCTAGCATTGATCCGAGACT
    np.set_printoptions(formatter={'float': lambda x: "{0:0.1f}".format(x)})
    print(np.array_repr(signal).replace('\n', '').replace(' ', ''), read)
    segmentprobs = feedSegmentation(signal, read, pipe)[1:-1]
    stopFeeding(pipe)

    fig, ax1 = plt.subplots(figsize=(12,8), dpi=150)
    ax1.plot(signal, color='red', label='signal')
    ax1.set_ylim((0, max(signal)+10))
    ax1.set_ylabel('Signal pico Ampere')
    ax1.grid(True, 'both', 'x')
    ax2 = ax1.twinx()
    ax2.plot(np.arange(0, len(segmentprobs)), segmentprobs, '.', label='probabilities', color='grey')
    ax2.fill_between(np.arange(0, len(segmentprobs)), min(segmentprobs), segmentprobs, color='grey', alpha=0.6)
    ax2.set_ylim((-20, 40))
    ax2.set_xlim((0, len(signal)))
    ax2.set_xlabel('Interspace Between Signal Points')
    ax2.set_ylabel('Transition Probability')
    ax2.set_xticks(np.arange(0, len(segmentprobs), l))
    h1, l1 = ax1.get_legend_handles_labels()
    h2, l2 = ax2.get_legend_handles_labels()
    plt.legend([*h1, *h2], [*l1, *l2])
    plt.tight_layout()
    plt.savefig('output.pdf')

    with open('output.txt', 'w') as w:
        w.write(f"Segmentpobabilies (log):\n{segmentprobs}\nRead:\n{read}")

if __name__ == '__main__':
    main()