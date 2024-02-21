#!/usr/bin/env python
# author: Jannes Spangenberg
# e-mail: jannes.spangenberg@uni-jena.de
# github: https://github.com/JannesSP
# website: https://jannessp.github.io

TERM_STRING = "$"

from subprocess import PIPE, Popen
import re
import numpy as np
from pathlib import Path
import gzip
# import pandas as pd
from Bio import SeqIO
from os.path import splitext

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

def loadFastx(path : str) -> dict:
    '''
    Returns
    -------
    readDict : dict
        {readid : read}
    '''
    readDict = {}
    if path.endswith(".gz"):
        path = gzip.open(path, "rt")
    for record in SeqIO.parse(path, splitext(path)[1]):
        readDict[record.id] = str(record.seq)
    return readDict

# def readNanoPolyA(file : str) -> pd.DataFrame:
#     return pd.read_csv(file, sep='\t')

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

def openCPPScriptParamsTrain(cpp_script : str, params : dict) -> Popen:
    '''
    Open cpp script with Popen.

    Parameters
    ----------
    cpp_script : str
        Path of script
    params : dict
        {str : float}

    Returns
    -------
    subprocess : Popen
    '''



    script = [cpp_script]
    for param in params:
        script.extend([f"-{param}", str(params[param])])
    script.append("--atrain")
    script=" ".join(script)
    # print("Popen call:", script)
    return Popen(script, shell=True, stdout=PIPE, stdin=PIPE, stderr=PIPE)

def openCPPScriptParams(cpp_script : str, params : dict) -> Popen:
    '''
    Open cpp script with Popen.

    Parameters
    ----------
    cpp_script : str
        Path of script
    params : dict
        {str : float}

    Returns
    -------
    subprocess : Popen
    '''
    script = [cpp_script]
    for param in params:
        script.extend([f"-{param}", str(params[param])])
    script=" ".join(script)
    print("Popen call:", script)
    return Popen(script, shell=True, stdout=PIPE, stdin=PIPE, stderr=PIPE)

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
    cookie = f"{str(signal.tolist()).replace(' ', '').replace('[', '').replace(']', '')} {read}\n"
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
    print("Output length", len(output), len(probs))
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