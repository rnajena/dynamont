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
from os.path import join, dirname
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

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
    for record in SeqIO.parse(path, splitext(path)[1][1:]):
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

# https://stackoverflow.com/questions/32570029/input-to-c-executable-python-subprocess
def trainSegmentation(signal : np.ndarray, read : str, params : dict, script : str) -> tuple:
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
    params : dict
        {str : float}
    Z : float
    '''
    pipe = openCPPScriptParamsTrain(script, params)
    # prepare cookie for segmentation
    cookie = f"{str(signal.tolist()).replace(' ', '').replace('[', '').replace(']', '')} {read}\n"
    c = open(join(dirname(__file__), 'last_cookie.txt'), 'w')
    c.write(cookie)
    c.close()
    # transfer data to bytes - needed in Python 3
    cookie = bytes(cookie, 'UTF-8')
    # feed cookie to segmentation
    pipe.stdin.write(cookie)
    pipe.stdin.flush()
    output = pipe.stdout.readline().strip().decode('UTF-8')
    try:
        params = {param.split(":")[0] : float(param.split(":")[1]) for param in output.split(";")}
    except Exception as e:
        print(output)
        print(pipe.stdout.readlines())
        print(e)
        exit(1)

    Z = float(pipe.stdout.readline().strip().decode('UTF-8').split(':')[1])
    stopFeeding(pipe)

    return params, Z

def openCPPScriptParamsCalcZ(cpp_script : str, params : dict) -> Popen:
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
    script.append("--calcZ")
    script=" ".join(script)
    # print("Popen call:", script)
    return Popen(script, shell=True, stdout=PIPE, stdin=PIPE, stderr=PIPE)

def calcZ(signal : np.ndarray, read : str, params : dict, script : str) -> float:
    pipe = openCPPScriptParamsCalcZ(script, params)
    # prepare cookie for segmentation
    cookie = f"{str(signal.tolist()).replace(' ', '').replace('[', '').replace(']', '')} {read}\n"
    c = open(join(dirname(__file__), 'last_cookie.txt'), 'w')
    c.write(cookie)
    c.close()
    # transfer data to bytes - needed in Python 3
    cookie = bytes(cookie, 'UTF-8')
    # feed cookie to segmentation
    pipe.stdin.write(cookie)
    pipe.stdin.flush()
    Z = float(pipe.stdout.readline().strip().decode('UTF-8'))
    stopFeeding(pipe)
    return Z

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
def feedSegmentation(signal : np.ndarray, read : str, pipe : Popen) -> tuple:
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
        numpy array of marking the borders
        format: [[start : int, end : int, readpos : int, state : str] ...]
    probabilities : np.ndarray
        in 3' -> 5' orientation
        numpy array of floats showing the border probability
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
    assert segments.size, f"segments: {segments}\n==\nsignal: {signal}\n==\nread: {read}\n"
    return segments, probs

def stopFeeding(pipe : Popen) -> None:
    pipe.stdin.write(bytes(f'{TERM_STRING} {TERM_STRING}\n', 'UTF-8'))
    pipe.stdin.close()
    pipe.stdout.close()

def readPolyAEnd(file : str) -> dict:
    df = pd.read_csv(file, usecols=['readname', 'transcript_start'], sep='\t')
    df = df.astype({'readname' : str, 'transcript_start' : int})
    # df.set_index('readname', inplace=True)
    return pd.Series(df.transcript_start.values, index=df.readname).to_dict()

def readPolyAStartEnd(file : str) -> dict:
    df = pd.read_csv(file, usecols=['readname', 'polya_start', 'transcript_start'], sep='\t')
    df = df.astype({'readname' : str, 'polya_start' : int, 'transcript_start' : int})
    return pd.Series([[a,b] for a,b in zip(df.polya_start.values, df.transcript_start.values)], index=df.readname).to_dict()

def plotParameters(param_file : str, outdir : str) -> None:
    df = pd.read_csv(param_file, sep=',')
    for column in df:
        if column in ['epoch', 'batch']:
            continue
        sns.set_theme()
        sns.lineplot(data=df, x="batch", y=column, hue='epoch')
        plt.title(f"{column} parameter change during training")
        plt.savefig(join(outdir, f"{column}.pdf"))
        plt.cla()
        plt.close()
