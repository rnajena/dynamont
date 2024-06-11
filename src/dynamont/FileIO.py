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
from os.path import join
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

class SegmentationError(Exception):
    """Raised when no segmentation was calculated for a read"""
    def __init__(self, read) -> None:
        self.read = read
        self.message = f"No segmentation calculated for {read}"
        super().__init__(self.message)

def readKmerModels(filepath : str) -> dict:
    '''
    Returns
    -------
    models : dict
        format: {kmer : (mean, stdev)}
    '''
    df = pd.read_csv(filepath, sep='\t')
    models = pd.Series(zip(df.level_mean.values, df.level_stdv.values), index=df.kmer).to_dict()
    return models

def writeKmerModels(filepath : str, kmerModels : dict) -> dict:
    '''
    Writes kmer models to a file
    '''
    with open(filepath, 'w') as w:
        w.write('kmer\tlevel_mean\tlevel_stdv\n')
        for kmer in kmerModels:
            mean = kmerModels[kmer][0]
            # safety for simulated data with stdev = 0
            stdev = kmerModels[kmer][1] if kmerModels[kmer][1] != 0 else 0.0001
            w.write(f'{kmer}\t{mean}\t{stdev}\n')

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

def loadFastx(path : str, quality : int = None) -> dict:
    '''
    Returns
    -------
    readDict : dict
        {readid : read}
    quality : int
        add a quality value to filter out reads that have a mean quality value below the given threshold
    '''
    readDict = {}
    if path.endswith(".gz"):
        path = gzip.open(path, "rt")
    for record in SeqIO.parse(path, splitext(path)[1][1:]):
        if quality is not None:
            if not np.mean(record.letter_annotations["phred_quality"]) >= quality:
                continue
        readDict[record.id] = str(record.seq)
    return readDict

# def readNanoPolyA(file : str) -> pd.DataFrame:
#     return pd.read_csv(file, sep='\t')

def openCPPScript(cpp_script : str) -> Popen:
    '''
    Popen([cpp_script], shell=True, stdout=PIPE, stdin=PIPE, stderr=PIPE)
    '''
    print("Popen call:", cpp_script)
    return Popen(cpp_script, shell=True, stdout=PIPE, stdin=PIPE, stderr=PIPE)

def openCPPScriptATrain(cpp_script : str, params : dict) -> Popen:
    '''
    Open cpp script with transition params in transition-training mode

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
    return openCPPScript(script)

def openCPPScriptTrain(cpp_script : str, params : dict, model_file : str, minSegLen : int) -> Popen:
    '''
    Open cpp script with transition params in training mode

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
    script.append("--train")
    script.append(f"--model {model_file}")
    script.append(f"--minSegLen {minSegLen}")
    script=" ".join(script)
    return openCPPScript(script)

def openCPPScriptCalcZ(cpp_script : str, params : dict, model_file : str = None, minSegLen : int = None) -> Popen:
    '''
    Open cpp script with transition params in calculate-Z mode

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
    if model_file is not None:
        script.append(f"--model {model_file}")
    if minSegLen is not None:
        script.append(f"--minSegLen {minSegLen}")
    script=" ".join(script)
    return openCPPScript(script)

def calcZ(signal : np.ndarray, read : str, params : dict, script : str, model_file : str = None, minSegLen : int = None) -> float:
    pipe = openCPPScriptCalcZ(script, params, model_file, minSegLen)
    feedPipe(signal, read, pipe)
    Z = float(pipe.stdout.readline().strip().decode('UTF-8'))
    stopFeeding(pipe)
    return Z

def openCPPScriptParams(cpp_script : str, params : dict) -> Popen:
    '''
    Open cpp script with transition params

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
    return openCPPScript(script)

def feedPipe(signal : np.ndarray, read : str, pipe : Popen) -> None:
    '''
    Feeds signal and read to pipe

    Params
    ------
    signal : np.ndarray
    read : str
    pipe : Popen
    '''
    # prepare cookie for segmentation
    signal = str(np.around(signal.tolist(), 3).tolist()).replace(' ', '').replace('[', '').replace(']', '')
    cookie = f"{signal}\n{read}\n"
    # with open("lastCookie.txt", "w") as w:
    #     w.write(cookie)
    # print(cookie)
    cookie = bytes(cookie, 'UTF-8')
    pipe.stdin.write(cookie)
    pipe.stdin.flush()

# # https://stackoverflow.com/questions/32570029/input-to-c-executable-python-subprocess
# def trainTransitions(signal : np.ndarray, read : str, params : dict, script : str) -> tuple:
#     '''
#     Parse & feed signal & read to the C++ segmentation script.

#     Parameters
#     ----------
#     signal : np.ndarray
#     read : str
#         in 3' -> 5' direction
#     stream
#         Open stdin stream of the C++ segmentation algorithm

#     Returns
#     -------
#     params : dict
#         {str : float}
#     Z : float
#     '''
#     pipe = openCPPScriptATrain(script, params)
#     feedPipe(signal, read, pipe)
#     output = pipe.stdout.readline().strip().decode('UTF-8')
#     params = {param.split(":")[0] : float(param.split(":")[1]) for param in output.split(";")}

#     Z = float(pipe.stdout.readline().strip().decode('UTF-8').split(':')[1])
#     stopFeeding(pipe)

#     return params, Z

# https://stackoverflow.com/questions/32570029/input-to-c-executable-python-subprocess
def trainTransitionsEmissions(signal : np.ndarray, read : str, params : dict, script : str, model_file : str, minSegLen : int) -> tuple:
    '''
    Parse & feed signal & read to the C++ segmentation script.

    Parameters
    ----------
    signal : np.ndarray
    read : str
        in 3' -> 5' direction
    stream
        Open stdin stream of the C++ segmentation algorithm

    Returns
    -------
    segments : np.ndarray

    params : dict
        {str : float}
    Z : float
    '''
    pipe = openCPPScriptTrain(script, params, model_file, minSegLen)
    # print(read)
    feedPipe(signal, read, pipe)
    # print(pipe.stdout.readline().strip().decode('UTF-8'))
    # print(pipe.stdout.readline().strip().decode('UTF-8'))
    
    trainedParams = pipe.stdout.readline().strip().decode('UTF-8')
    # print(trainedParams)
    # print(pipe.stderr.readline())

    # first the transition parameters
    try:
        params = {param.split(":")[0] : float(param.split(":")[1]) for param in trainedParams.split(";")}
    except IndexError:
        print("ERROR", trainedParams)
        exit(1)

    # then updated emission updated
    trainedParams = pipe.stdout.readline().strip().decode('UTF-8')
    #               kmer                    mean                                        stdev
    newModels = {param.split(":")[0] : (float(param.split(":")[1].split(",")[0]), float(param.split(":")[1].split(",")[1])) for param in trainedParams.split(";")[:-1]}
    # print(trainedParams)
    # print(newModels)
    # then Z
    Z = float(pipe.stdout.readline().strip().decode('UTF-8').split(':')[1])
    
    # print(pipe.stderr.readline())
    # print(pipe.stderr.readline())

    # then segmentation
    segmentation = pipe.stdout.readline().strip().decode('UTF-8')
    segments = formatSegmentationOutput(segmentation, len(signal), len(read))

    stopFeeding(pipe)

    return segments, params, newModels, Z

# https://stackoverflow.com/questions/32570029/input-to-c-executable-python-subprocess
def feedSegmentation(signal : np.ndarray, read : str, script : str, params : dict = None) -> np.ndarray:
    '''
    Parse & feed signal & read to the C++ segmentation script.

    Parameters
    ----------
    segmentation : np.ndarray
        numpy array of marking the borders
        format: [[start : int, end : int, readpos : int, state : str] ...]
            start : signal idx
            end : signal idx
            readpos : base position within 5' - 3' read
            state: {'M', 'I', 'D'}
                'M' : match
                'I' : insertion of a base
                'D' : deletion of a base
    read : str
        in 3' -> 5' direction
    script : str
        script file name

    Returns
    -------
    segmentation : np.ndarray
        numpy array of marking the borders
        format: [[start : int, end : int, readpos : int, state : str] ...]
            start : signal idx
            end : signal idx
            readpos : base position within 5' - 3' read
            state: {'M', 'I', 'D'}
                'M' : match
                'I' : insertion of a base
                'D' : deletion of a base
    '''
    # probabilities : np.ndarray
    #     in 3' -> 5' orientation
    #     numpy array of floats showing the border probability
    
    if params is None:
        pipe = openCPPScript(script)
    else:
        pipe = openCPPScriptParams(script, params)
    # print(pipe.args, len(signal), len(read))
    feedPipe(signal, read, pipe)
    # receive segmentation result
    output = pipe.stdout.readline().strip().decode('UTF-8')
    # format output into np.ndarray
    segments = formatSegmentationOutput(output, len(signal), len(read))
    stopFeeding(pipe)

    return segments #, probs

def formatSegmentationOutput(output : str, signalLength : int, readLength : int) -> np.ndarray:
    '''
    Receives the segmentation output from CPP script and returns it as a np.ndarray
    
    Parameters
    ----------
    output : str
        Segmentation output from CPP script
    signalLength : int
        Length of the read signal
    readLength : int
        Length of the read or number of nucleotides
    
    Returns
    -------
    segmentation : np.ndarray
        numpy array of marking the borders
        format: [[start : int, end : int, readpos : int, state : str] ...]
            start : signal idx
            end : signal idx
            readpos : base position within 5' - 3' read
            state: {'M', 'I', 'D'}
                'M' : match
                'I' : insertion of a base
                'D' : deletion of a base
    '''
    output = re.findall('\D[\d,]+', output)
    segments = []
    for i in range(len(output)):
        state = output[i][0]
        basepos, start = output[i][1:].split(',')
        end = output[i+1][1:].split(',')[1] if i < len(output)-1 else signalLength
        # convert basepos to 5' -> 3' direction
        segments.append([int(start), int(end), readLength - int(basepos) - 1, state])
        
    return np.array(segments, dtype=object)

def genSegmentTable(signal : np.ndarray, read : str, script : str, readid : str) -> str:
    '''
    Transform segmentation into a string to write into a file

    Parameters
    ----------
    signal : np.ndarray
    read : str
        in 3' -> 5' direction
    script : str
        script file name
    readid : str

    Returns
    -------
    table : str
        string to write in a segmentation result file
    '''
    segmentation = feedSegmentation(signal, read, script)
    table = ''
    for segment in segmentation:
        table += readid + ',' + ','.join(segment) + '\n'
    return table

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
    '''
    Returns
    -------
    polya : dict
        {readid : [polyastart, polyaend]}
    '''
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
