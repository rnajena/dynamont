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
from os.path import join
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from multiprocessing import Queue

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
    # try:
    models = pd.Series(zip(df.level_mean.values, df.level_stdv.values), index=df.kmer).to_dict()
    # except AttributeError:
    #     df['level_stdv'] = 1.0
    #     models = pd.Series(zip(df.level_mean.values, df.level_stdv.values), index=df.kmer).to_dict()
    return models

def writeKmerModels(filepath : str, kmerModels : dict) -> dict:
    '''
    Writes kmer models to a file
    '''
    with open(filepath, 'w') as w:
        w.write('kmer\tlevel_mean\tlevel_stdv\n')
        for kmer in kmerModels:
            # happens in indel model, when kmer is deleted in all paths
            mean = kmerModels[kmer][0] # if not np.isnan(kmerModels[kmer][0]) else backUpModel[kmer][0]
            # safety to very small stdev TODO find another way to fix this? happens in real data, when segment very small - nucleotides do not match signal
            stdev = kmerModels[kmer][1] # if kmerModels[kmer][1] > 0 else backUpModel[kmer][1]

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

# def loadBasecalls(bamFile : str, quality : int = None) -> dict:
#     '''
#     Parameters
#     ----------
#     bamfile : str
#         basecalls in .bam format from the dorado basecaller
#     quality : int
#         add a quality value to filter out reads that have a mean quality value below the given threshold

#     Returns
#     -------
#     readDict : dict
#         {readid : (read, transcript_start)}
#     '''
#     readDict = {}
#     with pysam.AlignmentFile(bamFile, "r" if bamFile.endswith('.sam') else "rb", check_sq=False) as samfile:
#         for read in samfile.fetch(until_eof=True):
#             seq = read.query_sequence
#             rid = read.query_name
#             ts = read.get_tag("ts")
            
#             # skip reads below quality threshold
#             if quality:
#                 avg_qual = sum(read.query_qualities) / len(seq)
#                 if avg_qual < quality:
#                     continue

#             readDict[rid] = (seq, ts)
#     return readDict

# def indexBasecalls(bamFile : str, rawFiles : list) -> dict:
#     idx = {}
#     with pysam.AlignmentFile(bamFile, "r" if bamFile.endswith('.sam') else "rb", check_sq=False) as samfile:
#         for basecalled_read in samfile.fetch(until_eof=True):
#             rid = basecalled_read.query_name
#             for file in rawFiles:
#                 if rid in read5.read(file).getReads():
#                     idx[rid] = file
#     return idx

def openCPPScript(cpp_script : str) -> Popen:
    '''
    Popen([cpp_script], shell=True, stdout=PIPE, stdin=PIPE)
    '''
    # print("Popen call:", cpp_script)
    return Popen(cpp_script, stdout=PIPE, stdin=PIPE, text=True)

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
    script.extend(["--train", "--model", model_file, "--minSegLen", str(minSegLen)])
    return openCPPScript(script)

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
        if param == 'p':
            script.extend(['-p'])
        else:
            script.extend([f"-{param}", str(params[param])])
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
        script.extend(["--model", model_file])
    if minSegLen is not None:
        script.extend(["--minSegLen",  str(minSegLen)])
    return openCPPScript(script)

def calcZ(signal : np.ndarray, read : str, params : dict, script : str, model_file : str = None, minSegLen : int = None, readid : str = None) -> float:
    pipe = openCPPScriptCalcZ(script, params, model_file, minSegLen)
    try:
        Z = float(feedPipe(signal, read, pipe).strip())
    except:
        return readid
    return Z

def feedPipe(signal : np.ndarray, read : str, pipe : Popen) -> str:
    '''
    Feeds signal and read to pipe with an immediate stop signal afterwards.
    Uses communicate.

    Params
    ------
    signal : np.ndarray
    read : str
    pipe : Popen

    Returns
    -------
    result : str
    '''
    # print("Debug", "Signal length", len(signal), "Read length", len(read))

    # prepare cookie for segmentation
    # signal = str(np.around(signal.tolist(), 5).tolist()).replace(' ', '').replace('[', '').replace(']', '')
    signal = str(signal.tolist()).replace(' ', '').replace('[', '').replace(']', '')
    cookie = f"{signal}\n{read}\n{TERM_STRING}\n{TERM_STRING}\n"
    result = pipe.communicate(input=cookie)[0]
    stopFeeding(pipe)
    # with open("lastCookie.txt", "w") as w:
    #     w.write(cookie)
    return result

# https://stackoverflow.com/questions/32570029/input-to-c-executable-python-subprocess
def trainTransitionsEmissions(signal : np.ndarray, read : str, params : dict, script : str, model_file : str, minSegLen : int, readid : str) -> tuple|str:
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
    try:
        result = feedPipe(signal, read, pipe).split('\n')
        transitionParams, modelParams, Z, _ = result
    except:
        print(f"ERROR while extracting result in {readid}, {result}")
        # with open("failed_input.txt", 'w') as w:
        #     w.write(str(signal.tolist()).replace(' ', '').replace('[', '').replace(']', '') + '\n' + read + '\n')
        # raise SegmentationError(readid)
        return readid
    try:
        params = {param.split(":")[0] : float(param.split(":")[1]) for param in transitionParams.split(";")}
    except:
        print(f"ERROR while extracting transitions params in {readid}", transitionParams)
        # with open("failed_input.txt", 'w') as w:
        #     w.write(str(signal.tolist()).replace(' ', '').replace('[', '').replace(']', '') + '\n' + read + '\n')
        # raise SegmentationError(readid)
        return readid
    # then updated emission updated
    #               kmer                    mean                                        stdev
    newModels = {param.split(":")[0] : (float(param.split(":")[1].split(",")[0]), float(param.split(":")[1].split(",")[1])) for param in modelParams.split(";")[:-1]}
    Z = float(Z.split(':')[-1])
    return params, newModels, Z

# https://stackoverflow.com/questions/32570029/input-to-c-executable-python-subprocess
def feedSegmentation(signal : np.ndarray, read : str, script : str, params : dict = None) -> tuple:
    '''
    Parse & feed signal & read to the C++ segmentation script.
    Opens and closes a pipe to the given script.

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
    params : dict
        params for the c++ script in short form: e.g. {-m : "model_path"}

    Returns
    -------
    segmentation : np.ndarray
        numpy array of marking the borders
        format: [[start : int, end : int, readpos : int, state : str] ...]
            start : signal idx
            end : signal idx
            readpos : base position within 5' - 3' read
            state: states of the HMM used
            polish: if available a polish of the motif
    '''
    # probabilities : np.ndarray
    #     in 3' -> 5' orientation
    #     numpy array of floats showing the border probability

    if params is None:
        pipe = openCPPScript(script)
    else:
        pipe = openCPPScriptParams(script, params)
    result = feedPipe(signal, read, pipe)

    try:
        output, probs, _ = result.split('\n')
        probs = np.array(list(map(float, probs.split(',')[:-1])))[1:]
    except:
        try:
            output, _  = result.split('\n')
            probs = None
        except:
            print(signal)
            print(read)
            with open("failed_input.txt", "w") as w:
                w.write(str(signal.tolist()).replace(' ', '').replace('[', '').replace(']', ''))
                w.write('\n')
                w.write(read)
                w.write('\n')
            exit(1)

    # receive segmentation result
    # format output into np.ndarray
    # try:
    segments = formatSegmentationOutput(output, len(signal), read[::-1])
    # except:
        # some error during segmentation - no proper output
        # return np.array([])
    return segments , probs

globalPipe : Popen = None

def feedSegmentationAsynchronous(CPP_SCRIPT : str, params : dict, signal : np.ndarray, read : str, readid : str, queue : Queue) -> None:
    '''
    Parse & feed signal & read to the C++ segmentation script.
    Needs an open pipe for communication.

    Parameters
    ----------
    signal : np.ndarray
        in 4' -> 5' direction
    read : str
        in 3' -> 5' direction
    readid : str
    pipe : Popen
    queue : Queue
    '''
    global globalPipe
    # terminate
    if read == TERM_STRING:
        stopFeeding(globalPipe)
        return None
    if globalPipe is None:
        globalPipe = openCPPScriptParams(CPP_SCRIPT, params)

    output = feedPipeAsynchronous(signal, read, globalPipe)
    try:
        segments = formatSegmentationOutput(output, len(signal), read)
    except Exception as e:
        print("Exception in readid", readid, e)
        return None
    queue.put(formatSegmentation(readid, segments))

def feedPipeAsynchronous(signal : np.ndarray, read : str, pipe : Popen) -> str:
    '''
    Feeds signal and read to pipe without a stop signal.
    Uses stdin.write
    Then reads and returns one line with stdout.readline().

    Params
    ------
    signal : np.ndarray
    read : str
    pipe : Popen

    Returns
    -------
    result : str
    '''
    # prepare cookie for segmentation
    # signal = str(np.around(signal.tolist(), 5).tolist()).replace(' ', '').replace('[', '').replace(']', '')
    signal = str(signal.tolist()).replace(' ', '').replace('[', '').replace(']', '')
    cookie = f"{signal}\n{read}\n"
    # with open("lastCookie.txt", "w") as w:
    #     w.write(cookie)
    pipe.stdin.write(cookie)
    pipe.stdin.flush()
    result = pipe.stdout.readline()
    return result

def formatSegmentationOutput(output : str, sigLen : int, read : str) -> np.ndarray:
    '''
    Receives the segmentation output from CPP script and returns it as a np.ndarray

    Parameters
    ----------
    output : str
        Segmentation output from CPP script
    sigLen : int
        Length of the read signal
    read : str
        5' -> 3' direction

    Returns
    -------
    segmentation : np.ndarray
        numpy array of marking the borders
        format: [[start : int, end : int, readpos : int, state : str] ...]
            start : signal idx
            end : signal idx
            readpos : base position within 5' - 3' read
            state : match, extend, etc
    '''
    segments = []
    output = output.strip().split(';')[:-1]
    for i, segment in enumerate(output):
        print(segment)
        # split segment state
        state = segment[0]
        segment = segment[1:]

        # get basepos, signal start and polish if existing
        try:
            basepos, start, prob, polish = segment.split(',')
        except ValueError:
            basepos, start, prob = segment.split(',')
            polish = 'NA'

        end = output[i+1][1:].split(',')[1] if i < len(output)-1 else sigLen

        # convert basepos to 5' -> 3' direction
        pos = len(read) - int(basepos) - 1
        segments.append([int(start), int(end), pos, read[pos], read[max(0, pos-2):min(len(read), pos+3)], state, prob, polish])

    return np.array(segments, dtype=object)

def formatSegmentation(readid : str, segmentation : np.ndarray) -> str:
    '''
    Transform segmentation into a string to write into a file

    Parameters
    ----------
    segmentation : np.ndarray
    readid : str

    Returns
    -------
    table : str
        string to write in a segmentation result file
    '''
    # if model is not None:
    #     params = {"m" : model}
    # segmentation = feedSegmentation(signal, read, script, params)
    table = ''
    for segment in segmentation:
        table += readid + ',' + ','.join(list(map(str, segment))) + '\n'
    return table

def stopFeeding(pipe : Popen) -> None:
    try:
        pipe.kill()
    except:
        pass
    return None

    # if pipe is None:
    #     return None
    # if pipe.poll() == 0:
    #     return None
    # else:
    #     try:
    #         pipe.communicate(input=f'{TERM_STRING}\n{TERM_STRING}\n')
    #     except:
    #         pass
    #     pipe.kill()
    #     return None

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
    df = pd.read_csv(file, usecols=['readname', 'polya_start', 'transcript_start', 'qc_tag'], sep='\t')
    # filter for good quality reads, set by nanopolish
    df = df[df['qc_tag'] == "PASS"]
    df = df.astype({'readname' : str, 'polya_start' : int, 'transcript_start' : int})
    return pd.Series([[a,b] for a,b in zip(df.polya_start.values, df.transcript_start.values)], index=df.readname).to_dict()

def plotParameters(param_file : str, outdir : str) -> None:
    import matplotlib
    matplotlib.use('Agg') # Use non-interactive backend (no GUI)
    df = pd.read_csv(param_file, sep=',')
    for column in df:
        if column in ['epoch', 'batch', 'read']:
            continue
        sns.set_theme()        
        sns.lineplot(data=df, x="batch", y=column, hue='epoch')
        plt.title(f"{column} parameter change during training")
        plt.ylabel("Parameter Value")
        print("Savefig: ", join(outdir, f"{column}.pdf"))
        plt.savefig(join(outdir, f"{column}.pdf"))
        plt.close()
