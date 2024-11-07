#!/usr/bin/env python
# author: Jannes Spangenberg
# e-mail: jannes.spangenberg@uni-jena.de
# github: https://github.com/JannesSP
# website: https://jannessp.github.io

from subprocess import PIPE, Popen
import numpy as np
from pathlib import Path
from os.path import join
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from multiprocessing import Queue

def hampelFilter(signal, wSize=3, nSigmas=3):
    """
    Apply Hampel filter to detect and replace outliers in a signal.

    Parameters:
        signal (np.ndarray): The input signal (1D array).
        window_size (int): The size of the sliding window. Defaults to 3.
        n_sigmas (float): The threshold in terms of standard deviations (MAD).
                          Defaults to 3, which is commonly used.

    Returns:
        np.ndarray: The filtered signal with outliers replaced.
    """
    # Copy the signal to avoid modifying the original
    filtSignal = signal.copy()
    
    # Half window size for sliding
    k = 1.4826  # Constant to convert MAD to standard deviation
    hwSize = wSize // 2
    
    # Loop over the signal with a sliding window
    for i in range(hwSize, len(signal) - hwSize):
        # Define the current window around the ith element
        window = signal[i - hwSize : i + hwSize + 1]
        
        # Calculate the median and MAD of the window
        median = np.median(window)
        mad = k * np.median(np.abs(window - median))
        
        # Identify and replace outliers
        if np.abs(signal[i] - median) > nSigmas * mad:
            filtSignal[i] = median  # Replace with the median value
    
    return filtSignal

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

def writeKmerModels(filepath : str, kmerModels : dict) -> None:
    """Writes the kmer models to a file.

    Args:
        filepath (str):
        kmerModels (dict):
    """
    with open(filepath, 'w') as w:
        w.write('kmer\tlevel_mean\tlevel_stdv\n')
        for kmer in kmerModels:
            mean = kmerModels[kmer][0]
            stdev = kmerModels[kmer][1]
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

def openCPPScript(script : str) -> Popen:
    '''
    Opens a subprocess of the given script.

    Parameters
    ----------
    script : str
        Path of script

    Returns
    -------
    subprocess : Popen
    '''
    # print("Popen call:", ' '.join(script))
    return Popen(script, stdout=PIPE, stdin=PIPE, stderr=PIPE, text=True)

def openCPPScriptATrain(script : str, params : dict) -> Popen:
    '''
    Open cpp script with transition params in transition-training mode

    Parameters
    ----------
    script : str
        Path of script
    params : dict
        {str : float}

    Returns
    -------
    subprocess : Popen
    '''
    script = [script]
    for param in params:
        script.extend([f"-{param}", str(params[param])])
    script.append("--atrain")
    return openCPPScript(script)

def openCPPScriptTrain(script : str, params : dict, model_file : str) -> Popen:
    '''
    Open cpp script with transition params in training mode

    Parameters
    ----------
    script : str
        Path of script
    params : dict
        {str : float}

    Returns
    -------
    subprocess : Popen
    '''
    script = [script]
    for param in params:
        script.extend([f"-{param}", str(params[param])])
    script.extend(["--train", "--model", model_file])
    return openCPPScript(script)

def openCPPScriptParams(script : str, params : dict) -> Popen:
    '''
    Open cpp script with transition params

    Parameters
    ----------
    script : str
        Path of script
    params : dict
        {str : float}

    Returns
    -------
    subprocess : Popen
    '''
    script = [script]
    for param in params:
        if param == 'p':
            script.extend(['-p'])
        else:
            script.extend([f"-{param}", str(params[param])])
    return openCPPScript(script)

def openCPPScriptCalcZ(script : str, params : dict, model : str = None) -> Popen:
    '''
    Open cpp script with transition params in calculate-Z mode

    Parameters
    ----------
    scrip : str
        Path of script
    params : dict
        {str : float}

    Returns
    -------
    subprocess : Popen
    '''
    script = [script]
    for param in params:
        script.extend([f"-{param}", str(params[param])])
    script.append("--calcZ")
    if model is not None:
        script.extend(["--model", model])
    return openCPPScript(script)

def calcZ(signal : np.ndarray, read : str, params : dict, script : str, model : str = None, readid : str = None) -> float:
    pipe = openCPPScriptCalcZ(script, params, model)
    result, errors = feedPipe(signal, read, pipe)
    if errors:
        print(errors)
        return readid
    Z = float(result)
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
    signal = ",".join([f'{x:.5f}' for x in signal])
    cookie = f"{signal}\n{read}\n"
    # with open("cookie.txt", 'w') as w:
    #     w.write(cookie)
    # print("Written Cookie")
    results, errors = pipe.communicate(input=cookie)
    # print(result)
    return results.strip(), errors.strip()

# https://stackoverflow.com/questions/32570029/input-to-c-executable-python-subprocess
def trainTransitionsEmissions(signal : np.ndarray, read : str, params : dict, script : str, model : str, readid : str) -> tuple|str:
    '''
    Parse & feed signal & read to the C++ segmentation script.

    Parameters
    ----------
    signal : np.ndarray
        orientation must match with read
    read : str
        orientation must match with signal
    stream
        Open stdin stream of the C++ segmentation algorithm

    Returns
    -------
    segments : np.ndarray

    params : dict
        {str : float}
    Z : float
    '''
    pipe = openCPPScriptTrain(script, params, model)
    
    # print(len(signal), len(read))

    result, errors = feedPipe(signal, read, pipe)
    if errors:
        print(errors)
        print(f"ERROR while extracting result in {readid}")
        # with open("failed_input.txt", 'w') as w:
        #     w.write(str(signal.tolist()).replace(' ', '').replace('[', '').replace(']', '') + '\n' + read + '\n')
        # raise SegmentationError(readid)
        return readid
    transitionParams, modelParams, Z = result.split('\n')

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
def feedSegmentation(signal : np.ndarray, read : str, script : str, sigOffset : int, params : dict = None) -> tuple:
    '''
    Parse & feed signal & read to the C++ segmentation script.
    Opens and closes a pipe to the given script.

    Parameters
    ----------
    signal : np.ndarray
        orientation must match with read
    read : str
        orientation must match with signal
    script : str
        script file name
    signal_offset : int
    params : dict
        params for the c++ script in short form: e.g. {-m : "model"}

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
    result, errors = feedPipe(signal, read, pipe)
    
    if errors:
        print(f"ERROR appeared for {read}")
        print(signal)
        print(read)
        print(errors)
        with open("failed_input.txt", "w") as w:
            w.write(str(np.around(signal, 4).tolist()).replace(' ', '').replace('[', '').replace(']', ''))
            w.write('\n')
            w.write(read)
            w.write('\n')
        exit(1)

    try:
        output, probs = result.split('\n')
        probs = np.array(list(map(float, probs.split(',')[:-1])))[1:]
    except:
        try:
            output  = result.split('\n')
            probs = None
        except:
            print("ERROR while extracting result in {read}")
            print(signal)
            print(read)
            with open("failed_input.txt", "w") as w:
                w.write(str(np.around(signal, 4).tolist()).replace(' ', '').replace('[', '').replace(']', ''))
                w.write('\n')
                w.write(read)
                w.write('\n')
            exit(1)

    # receive segmentation result and format output into np.ndarray
    segments = formatSegmentationOutput(output, sigOffset, len(signal) + sigOffset, read[::-1])
    return segments, probs

def feedSegmentationAsynchronous(CPP_SCRIPT : str, params : dict, signal : np.ndarray, read : str, signal_offset : int, readid : str, signalid : str, queue : Queue) -> None:
    '''
    Parse & feed signal & read to the C++ segmentation script.
    Needs an open pipe for communication.

    Parameters
    ----------
    signal : np.ndarray
        in 3 -> 5' direction
    read : str
        in 3' -> 5' direction
    signal_offset : int
        offset index of the signal
    readid : str
        id within the basecall .bam file of dorado
    signalid : str
        "pi" flag within the basecall .bam file of dorado
    pipe : Popen
    queue : Queue
    '''
    pipe = openCPPScriptParams(CPP_SCRIPT, params)
    output, errors = feedPipe(signal, read, pipe)
    if errors:
        queue.put("error: " + errors)
        return
    segments = formatSegmentationOutput(output, signal_offset, len(signal) + signal_offset, read[::-1])
    out = formatSegmentation(readid, signalid, segments)
    if out:
        queue.put(out)

def formatSegmentationOutput(output : str, sigOffset : int, lastIndex : int, read : str) -> np.ndarray:
    '''
    Receives the segmentation output from CPP script and returns it as a np.ndarray

    Parameters
    ----------
    output : str
        Segmentation output from CPP script
    signal_offset : int
        starting index offset
    lastIndex : int
        ending index with offset already added
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
    output = output.split(';')[:-1]
    for i, segment in enumerate(output):
        # split segment state
        state = segment[0]
        segment = segment[1:]

        # get basepos, signal start and polish if existing
        try:
            basepos, start, prob, polish = segment.split(',')
        except ValueError:
            basepos, start, prob = segment.split(',')
            polish = 'NA'

        start = int(start) + sigOffset
        end = int(output[i+1][1:].split(',')[1]) + sigOffset if i < len(output)-1 else lastIndex

        # convert basepos to 5' -> 3' direction
        pos = len(read) - int(basepos) - 1
        segments.append([int(start), int(end), pos, read[pos], read[max(0, pos-2):min(len(read), pos+3)], state, float(prob), polish])

    return np.array(segments, dtype=object)

def formatSegmentation(readid : str, signalid : str, segmentation : np.ndarray) -> str:
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
    lines = []
    for segment in segmentation:
        line = f"{readid},{signalid}," + ','.join(map(str, segment))
        lines.append(line)
    return '\n'.join(lines) + '\n'

def stopFeeding(pipe : Popen) -> None:
    '''
    Stop feeding a pipe with data

    Parameters
    ----------
    pipe : Popen
        the pipe to stop feeding
    '''
    pipe.kill()
    return None

def plotParameters(param_file : str, outdir : str) -> None:
    '''
    Generate line plots for each parameter in the given CSV file over training batches, 
    saving them as PDF files in the specified output directory.

    Parameters
    ----------
    param_file : str
        Path to the CSV file containing parameter data with columns for 'epoch', 'batch', and parameters.
    outdir : str
        Directory to save the output PDF files.
    '''
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
