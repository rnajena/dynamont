#!/usr/bin/env python
# author: Jannes Spangenberg
# github: https://github.com/JannesSP
# website: https://jannessp.github.io

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns
import sys
from os.path import join, dirname
from dynamont import Aligner

def hampel(signal : np.ndarray, WINDOW : int = 3, n_sigmas : float = 3.0) -> None:
    """
    Apply Hampel filter in place to detect and replace outliers in a signal.

    Parameters:
        signal (np.ndarray): The input signal (1D array).
        WINDOW (int): The size of the sliding window. Defaults to 3.
        n_sigmas (float): The threshold in terms of standard deviations (MAD).
                          Defaults to 3, which is commonly used.
    """
    if signal.size <= WINDOW:
        return

    original = signal.copy()
    windows = np.lib.stride_tricks.sliding_window_view(original, WINDOW)

    n_values = len(signal) - WINDOW - (WINDOW % 2 == 0)
    windows = windows[:n_values]

    median = np.median(windows, axis=1)
    mad = np.median(np.abs(windows - median[:, None]), axis=1)
    sigma = 1.4826 * mad

    center = original[WINDOW // 2:len(signal) - WINDOW // 2 - 1]
    replace = np.abs(center - median) > n_sigmas * sigma
    signal[WINDOW // 2:len(signal) - WINDOW // 2 - 1][replace] = median[replace]

    return

def cnt_nts(sequence):
    """
    Count the occurrences of A, C, G, and T in a DNA sequence.

    Parameters:
        sequence (str): DNA sequence as a string.

    Returns:
        dict: Dictionary with counts of A, C, G, and T.
    """
    return {
        "A": sequence.count("A"),
        "C": sequence.count("C"),
        "G": sequence.count("G"),
        "T": sequence.count("T"),
    }

def cnt_nts_ratios(sequence):
    """
    Count the occurrences of A, C, G, and T in a DNA sequence.

    Parameters:
        sequence (str): DNA sequence as a string.

    Returns:
        dict: Dictionary with counts of A, C, G, and T.
    """
    return {
        "A": sequence.count("A") / len(sequence),
        "C": sequence.count("C") / len(sequence),
        "G": sequence.count("G") / len(sequence),
        "T": sequence.count("T") / len(sequence),
    }

class SegmentationError(Exception):
    """Raised when no segmentation was calculated for a read"""
    def __init__(self, read) -> None:
        self.read = read
        self.message = f"No segmentation calculated for {read}"
        super().__init__(self.message)

def _decode_native_state(state) -> str:
    if isinstance(state, bytes):
        return state.decode("ascii")
    if isinstance(state, str):
        return state
    if isinstance(state, (int, np.integer)):
        return chr(int(state))
    return str(state)


def get_model(pore : str) -> str:
    """
    Return the path to the kmer model file for a given pore type.

    Parameters
    ----------
    pore : str
        Pore generation used to sequence the data.

    Returns
    -------
    str
        Path to the kmer model file.
    """
    MODELS = {
        "rna002" : "models/rna/rna002/rna002_5mer.model",
        "rna004" : "models/rna/rna004/rna004_9mer.model",
        "dna_r10_260bps" : "models/dna/r10.4.1/dna_r10.4.1_e8.2_260bps.model",
        "dna_r10_400bps" : "models/dna/r10.4.1/dna_r10.4.1_e8.2_400bps.model",
    }
    base_dir = dirname(dirname(dirname(__file__)))
    return join(base_dir, MODELS.get(pore, pore))

def read_kmer_model(file : str) -> dict:
    """
    Reads the kmer models from a file.
    Parameters
    ----------
    file : str
        path to the kmer model file

    Returns
    -------
    models : dict
        format: {kmer : (mean, stdev)}
    """
    df = pd.read_csv(file, sep='\t')
    models = pd.Series(zip(df.level_mean.values, df.level_stdv.values), index=df.kmer).to_dict()
    return models

def write_kmer_model(file : str, kmer_model : dict) -> None:
    """
    Writes the kmer models to a file.

    Parameters
    ----------
    file (str):
        Path to the output file.
    kmer_model (dict):
        Dictionary containing kmer models.
    """
    with open(file, 'w') as w:
        w.write('kmer\tlevel_mean\tlevel_stdv\n')
        for kmer in kmer_model:
            mean = kmer_model[kmer][0]
            stdev = kmer_model[kmer][1]
            w.write(f'{kmer}\t{mean}\t{stdev}\n')

def _make_native_aligner(model: str, params: dict, mode: str) -> Aligner:
    pore = params.get("r") or params.get("pore")
    if pore is None:
        raise ValueError("Missing pore type in segmentation parameters")

    threads = int(params.get("t", 1))
    band = int(params.get("band", 400))
    return Aligner(model, pore, mode=mode, threads=threads, band=band)

def train_transition_emission(signal : np.ndarray, read : str, params : dict, script : str, model : str, signalid : str):
    try:
        aligner = _make_native_aligner(model, params, script)
        transitionParams = {key: float(value) for key, value in params.items() if key not in {"r", "t", "band"}}

        if script == "basic":
            result = aligner.train(signal, read)
            trainedParams = {key: float(value) for key, value in result["transition_params"].items()}
            modelKeys = list(read_kmer_model(model).keys())
            newModels = {
                kmer: (float(entry["mean"]), float(entry["stdev"]))
                for kmer, entry in zip(modelKeys, result["emission_model"])
            }
            return trainedParams, newModels, float(result["Z"])

        result = aligner.align(signal, read, calc_probabilities=False)
        return transitionParams, read_kmer_model(model), float(result["Z"])
    except Exception as error:
        print(f"error: native, {error} T: {len(signal)} N: {len(read)} Sid: {signalid}", file=sys.stderr)
        return signalid

def calcZ(signal : np.ndarray, read : str, params : dict, script : str, model : str, signalid : str):
    try:
        aligner = _make_native_aligner(model, params, script)
        result = aligner.align(signal, read, calc_probabilities=False)
        return float(result["Z"])
    except Exception as error:
        print(f"error: native, {error} T: {len(signal)} N: {len(read)} Sid: {signalid}", file=sys.stderr)
        return signalid

def segmentation_to_string(
    result: dict,
    readid: str,
    signalid: str,
    sigOffset: int,
    lastIndex: int,
    read: str,
    kmerSize: int,
    rna: bool,
) -> bytes:
    seq_pos = result["sequence_positions"]
    sig_pos = result["signal_positions"]
    probs = result["probabilities"]
    states = result["states"]
    polishes = result.get("polishes")
    n_segments = len(seq_pos)
    lines = []

    for i in range(n_segments):
        basepos = int(seq_pos[i])
        start = int(sig_pos[i]) + sigOffset
        end = int(sig_pos[i + 1]) + sigOffset if i < n_segments - 1 else lastIndex
        state = _decode_native_state(states[i])
        polish = str(polishes[i]) if polishes is not None and polishes[i] else "NA"
        motif = read[
            max(0, basepos - (kmerSize // 2)):
            min(len(read), basepos + (kmerSize // 2) + 1)
        ]
        base = read[basepos]

        if rna:
            motif = motif[::-1]
            basepos = len(read) - basepos - 1

        lines.append(
            f"{readid},{signalid},{start},{end},{basepos},{base},{motif},"
            f"{state},{float(probs[i]):.6f},{polish}\n"
        )

    return "".join(lines).encode("utf-8")

def plt_parameters(file : str, outdir : str) -> None:
    '''
    Generate line plots for each parameter in the given CSV file over training batches, 
    saving them as PDF files in the specified output directory.

    Parameters
    ----------
    file : str
        Path to the CSV file containing parameter data with columns for 'epoch', 'batch', and parameters.
    outdir : str
        Directory to save the output PDF files.
    '''
    df = pd.read_csv(file, sep=',')
    for column in df:
        if column in ['epoch', 'batch', 'read']:
            continue
        sns.set_theme()        
        sns.lineplot(data=df, x="batch", y=column, hue='epoch')
        plt.title(f"{column} parameter change during training")
        plt.ylabel("Parameter Value")
        print("Savefig: ", join(outdir, f"{column}.pdf"), file=sys.stderr)
        plt.savefig(join(outdir, f"{column}.pdf"))
        plt.close()