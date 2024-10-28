#!/usr/bin/env python
# author: Jannes Spangenberg
# e-mail: jannes.spangenberg@uni-jena.de
# github: https://github.com/JannesSP
# website: https://jannessp.github.io

from argparse import ArgumentDefaultsHelpFormatter, ArgumentParser, Namespace
import h5py
import numpy as np
import read5
from hampel import hampel
import pysam
import pywt
from os.path import join

def parse() -> Namespace:
    parser = ArgumentParser(
        formatter_class=ArgumentDefaultsHelpFormatter
    )
    parser.add_argument("-r", "--raw", type=str, required=True, metavar="FAST5|SLOW5|POD5", help="Input raw file format (FAST5, SLOW5, POD5)")
    parser.add_argument("-b", "--basecalls", type=str, required=True, metavar="BAM", help="Basecalls of ONT training data as.bam file")
    parser.add_argument("-o", "--output", type=str, required=True, metavar="HDF5", help="Output HDF5 file")
    return parser.parse_args()

def waveletPeaks(signal: np.ndarray, wavelet: str = 'gaus1', threshold: float = 1.1) -> np.ndarray:
    """
    Detects peaks in a signal using continuous wavelet transform.

    Parameters
    ---
    - signal (np.ndarray): The input signal as a numpy array.
    - wavelet (str): The type of wavelet to use for the transform. Default is 'gaus1'.
    - threshold (float): The threshold value for peak detection. Peaks with values above this threshold are considered.

    Returns
    ---
    np.ndarray: An array of indices where peaks are detected in the input signal.
    """
    # Perform continuous wavelet transform on the signal with scale 6
    coef, freqs = pywt.cwt(signal, [6], wavelet)
    c = np.abs(coef[0])
    peaks = []
    # Iterate over the wavelet coefficients to find local maxima
    for i in range(len(c)):
        # Check if the current point is a local maximum and the standard deviation around it is low
        if c[i] >= threshold and c[i] == np.max(c[max(0, i-3):i+3]) and np.std(c[max(0, i-3):i+3]) < threshold/2:
            peaks.append(i)
    peaks = np.array(peaks)
    return peaks

def wavelet(raw : str, basecalls : str, output_file: str) -> None:
    """
    Processes raw signal data and basecalls to detect edges using the wavelet transform and stores the results in an HDF5 file.

    Parameters
    ---
    - raw (str): Path to the directory containing raw signal files.
    - basecalls (str): Path to the basecalls file in BAM format.
    - output_file (str): Path to the output HDF5 file where the detected wavelet edges will be stored.

    Returns
    ---
    None: This function does not return any value. It writes the detected wavelet edges to the specified HDF5 file.
    """
    with h5py.File(output_file, 'w') as hdf:
        with pysam.AlignmentFile(basecalls, "r" if basecalls.endswith('.sam') else "rb", check_sq=False) as samfile:

            for basecalled_read in samfile.fetch(until_eof=True):
                readid = basecalled_read.query_name
                signalid = basecalled_read.get_tag("pi") if basecalled_read.has_tag("pi") else readid
                sp = basecalled_read.get_tag("sp") if basecalled_read.has_tag("sp") else 0 # start sample of split read by the basecaller
                ns = basecalled_read.get_tag("ns") # numbers of samples used in basecalling
                ts = basecalled_read.get_tag("ts") # start sample of basecalling
                rawFile = join(raw, basecalled_read.get_tag("fn"))
                r5 = read5.read(rawFile)

                shift = basecalled_read.get_tag("sm")
                scale = basecalled_read.get_tag("sd")
                signal = r5.getpASignal(signalid).astype(np.float32)
                signal = (signal - shift) / scale
                filtered = hampel(signal, 6, 5.).filtered_data

                waveletEdges = waveletPeaks(filtered[sp+ts:sp+ns]) + (ts + sp)

                key = f'{signalid}/waveletEdge'
                # save to hdf5 file
                if not signalid in hdf.keys():
                    hdf.create_dataset(key, data=waveletEdges, maxshape=(None,))
                else:
                    hdf[key].resize(hdf[key].shape[0] + len(waveletEdges))
                    hdf[key][hdf[key].shape[0]:] = waveletEdges

def main() -> None:
    args = parse()
    wavelet(args.raw, args.basecalls, args.output)

if __name__ == '__main__':
    main()