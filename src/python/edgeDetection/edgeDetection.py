#!/usr/bin/env python
# author: Jannes Spangenberg
# e-mail: jannes.spangenberg@uni-jena.de
# github: https://github.com/JannesSP
# website: https://jannessp.github.io

import h5py
import numpy as np
import read5_ont
import pysam
import pywt
import multiprocessing as mp
from argparse import ArgumentDefaultsHelpFormatter, ArgumentParser, Namespace
from os.path import join
from src.python.segmentation.FileIO import hampelFilter

def parse() -> Namespace:
    parser = ArgumentParser(
        formatter_class=ArgumentDefaultsHelpFormatter
    )
    parser.add_argument("-r", "--raw", type=str, required=True, metavar="FAST5|POD5", help="Input raw file format (FAST5, POD5)")
    parser.add_argument("-b", "--basecalls", type=str, required=True, metavar="BAM", help="Basecalls of ONT training data as.bam file")
    parser.add_argument("-o", "--output", type=str, required=True, metavar="HDF5", help="Output HDF5 file")
    parser.add_argument("--pore", type=str, required=True, metavar="STR", choices=["rna_r9", "dna_r9", "rna_rp4", "dna_r10_260bps", "dna_r10_400bps"])
    parser.add_argument("-p", "--processes", type=int, default=mp.cpu_count(), metavar="INT", help="Number of processes to use for parallel processing (default: all available CPUs)")
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
    from scipy.signal import find_peaks
    halfThreshold = threshold / 2
    # Perform continuous wavelet transform on the signal with scale 6
    coef, _ = pywt.cwt(signal, [6], wavelet)
    coef = np.abs(coef[0])

    # Use SciPy's find_peaks to locate peaks above the threshold
    peaks, _ = find_peaks(coef, height=threshold, distance=3)

    # Further filter the peaks based on a local standard deviation threshold
    final_peaks = []
    for peak in peaks:
        # Calculate standard deviation in a 6-element window around each peak
        window = coef[max(0, peak-3):peak+3]
        if np.std(window) < halfThreshold:
            final_peaks.append(peak)

    return np.array(final_peaks)

def writer(h5file : str, q : mp.Queue) -> None:
    with h5py.File(h5file, 'w') as hdf:
        i = 0
        while True:
            signalid, waveletEdges = q.get()

            if signalid == 'kill':
                break

            i+=1

            print(f"Extracting edges for {i} reads", end='\r')
        
            key = f'{signalid}/waveletEdge'
            # save to hdf5 file
            if not signalid in hdf.keys():
                hdf.create_dataset(key, data=waveletEdges, maxshape=(None,), dtype = 'f')
            else:
                current_shape = hdf[key].shape[0]
                new_shape = current_shape + len(waveletEdges)
                hdf[key].resize(new_shape, axis=0)
                hdf[key][current_shape:new_shape] = waveletEdges

        print()
        
def extractingEdges(signalid : str, rawFile : str, start : int, end : int, threshold : float, shift : float, scale : float, pore : str, queue : mp.Queue) -> None:
    r5 = read5_ont.read(rawFile)
    if pore in ["rna_r9", "dna_r9"]:
        signal = r5.getpASignal(signalid)
    else:
        signal = r5.getSignal(signalid)
    signal = (signal - shift) / scale
    hampelFilter(signal, 6, 5.)
    waveletEdges = waveletPeaks(signal[start:end], 'gaus1', threshold) + start
    queue.put((signalid, waveletEdges))
        
def wavelet(raw : str, basecalls : str, outfile: str, processes : int, pore : str) -> None:
    """
    Processes raw signal data and basecalls to detect edges using the wavelet transform and stores the results in an HDF5 file.

    Parameters
    ---
    - raw (str): Path to the directory containing raw signal files.
    - basecalls (str): Path to the basecalls file in BAM format.
    - outfile (str): Path to the output HDF5 file where the detected wavelet edges will be stored.

    Returns
    ---
    None: This function does not return any value. It writes the detected wavelet edges to the specified HDF5 file.
    """
    # processes = 40
    pool = mp.Pool(processes)
    queue = mp.Manager().Queue()
    watcher = pool.apply_async(writer, (outfile, queue))
    jobs = [None for _ in range(processes - 1)]
    
    with pysam.AlignmentFile(basecalls, "r" if basecalls.endswith('.sam') else "rb", check_sq=False) as samfile:
        for i, basecalled_read in enumerate(samfile.fetch(until_eof=True)):
            readid = basecalled_read.query_name
            signalid = basecalled_read.get_tag("pi") if basecalled_read.has_tag("pi") else readid
            sp = basecalled_read.get_tag("sp") if basecalled_read.has_tag("sp") else 0 # if split read get start offset of the signal
            ns = basecalled_read.get_tag("ns") # ns:i: 	the number of samples in the signal prior to trimming
            ts = basecalled_read.get_tag("ts") # ts:i: 	the number of samples trimmed from the start of the signal
            rawFile = join(raw, basecalled_read.get_tag("fn"))
            shift = basecalled_read.get_tag("sm")
            scale = basecalled_read.get_tag("sd")
            # signal = (signal - shift) / scale
    
            jobs[i%(processes-1)] = pool.apply_async(extractingEdges, (signalid, rawFile, sp+ts, sp+ns, 1.1, shift, scale, pore, queue))
        
    for job in jobs:
        job.get()
    
    queue.put(('kill', None))
    watcher.get()
        
    pool.close()
    pool.join()

    print(f'Done extracting edges for {i} reads')

def main() -> None:
    args = parse()
    print('Start extracting')
    wavelet(args.raw, args.basecalls, args.output, args.processes)

if __name__ == '__main__':
    main()