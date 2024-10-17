#!/usr/bin/env python
# author: Jannes Spangenberg
# e-mail: jannes.spangenberg@uni-jena.de
# github: https://github.com/JannesSP
# website: https://jannessp.github.io

from argparse import ArgumentDefaultsHelpFormatter, ArgumentParser, Namespace
import numpy as np
from os.path import exists, join, dirname
from os import makedirs, name
from FileIO import calcZ, plotParameters, trainTransitionsEmissions, readKmerModels, writeKmerModels
import read5
import multiprocessing as mp
from hampel import hampel
from datetime import datetime
from collections import deque
import pysam

class ManagedList:
    def __init__(self, values, max_size=100, oldest_slice_size=10, outlierPercentile = 10):
        # self.type = type # Emission or Transition
        self.max_size : int = max_size
        self.oldest_slice_size : int = oldest_slice_size  # The size of the slice of the oldest values to check for outliers
        self.values = deque(values, maxlen=max_size)  # Deque to efficiently store values
        self.outlierPercentile = outlierPercentile

    def _iqr(self):
        """Helper method to calculate IQR (Interquartile Range) and determine outliers."""
        if len(self.values) < 4:
            return None, None, None  # Not enough data for IQR
        q1 = np.percentile(self.values, self.outlierPercentile)
        q3 = np.percentile(self.values, 100-self.outlierPercentile)
        iqr = q3 - q1
        return q1, q3, iqr

    def _is_outlier(self, value, q1, q3, iqr):
        """Determines if a value is an outlier based on the current IQR."""
        lower_bound = q1 - 1.5 * iqr
        upper_bound = q3 + 1.5 * iqr
        return value < lower_bound or value > upper_bound

    def _find_strongest_outlier_in_slice(self):
        """Finds the strongest outlier in the oldest slice of values with respect to the entire list."""
        
        # Calculate IQR for the entire list
        q1, q3, iqr = self._iqr()

        # Take a slice of the oldest values
        oldest_slice = list(self.values)[:self.oldest_slice_size]

        # Find outliers in the oldest slice
        outliers = [(value, i) for i, value in enumerate(oldest_slice) if self._is_outlier(value, q1, q3, iqr)]
        
        if not outliers:
            return None  # No outliers in the slice
        
        # Find the strongest outlier (farthest from the nearest bound)
        def outlier_strength(value):
            return max(abs(value - (q1 - 1.5 * iqr)), abs(value - (q3 + 1.5 * iqr)))

        # Sort by strength of the outlier and return the index of the strongest outlier
        strongest_outlier = max(outliers, key=lambda x: outlier_strength(x[0]))
        return strongest_outlier[1]  # Return the index within the slice

    def _remove_oldest_or_outlier(self):
        """Removes the strongest outlier from the oldest slice, or the oldest value if no outliers exist."""
        outlier_index = self._find_strongest_outlier_in_slice()
        
        if outlier_index is not None:
            # Remove the strongest outlier from the oldest slice
            del self.values[outlier_index]
        else:
            # If no outliers, remove the oldest value
            self.values.popleft()

    def mean(self):
        """Returns the mean of the values in the managed list."""
        if len(self.values) == 0:
            return None  # Handle the case where the list is empty
        return np.mean(self.values)

    def median(self):
        """Returns the median of the values in the managed list."""
        if len(self.values) == 0:
            return None  # Handle the case where the list is empty
        return np.median(self.values)

    def add(self, value):
        """Add a new value to the list, managing overflow and outliers."""
        # If the list is full
        if len(self.values) >= self.max_size:
            # First, remove the strongest outlier from the oldest slice, or just the oldest
            # if self.type == "Emission":
            self._remove_oldest_or_outlier()

        # Add the new value to the list
        self.values.append(value)

    def get_list(self):
        """Returns the current managed list."""
        return list(self.values)
    
    def clear(self):
        """Removes all elements of the managed list."""
        self.values = deque()


# class ManagedList:
#     def __init__(self, values, max_size=100):
#         """Initialize the ManagedList with a maximum size."""
#         self.max_size = max_size
#         self.values = deque(values, maxlen=max_size)  # deque automatically limits size

#     def add(self, value):
#         """Add a new value to the list, removing the oldest if necessary."""
#         self.values.append(value)

#     def get_list(self):
#         """Returns the current list of values."""
#         return list(self.values)

#     def __repr__(self):
#         """String representation of the ManagedList."""
#         return f"ManagedList({list(self.values)})"
    
#     def mean(self):
#         """Returns the mean of the values in the managed list."""
#         if len(self.values) == 0:
#             return None  # Handle the case where the list is empty
#         return np.mean(self.values)

#     def median(self):
#         """Returns the median of the values in the managed list."""
#         if len(self.values) == 0:
#             return None  # Handle the case where the list is empty
#         return np.median(self.values)

def parse() -> Namespace:
    parser = ArgumentParser(
        formatter_class=ArgumentDefaultsHelpFormatter
    )
    # required
    parser.add_argument('-r', '--raw',   type=str, required=True, metavar="PATH", help='Path to raw ONT data. [POD5|FAST5|SLOW5]')
    parser.add_argument('-b', '--basecalls', type=str, required=True, metavar="BAM", help='Basecalls of ONT training data as .bam file')
    parser.add_argument('-o', '--outdir',   type=str, required=True, metavar="PATH", help='Outpath to write files')
    parser.add_argument('-p', '--pore',  type=str, required=True, choices=["rna_r9", "dna_r9", "rna_rp4", "dna_r10_260bps", "dna_r10_400bps"], help='Pore generation used to sequence the data')
    parser.add_argument('--mode',  type=str, required=True, choices=['basic', 'banded', 'resquiggle'], help='Segmentation algorithm used for segmentation')
    # optional
    parser.add_argument('--model_path', type=str, help='Which initial kmer models to use for training')
    parser.add_argument('--batch_size', type=int, default=24, help='Number of reads to train before updating')
    parser.add_argument('--max_batches', type=int, default=None, help='Numbers of batches to train each epoch')
    parser.add_argument('-e', '--epochs', type=int, default=1, help='Number of training epochs')
    parser.add_argument('-c', '--minSegLen', type=int, default=1, help='Minmal allowed segment length')
    parser.add_argument('-q', '--qscore', type=int, default=0, help='Minmal allowed quality score')
    return parser.parse_args()

def train(dataPath : str, basecalls : str, batch_size : int, epochs :int, param_file : str, mode : str, model_path : str, minSegLen : int, max_batches : int, pore : str, minQual : float = None) -> None:

    kmerModels = readKmerModels(model_path)
    trainedModels = join(dirname(param_file), 'trained_0_0.model')
    writeKmerModels(trainedModels, kmerModels)
    
    param_writer = open(param_file, 'w')

    if mode == 'basic':
        mode = 'dynamont_NT'
        transitionParams = {
            'e1': 1.0,
            'm1': 0.03,
            'e2': 0.97
            }
    elif mode == 'banded':
        mode = 'dynamont_NT_banded'
        transitionParams = {
            'e1': 1.0,
            'm1': 0.03,
            'e2': 0.97
            }
    elif mode == 'resquiggle':
        mode = 'dynamont_NTK'
        transitionParams = {
            'a1': 0.030570760232128354,
            'a2': 0.6101330423858251,
            'p1': 0.0383208584222504,
            'p2': 0.08685413259505585,
            'p3': 0.015545502978626788,
            's1': 0.0010061921517067177,
            's2': 0.001534865889565432,
            's3': 0.13149698483529376,
            'e1': 1.0,
            'e2': 0.9989938127816652,
            'e3': 0.9616791335752949,
            'e4': 0.8801771554834845,
            'i1': 0.0008630902783022252,
            'i2': 0.24282445342753348
            }
    else:
        print(f'Mode {mode} not implemented')
        exit(1)

    CPP_SCRIPT = join(dirname(__file__), f'{mode}')
    if name == 'nt': # check for windows
        CPP_SCRIPT+='.exe'

    # collect trained parameters to get an ensemble training in the end to reduce outlier trainings
    paramCollector = {kmer:(ManagedList([kmerModels[kmer][0]]), ManagedList([kmerModels[kmer][1]])) for kmer in kmerModels}
    paramCollector = paramCollector | {param : ManagedList([transitionParams[param]]) for param in transitionParams}
    kmerSeen = set()
    param_writer.write("epoch,batch,read,")
    for param in transitionParams:
        param_writer.write(param+',')
    param_writer.write("Zchange\n")
    i = 0
    qualSkipped = 0
    noMatchingReadid = 0

    with mp.Pool(batch_size) as p:

        for e in range(epochs):
            mp_items = []
            training_readids = []
            batch_num = 0
            oldFile = None

            with pysam.AlignmentFile(basecalls, "r" if basecalls.endswith('.sam') else "rb", check_sq=False) as samfile:
                for basecalled_read in samfile.fetch(until_eof=True):
                    
                    # skip low qual reads if activated
                    qual = basecalled_read.get_tag("qs")
                    if minQual and qual < minQual:
                        qualSkipped+=1
                        continue

                    # init read
                    seq = basecalled_read.query_sequence
                    readid = basecalled_read.query_name
                    ts = basecalled_read.get_tag("ts")
                    rawFile = join(dataPath, basecalled_read.get_tag("fn"))

                    if oldFile != rawFile:
                        oldFile = rawFile
                        r5 = read5.read(rawFile)

                    # fill batch
                    if len(mp_items) < batch_size:
                        # saw more consistency for short reads when using the mean
                        try:
                            signal = r5.getZNormSignal(readid, "mean")[ts:]
                        except:
                            noMatchingReadid+=1
                        signal = hampel(signal, 6, 5.).filtered_data # small window and high variance allowed: just to filter outliers that result from sensor errors, rest of the original signal should be kept
                        mp_items.append([signal, seq[::-1], transitionParams | {'r' : pore}, CPP_SCRIPT, trainedModels, minSegLen, readid])
                        training_readids.append(readid)

                    if len(mp_items) == batch_size:
                        print("============================")
                        print(f"{datetime.now().strftime('%Y-%m-%d_%H-%M-%S')}: Training epoch: {e}, reads: {i}, batch: {batch_num}\n{transitionParams}")
                        print("Training with read:", training_readids)
                        batch_num += 1
                        Zs = []

                        for readid, result in enumerate(p.starmap(trainTransitionsEmissions, mp_items)):
                            if isinstance(result, str):
                                print(f"No segmentation calculated for {result} in {e}: {trainedModels}.")
                                del mp_items[readid]
                                continue

                            trainedParams, newModels, Z = result
                            i += 1
                            Zs.append(Z)

                            for param in trainedParams:
                                paramCollector[param].add(trainedParams[param])

                            for kmer in newModels:
                                kmerSeen.add(kmer)
                                paramCollector[kmer][0].add(newModels[kmer][0])
                                paramCollector[kmer][1].add(newModels[kmer][1])

                        print(f"Zs: {Zs}")

                        # update parameters
                        param_writer.write(f'{e},{batch_num},{i},') # log
                        for param in transitionParams:
                            try:
                                transitionParams[param] = paramCollector[param].mean()
                            except:
                                print(param, paramCollector[param].get_list())
                                exit(1)
                            param_writer.write(f'{transitionParams[param]},') # log

                        for kmer in kmerSeen:
                            kmerModels[kmer] = [paramCollector[kmer][0].mean(), paramCollector[kmer][1].mean()]

                        trainedModels = join(dirname(trainedModels), f"trained_{e}_{batch_num}.model")
                        writeKmerModels(trainedModels, kmerModels)
                        param_writer.flush()

                        # rerun with new parameters to compare Zs
                        for j in range(len(mp_items)):
                            mp_items[j][2] = transitionParams
                            mp_items[j][4] = trainedModels
                        Zdiffs = []
                        for j, result in enumerate(p.starmap(calcZ, mp_items)):
                            if isinstance(result, str):
                                print(f"No segmentation calculated for {result} in {e} calcZ.")
                                Zdiffs.append(0)
                                continue
                            Z = result
                            Zdiffs.append(Z - Zs[j])

                        print(f"Z changes: {Zdiffs}")
                        deltaZ = np.mean(Zdiffs)
                        param_writer.write(f'{deltaZ}\n') # log
                        param_writer.flush() # log
                        # initialize new batch
                        kmerSeen = set()
                        mp_items = []
                        training_readids = []

                        if max_batches is not None and batch_num >= max_batches:
                            break

        param_writer.close()
        print("Done training")
        print(f"Skipped reads due to low quality: {qualSkipped}")

def main() -> None:
    args = parse()

    # due to algorithm min len is 1 by default, minSegLen stores number of positions to add to that length
    assert args.minSegLen>=1, "Please choose a minimal segment length greater than 0"

    # init outdir
    outdir = args.outdir + f'_{datetime.now().strftime("%Y-%m-%d_%H-%M-%S")}'
    if not exists(outdir):
        makedirs(outdir)

    param_file = join(outdir, 'params.csv')

    train(args.raw, args.basecalls, args.batch_size, args.epochs, param_file, args.mode, args.model_path, args.minSegLen - 1, args.max_batches, args.pore, args.qscore)
    plotParameters(param_file, outdir)

if __name__ == '__main__':
    main()