#!/usr/bin/env python
# author: Jannes Spangenberg
# e-mail: jannes.spangenberg@uni-jena.de
# github: https://github.com/JannesSP
# website: https://jannessp.github.io

from argparse import ArgumentDefaultsHelpFormatter, ArgumentParser, Namespace
import numpy as np
from os.path import exists, join, dirname
from os import makedirs, name
from FileIO import getFiles, loadFastx, calcZ, readPolyAStartEnd, plotParameters, trainTransitionsEmissions, readKmerModels, writeKmerModels
from read5 import read
import multiprocessing as mp
from hampel import hampel
import random
from datetime import datetime
from collections import deque


class ManagedList:
    def __init__(self, values, max_size=100, oldest_slice_size=10):
        self.max_size : int = max_size
        self.oldest_slice_size : int = oldest_slice_size  # The size of the slice of the oldest values to check for outliers
        self.values = deque(values)  # Deque to efficiently store values

    def _iqr(self):
        """Helper method to calculate IQR (Interquartile Range) and determine outliers."""
        if len(self.values) < 4:
            return None, None, None  # Not enough data for IQR
        q1 = np.percentile(self.values, 25)
        q3 = np.percentile(self.values, 75)
        iqr = q3 - q1
        return q1, q3, iqr

    def _is_outlier(self, value, q1, q3, iqr):
        """Determines if a value is an outlier based on the current IQR."""
        lower_bound = q1 - 1.5 * iqr
        upper_bound = q3 + 1.5 * iqr
        return value < lower_bound or value > upper_bound

    def _find_strongest_outlier_in_slice(self):
        """Finds the strongest outlier in the oldest slice of values with respect to the entire list."""
        if len(self.values) < 4:
            return None  # Not enough values to calculate outliers
        
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
            self._remove_oldest_or_outlier()

        # Add the new value to the list
        self.values.append(value)

    def get_list(self):
        """Returns the current managed list."""
        return list(self.values)
    
    def clear(self):
        """Removes all elements of the managed list."""
        self.values = deque()

def parse() -> Namespace:
    parser = ArgumentParser(
        formatter_class=ArgumentDefaultsHelpFormatter
    )
    parser.add_argument('--raw', type=str, required=True, help='Raw ONT training data')
    parser.add_argument('--fastx', type=str, required=True, help='Basecalls of ONT training data')
    parser.add_argument('--polya', type=str, required=True, help='Poly A table from nanopolish polya containing the transcript starts')
    parser.add_argument('--out', type=str, required=True, help='Outpath to write files')
    parser.add_argument('--mode', type=str, choices=['basic', 'indel', '3d', '3d_reduced'], required=True, help='Segmentation algorithm used for segmentation')
    parser.add_argument('--model_path', type=str, default=join(dirname(__file__), '..', '..', 'data', 'norm_models', 'rna_r9.4_180mv_70bps_extended_stdev0_25.model'), help='Which models to train')
    parser.add_argument('--pore', type=int, choices=[9, 10], default=9, help='Pore generation used to sequence the data')
    parser.add_argument('--batch_size', type=int, default=24, help='Number of reads to train before updating')
    parser.add_argument('--epochs', type=int, default=1, help='Number of training epochs')
    parser.add_argument('--same_batch', action='store_true', help='Always train on first batch of signals in the file')
    parser.add_argument('--random_batch', type=int, metavar='seed', default=None, help='Randomizes the batches of each file with a given seed')
    # parser.add_argument('--stdev', type=float, default=None, help='Set a stdev value instead of using the calculated one for the models')
    parser.add_argument('--minSegLen', type=int, default=1, help='Minmal allowed segment length')
    parser.add_argument('--unnormalised', action="store_true", help="Use polyA standardized signal in instead of median normalised signal. Make sure the model matches this mode!")
    return parser.parse_args()

def train(rawdatapath : str, fastxpath : str, polya : dict, batch_size : int, epochs :int, param_file : str, mode : str, same_batch : bool, random_batch : bool, kmerModels : dict, trainedModels : str, minSegLen : int, errorfile : str, unnormalised : bool) -> None:
    # print("minSegLen", minSegLen)
    files = getFiles(rawdatapath, True)
    print(f'ONT Files: {len(files)}')
    basecalls = loadFastx(fastxpath)
    print(f'Training segmentation parameters with {len(basecalls)} reads')
    # write first model file with initialized models
    baseName = trainedModels.split('.')[0]
    trainedModels = baseName + '_0_0.model'
    writeKmerModels(trainedModels, kmerModels)
    
    param_writer = open(param_file, 'w')
    error_writer = open(errorfile, 'w')    

    if mode == 'basic':
        mode = 'basic_sparsed'
        transitionParams = {
            'e1': 1.0,
            'm1': 0.03189915859979101,
            'e2': 0.9681008434763126
            }
    elif mode == 'indel':
        transitionParams = {
            'e1': 1.0,
            'm1': 0.047923630044895006,
            'd1': 0.0795918281369677,
            'e2': 0.8719994721033932,
            'i1': 0.00048506999577316453,
            'm2': 0.006990024477599084,
            'i2': 0.019397852195642322,
            'm3': 0.9806021470164152,
            'd2': 0.9930099746422482
            }
    elif mode == '3d' or mode == '3d_reduced':
        mode = '3d_sparsed_reduced' if 'reduced' in mode else '3d_sparsed'
        # a1, a2, p1, p2, p3, s1, s2, s3, e1, e2, e3, e4, i1, i2
        transitionParams = {'a1': 0.030570760232128354, 'a2': 0.6101330423858251, 'p1': 0.0383208584222504, 'p2': 0.08685413259505585, 'p3': 0.015545502978626788, 's1': 0.0010061921517067177, 's2': 0.001534865889565432, 's3': 0.13149698483529376, 'e1': 1.0, 'e2': 0.9989938127816652, 'e3': 0.9616791335752949, 'e4': 0.8801771554834845, 'i1': 0.0008630902783022252, 'i2': 0.24282445342753348}
    else:
        print(f'Mode {mode} not implemented')
        exit(1)

    CPP_SCRIPT = join(dirname(__file__), f'segmentation_{mode}')
    if name == 'nt': # check for windows
        CPP_SCRIPT+='.exe'

    # collect all trained parameters to get an ensemble training in the end
    # paramCollector = {kmer:([kmerModels[kmer][0]], [kmerModels[kmer][1]]) for kmer in kmerModels}
    # paramCollector = paramCollector | {param : [transitionParams[param]] for param in transitionParams}
    paramCollector = {kmer:(ManagedList([kmerModels[kmer][0]]), ManagedList([kmerModels[kmer][1]])) for kmer in kmerModels}
    paramCollector = paramCollector | {param : ManagedList([transitionParams[param]]) for param in transitionParams}
    kmerSeen = set()
    param_writer.write("epoch,batch,read,")
    for param in transitionParams:
        param_writer.write(param+',')
    param_writer.write("Zchange\n")
    i = 0

    with mp.Pool(batch_size) as p:

        for e in range(epochs):
            mp_items = []
            training_readids = []
            batch_num = 0

            for file in files:
                r5 = read(file)
                readids = r5.getReads()
                if random_batch is not None:
                    random.Random(random_batch).shuffle(readids)
                    random_batch+=1 # to avoid that every file is shuffled in the same way

                for readid in readids:
                    if not readid in basecalls: # or (readid not in polyAIndex)
                        continue
                    # skip reads with undetermined transcript start, only take really good reads for training
                    if not readid in polya or polya[readid][1] == -1:
                        continue
                    if polya[readid][1] - polya[readid][0] < 30:
                        continue

                    # fill batch
                    if len(mp_items) < batch_size:
                        if unnormalised:
                            signal = r5.getPolyAStandardizedSignal(readid, polya[readid][0], polya[readid][1])[polya[readid][1]:]
                        else:
                            # signal = r5.getZNormSignal(readid, "median")[polya[readid][1]:]
                            signal = r5.getZNormSignal(readid, "mean")[polya[readid][1]:] # saw more consistency for short reads when using the mean
                        signal = hampel(signal, 6, 5.).filtered_data # small window and high variance allowed: just to filter outliers that result from sensor errors, rest of the original signal should be kept
                        mp_items.append([signal, basecalls[readid][::-1], transitionParams, CPP_SCRIPT, trainedModels, minSegLen, readid])
                        training_readids.append(readid)

                    if len(mp_items) == batch_size:
                        print("============================")
                        print(f"{datetime.now().strftime('%Y-%m-%d_%H-%M-%S')}: Training epoch: {e}, reads: {i}, batch: {batch_num}\n{transitionParams}")
                        print("Training with read:", training_readids)
                        batch_num += 1
                        Zs = []

                        for rid, result in enumerate(p.starmap(trainTransitionsEmissions, mp_items)):
                            if isinstance(result, str):
                                print(f"No segmentation calculated for {result} in {e}: {trainedModels}.")
                                error_writer.write(f"No segmentation calculated for {result} in {e}: {trainedModels}.\n")
                                error_writer.flush()
                                del mp_items[rid]
                                continue

                            trainedParams, newModels, Z = result
                            i += 1
                            Zs.append(Z)

                            for param in trainedParams:
                                # paramCollector[param].append(trainedParams[param])
                                paramCollector[param].add(trainedParams[param])

                            for kmer in newModels:
                                kmerSeen.add(kmer)
                                # paramCollector[kmer][0].append(newModels[kmer][0])
                                # paramCollector[kmer][1].append(newModels[kmer][1])
                                paramCollector[kmer][0].add(newModels[kmer][0])
                                paramCollector[kmer][1].add(newModels[kmer][1])

                        print(f"Zs: {Zs}")

                        # update parameters
                        param_writer.write(f'{e},{batch_num},{i},') # log
                        for param in transitionParams:
                            try:
                                # transitionParams[param] = np.mean(paramCollector[param])
                                transitionParams[param] = paramCollector[param].mean()
                            except:
                                # print(param, paramCollector[param])
                                print(param, paramCollector[param].get_list())
                                exit(1)
                            param_writer.write(f'{transitionParams[param]},') # log

                        for kmer in kmerSeen:
                            # kmerModels[kmer] = [np.mean(paramCollector[kmer][0]), np.mean(paramCollector[kmer][1])]
                            kmerModels[kmer] = [paramCollector[kmer][0].mean(), paramCollector[kmer][1].mean()]

                        trainedModels = baseName + f"_{e}_{batch_num}.model"
                        writeKmerModels(trainedModels, kmerModels)
                        param_writer.flush()

                        # rerun with new parameters to compare Zs
                        for j in range(len(mp_items)):
                            mp_items[j][2] = transitionParams
                            mp_items[j][4] = trainedModels
                        # print(mp_items)
                        Zdiffs = []
                        for j, result in enumerate(p.starmap(calcZ, mp_items)):
                            if isinstance(result, str):
                                print(f"No segmentation calculated for {result} in {e} calcZ.")
                                error_writer.write(f"No segmentation calculated for {result} in {e} calcZ.\n")
                                error_writer.flush()
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

                        if not same_batch:
                            mp_items = []
                            training_readids = []
                        
        error_writer.close()
        param_writer.close()

def main() -> None:
    args = parse()
    outdir = args.out + f'_{datetime.now().strftime("%Y-%m-%d_%H-%M-%S")}'
    # stdev = args.stdev
    if not exists(outdir):
        makedirs(outdir)
    param_file = join(outdir, 'params.csv')
    trainedModels = join(outdir, 'trained.model')
    errorfile = join(outdir, 'error.log')
    polya = readPolyAStartEnd(args.polya)
    # polyAEnd= readPolyAEnd(args.polya)
    # omvKmerModels = getTrainableModels(readKmerModels(args.model_path))
    assert args.minSegLen>=1, "Please choose a minimal segment length greater than 0"
    # due to algorithm min len is 1 by default, minSegLen stores number of positions to add to that length
    train(args.raw, args.fastx, polya, args.batch_size, args.epochs, param_file, args.mode, args.same_batch, args.random_batch, readKmerModels(args.model_path), trainedModels, args.minSegLen - 1, errorfile, args.unnormalised)
    plotParameters(param_file, outdir)

if __name__ == '__main__':
    main()