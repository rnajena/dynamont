#!/usr/bin/env python
# author: Jannes Spangenberg
# e-mail: jannes.spangenberg@uni-jena.de
# github: https://github.com/JannesSP
# website: https://jannessp.github.io

import numpy as np
import read5_ont
import multiprocessing as mp
import pysam
from argparse import ArgumentDefaultsHelpFormatter, ArgumentParser, Namespace
from os.path import exists, join, dirname
from os import makedirs, name
from datetime import datetime
from collections import deque
from src.python.segmentation.FileIO import calcZ, plotParameters, trainTransitionsEmissions, readKmerModels, writeKmerModels, hampelFilter, countNucleotideRatios

class ManagedList:
    def __init__(self, values, max_size=100):
        """Initialize the ManagedList with a maximum size."""
        self.values = deque(values, maxlen=max_size)  # deque automatically limits size

    def add(self, value):
        """Add a new value to the list, removing the oldest if necessary."""
        self.values.append(value)

    def get_list(self):
        """Returns the current list of values."""
        return list(self.values)

    def __repr__(self):
        """String representation of the ManagedList."""
        return f"ManagedList({list(self.values)})"
    
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

def parse() -> Namespace:
    parser = ArgumentParser(
        formatter_class=ArgumentDefaultsHelpFormatter
    )
    # required
    parser.add_argument('-r', '--raw',   type=str, required=True, metavar="PATH", help='Path to raw ONT data. [POD5|FAST5]')
    parser.add_argument('-b', '--basecalls', type=str, required=True, metavar="BAM", help='Basecalls of ONT training data as .bam file')
    parser.add_argument('-o', '--outdir',   type=str, required=True, metavar="PATH", help='Outpath to write files')
    parser.add_argument('-p', '--pore',  type=str, required=True, choices=["rna_r9", "dna_r9", "rna_rp4", "dna_r10_260bps", "dna_r10_400bps"], help='Pore generation used to sequence the data')
    parser.add_argument('--mode',  type=str, required=True, choices=['basic', 'resquiggle'], help='Segmentation algorithm used for segmentation')
    # optional
    parser.add_argument('--model_path', type=str, required=True, help='Which initial kmer models to use for training')
    parser.add_argument('--batch_size', type=int, default=24, help='Number of reads to train before updating')
    parser.add_argument('--max_batches', type=int, default=None, help='Numbers of batches to train each epoch')
    parser.add_argument('-e', '--epochs', type=int, default=1, help='Number of training epochs')
    parser.add_argument('-q', '--qscore', type=int, default=10, help='Minimal allowed quality score')
    return parser.parse_args()

def train(dataPath : str, basecalls : str, batch_size : int, epochs :int, param_file : str, mode : str, model_path : str, maxBatch : int, pore : str, minQual : float = None) -> None:

    kmerModels = readKmerModels(model_path)
    trainedModels = join(dirname(param_file), 'trained_0_0.model')
    writeKmerModels(trainedModels, kmerModels)
    
    paramWriter = open(param_file, 'w')

    if mode == 'basic':
        # CPP_SCRIPT = 'dynamont-NT'
        CPP_SCRIPT = 'dynamont-NT-banded'
        transitionParams = {
            'e1': 1.0,
            'm1': 0.03,
            'e2': 0.97
            }
    elif mode == 'resquiggle':
        CPP_SCRIPT = 'dynamont-NTC'
        transitionParams = {
            'a1': 0.012252440188168037,
            'a2': 0.246584724985145,
            'p1': 0.04477093133243305,
            'p2': 0.007687811003133089,
            'p3': 0.4469623669791557,
            's1': 0.05321209670114726,
            's2': 0.0007555035568187239,
            's3': 0.21999557711272136,
            'e1': 1.0,
            'e2': 0.9467879033992115,
            'e3': 0.9552290685034269,
            'e4': 0.9792321612614708,
            'i1': 7.208408117990252e-05,
            'i2': 0.08645733058947891
            }
    else:
        print(f'Mode {mode} not implemented')
        exit(1)

    if name == 'nt': # check for windows
        CPP_SCRIPT+='.exe'

    # collect trained parameters to get an ensemble training in the end to reduce outlier trainings
    paramCollector = {kmer:(ManagedList([kmerModels[kmer][0]]), ManagedList([kmerModels[kmer][1]])) for kmer in kmerModels}
    paramCollector = paramCollector | {param : ManagedList([transitionParams[param]]) for param in transitionParams}
    kmerSeen = set()
    paramWriter.write("epoch,batch,read,")
    for param in transitionParams:
        paramWriter.write(param+',')
    paramWriter.write("Zchange\n")
    i = 0
    qualSkipped = 0
    noMatchingReadid = 0

    with mp.Pool(batch_size) as p:

        for e in range(epochs):
            mpItems = []
            trainIDs = []
            batchNum = 0
            oldFile = None

            with pysam.AlignmentFile(basecalls, "r" if basecalls.endswith('.sam') else "rb", check_sq=False) as samfile:
                for basecalledRead in samfile.fetch(until_eof=True):
                    
                    # skip low qual reads if activated
                    qual = basecalledRead.get_tag("qs")
                    if minQual and qual < minQual:
                        qualSkipped+=1
                        continue

                    # init read
                    seq = basecalledRead.query_sequence
                    counts = countNucleotideRatios(seq)

                    # saw weird signals in rp4 data, a very homogenous signal jumping between two values
                    # does not look like a normal read, more like an artifact
                    # produces reads with high quality but the read consists only of repeats
                    if any(counts[base] >= 0.6 for base in counts.keys()):
                        continue

                    # init read, sometimes a read got split by the basecaller and got a new id
                    readid = basecalledRead.get_tag("pi") if basecalledRead.has_tag("pi") else basecalledRead.query_name
                    sp = basecalledRead.get_tag("sp") if basecalledRead.has_tag("sp") else 0 # if split read get start offset of the signal
                    ns = basecalledRead.get_tag("ns") # ns:i: 	the number of samples in the signal prior to trimming
                    ts = basecalledRead.get_tag("ts") # ts:i: 	the number of samples trimmed from the start of the signal
                    rawFile = join(dataPath, basecalledRead.get_tag("fn"))

                    if oldFile != rawFile:
                        oldFile = rawFile
                        r5 = read5_ont.read(rawFile)

                    # fill batch
                    if len(mpItems) < batch_size:
                        # saw more consistency for short reads when using the mean
                        try:
                            if pore in ["dna_r9", "rna_r9"]:
                                # for r9 pores, shift and scale are stored for pA signal in bam
                                signal = r5.getpASignal(readid)

                                # signal = (signal - basecalledRead.get_tag("sm")) / basecalledRead.get_tag("sd")

                            else:
                                # for new pores, shift and scale is directly applied to stored integer signal (DACs)
                                # this way the conversion from DACs to pA is skipped
                                signal = r5.getSignal(readid)
                                
                                # norm_signal = (signal - basecalledRead.get_tag("sm")) / basecalledRead.get_tag("sd")
                                # # slice signal, remove remaining adapter content until polyA starts
                                # start = np.argmax(norm_signal[sp+ts:] >= 0.8) + sp + ts + 100
                                # signal = signal[start : sp+ns]

                                #! normalize poly A region to median 0.9 (as in init models from ONT r9 and rp4) and scale to 0.15 (from training on r9 and rp4)
                                # shift = np.median(signal[:1000])
                                # scale = np.median(np.absolute(signal[:1000] - shift))
                                # signal = ((signal - shift) / scale) * 0.15 + 0.9


                            #! normalize whole signal
                            # mad = np.median(np.absolute(signal - np.median(signal)))
                            # if mad < 70: # check again for weird repetitive reads/signals with high quality
                            #     continue
                            signal = (signal - basecalledRead.get_tag("sm")) / basecalledRead.get_tag("sd")
                            signal = signal[sp+ts : sp+ns]
                        
                        except:
                            noMatchingReadid+=1
                            continue
                        # signal = hampel(signal, 6, 5.).filtered_data # small window and high variance allowed: just to filter outliers that result from sensor errors, rest of the original signal should be kept
                        hampelFilter(signal, 6, 5.) # small window and high variance allowed: just to filter outliers that result from sensor errors, rest of the original signal should be kept
                        if "rna" in pore:
                            seq = seq[::-1]
                            if not seq.startswith("AAAAAAAAA"):
                                seq = "AAAAAAAAA" + seq
                        mpItems.append([signal, seq, transitionParams | {'r' : pore, 't' : 4}, CPP_SCRIPT, trainedModels, readid])
                        trainIDs.append(readid)

                    if len(mpItems) == batch_size:
                        print("============================")
                        print(f"{datetime.now().strftime('%Y-%m-%d_%H-%M-%S')}: Training epoch: {e}, reads: {i}, batch: {batchNum}\n{transitionParams}")
                        print("Training with read:", trainIDs)
                        batchNum += 1
                        preZ = np.zeros(batch_size)

                        for readid, result in enumerate(p.starmap(trainTransitionsEmissions, mpItems)):
                            if isinstance(result, str):
                                print(f"No segmentation calculated for {result} in {e}: {trainedModels}.")
                                # del mp_items[readid]
                                continue

                            trainedParams, newModels, Z = result

                            i += 1
                            preZ[readid] = Z

                            for param in trainedParams:
                                paramCollector[param].add(trainedParams[param])

                            #! skip weird trainings
                            if ('AAAAAAAAA' in newModels and newModels['AAAAAAAAA'][0] < 0.5) or ('AAAAA' in newModels and newModels['AAAAA'][0] < 0.5):
                                continue

                            for kmer in newModels:
                                kmerSeen.add(kmer)
                                paramCollector[kmer][0].add(newModels[kmer][0])
                                paramCollector[kmer][1].add(newModels[kmer][1])

                        print(f"Zs: {preZ}")

                        # update parameters
                        paramWriter.write(f'{e},{batchNum},{i},') # log
                        for param in transitionParams:
                            try:
                                transitionParams[param] = paramCollector[param].mean()
                            except:
                                print(param, paramCollector[param].get_list())
                                exit(1)
                            paramWriter.write(f'{transitionParams[param]},') # log

                        for kmer in kmerSeen:
                            kmerModels[kmer] = [paramCollector[kmer][0].mean(), paramCollector[kmer][1].mean()]

                        trainedModels = join(dirname(trainedModels), f"trained_{e}_{batchNum}.model")
                        writeKmerModels(trainedModels, kmerModels)
                        paramWriter.flush()

                        # rerun with new parameters to compare Zs
                        for j in range(len(mpItems)):
                            mpItems[j][2] = transitionParams | {'r' : pore}
                            mpItems[j][4] = trainedModels
                        postZ = np.zeros(batch_size)
                        for j, result in enumerate(p.starmap(calcZ, mpItems)):
                            if isinstance(result, str):
                                print(f"No segmentation calculated for {result} in {e} calcZ.")
                                continue
                            Z = result
                            postZ[j] = Z

                        dZ = postZ - preZ

                        print(f"Z changes: {dZ}")
                        deltaZ = np.mean(dZ)
                        paramWriter.write(f'{deltaZ}\n') # log
                        paramWriter.flush() # log
                        # initialize new batch
                        kmerSeen = set()
                        mpItems = []
                        trainIDs = []

                        if maxBatch is not None and batchNum >= maxBatch:
                            break

        paramWriter.close()
        print("Done training")
        print(f"Skipped reads due to low quality: {qualSkipped}")

def main() -> None:
    args = parse()

    # init outdir
    outdir = args.outdir + f'_{datetime.now().strftime("%Y-%m-%d_%H-%M-%S")}'
    if not exists(outdir):
        makedirs(outdir)

    paramFile = join(outdir, 'params.csv')

    train(args.raw, args.basecalls, args.batch_size, args.epochs, paramFile, args.mode, args.model_path, args.max_batches, args.pore, args.qscore)
    plotParameters(paramFile, outdir)

if __name__ == '__main__':
    main()