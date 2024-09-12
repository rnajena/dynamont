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

INITMODEL = None

def parse() -> Namespace:
    parser = ArgumentParser(
        formatter_class=ArgumentDefaultsHelpFormatter
    )
    parser.add_argument('--raw', type=str, required=True, help='Raw ONT training data')
    parser.add_argument('--fastx', type=str, required=True, help='Basecalls of ONT training data')
    parser.add_argument('--polya', type=str, required=True, help='Poly A table from nanopolish polya containing the transcript starts')
    parser.add_argument('--out', type=str, required=True, help='Outpath to write files')
    parser.add_argument('--mode', type=str, choices=['basic', 'indel', '3d'], required=True, help='Segmentation algorithm used for segmentation')
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
    global INITMODEL
    writeKmerModels(trainedModels, kmerModels, INITMODEL)
    
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
    elif mode == '3d':
        mode = '3d_sparsed'
        # a1, a2, p1, p2, p3, s1, s2, s3, e1, e2, e3, e4, i1, i2
        transitionParams = {
            'a1': 0.007025917258946322,
            'a2': 0.02163382898644135,
            'p1': 0.9136072624254473,
            'p2': 0.7577089930417495,
            'p3': 0.031412750592942344,
            's1': 0.05633093742544732,
            's2': 0.003601331413021869,
            's3': 0.6537826547572912,
            'e1': 1.0,
            'e2': 0.9436690661033798,
            'e3': 0.08639273593439363,
            'e4': 0.23098202335984094,
            'i1': 0.0006817432661033797,
            'i2': 0.2931707726206201,
            }
    else:
        print(f'Mode {mode} not implemented')
        exit(1)

    CPP_SCRIPT = join(dirname(__file__), f'segmentation_{mode}')
    if name == 'nt': # check for windows
        CPP_SCRIPT+='.exe'

    # collect all trained parameters to get an ensemble training in the end
    paramCollector = {kmer:([kmerModels[kmer][0]], [kmerModels[kmer][1]]) for kmer in kmerModels}
    paramCollector = paramCollector | {param : [transitionParams[param]] for param in transitionParams}
    kmerSeen = set()
    param_writer.write("epoch,batch,read,")
    for param in transitionParams:
        param_writer.write(param+',')
    param_writer.write("Zchange\n")
    i = 0
    batch_num = 0

    with mp.Pool(batch_size) as p:

        for e in range(epochs):
            mp_items = []
            training_readids = []

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
                                paramCollector[param].append(trainedParams[param])

                            for kmer in newModels:
                                kmerSeen.add(kmer)
                                paramCollector[kmer][0].append(newModels[kmer][0])
                                paramCollector[kmer][1].append(newModels[kmer][1])

                        print(f"Zs: {Zs}")

                        # update parameters
                        param_writer.write(f'{e},{batch_num},{i},') # log
                        for param in transitionParams:
                            try:
                                transitionParams[param] = np.mean(paramCollector[param])
                            except:
                                print(param, paramCollector[param])
                                exit(1)
                            param_writer.write(f'{transitionParams[param]},') # log

                        for kmer in kmerSeen:
                            kmerModels[kmer] = [np.mean(paramCollector[kmer][0]), np.mean(paramCollector[kmer][1])]
                        
                        trainedModels = baseName + f"_{e}_{batch_num}.model"
                        writeKmerModels(trainedModels, kmerModels, INITMODEL)
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
    global INITMODEL
    INITMODEL = readKmerModels(args.model_path)
    assert args.minSegLen>=1, "Please choose a minimal segment length greater than 0"
    # due to algorithm min len is 1 by default, minSegLen stores number of positions to add to that length
    train(args.raw, args.fastx, polya, args.batch_size, args.epochs, param_file, args.mode, args.same_batch, args.random_batch, INITMODEL.copy(), trainedModels, args.minSegLen - 1, errorfile, args.unnormalised)
    plotParameters(param_file, outdir)

if __name__ == '__main__':
    main()