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
# from hampel import hampel
import random
# import OnlineMeanVar as omv

INITMODEL = None

def parse() -> Namespace:
    parser = ArgumentParser(
        formatter_class=ArgumentDefaultsHelpFormatter
    )
    parser.add_argument('--raw', type=str, required=True, help='Raw ONT training data')
    parser.add_argument('--fastx', type=str, required=True, help='Basecalls of ONT training data')
    parser.add_argument('--polya', type=str, required=True, help='Poly A table from nanopolish polya containing the transcript starts')
    parser.add_argument('--out', type=str, required=True, help='Outpath to write files')
    parser.add_argument('--mode', type=str, choices=['basic', 'indel'], required=True, help='Segmentation algorithm used for segmentation')
    parser.add_argument('--model_path', type=str, default=join(dirname(__file__), '..', '..', 'data', 'template_median69pA_extended.model'), help='Which models to train')
    parser.add_argument('--pore', type=int, choices=[9, 10], default=9, help='Pore generation used to sequence the data')
    parser.add_argument('--batch_size', type=int, default=24, help='Number of reads to train before updating')
    parser.add_argument('--epochs', type=int, default=1, help='Number of training epochs')
    parser.add_argument('--same_batch', action='store_true', help='Always train on first batch of signals in the file')
    parser.add_argument('--random_batch', type=int, metavar='seed', default=None, help='Randomizes the batches of each file with a given seed')
    parser.add_argument('--stdev', type=float, default=None, help='Set a stdev value instead of using the calculated one for the models')
    parser.add_argument('--minSegLen', type=int, default=1, help='Minmal allowed segment length')
    return parser.parse_args()

# def getTrainableModels(kmerModels : dict) -> dict:
#     # init online mean var
#     trainableModels = {kmer : omv.LocShift(120) for kmer in kmerModels}
#     # fill with drawn values of original paramters
#     for kmer in kmerModels:
#         trainableModels[kmer].append(np.random.normal(loc=kmerModels[kmer][0], scale=kmerModels[kmer][1], size=1000))
#     return trainableModels

def train(rawdatapath : str, fastxpath : str, polya : dict, batch_size : int, epochs :int, param_file : str, mode : str, same_batch : bool, random_batch : bool, kmerModels : dict, trainedModels : str, stdev : float, minSegLen : int) -> None:
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
    
    CPP_SCRIPT = join(dirname(__file__), f'segmentation_{mode}')
    if name == 'nt': # check for windows
        CPP_SCRIPT+='.exe'

    if mode == 'indel':
        # init
        transitionParams = {
            "e1":1.,
            "m2":.03333,
            "d1":.00001,
            "e2":.96664,
            "e3":.00001,
            "i1":.00001,
            "m3":.99,
            "i2":.01,
            "m4":.99,
            "d2":.01,
            # "s1":0.16,
            # "s2":1.0
        }
    elif mode == 'basic':
        transitionParams = {'e1': 1.0, 'm1': 0.111111, 'e2': 0.888889, 'e3': 0.0}
        # {
        #     "e1":1.,
        #     "m1":.025,
        #     "e2":.955,
        #     "e3":.02,
        # }
    elif mode == 'extended':
        # TODO
        print('Extended mode not implemented')
        exit(1)

    paramCollector = {param : [] for param in transitionParams}
    meanCollector = {kmer : [] for kmer in kmerModels}
    stdevCollector = {kmer : [] for kmer in kmerModels}
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

                    # if polya[readid] == -1:
                    #     polya[readid] = 0

                    # fill batch
                    if len(mp_items) < batch_size:
                        # signal = hampel(r5.getPolyAStandardizedSignal(readid, polya[readid][0], polya[readid][1])[polya[readid][1]:], 50, 2.).filtered_data
                        signal = r5.getPolyAStandardizedSignal(readid, polya[readid][0], polya[readid][1])[polya[readid][1]:]
                        mp_items.append([signal, basecalls[readid][::-1], transitionParams, CPP_SCRIPT, trainedModels, minSegLen])
                        training_readids.append(readid)

                    if len(mp_items) == batch_size:
                        print("============================")
                        print(f"Training epoch: {e}, reads: {i}, batch: {batch_num}\n{transitionParams}")
                        print("Training with read:", training_readids)
                        batch_num += 1
                        Zs = []
                        segmentations = []

                        for result in p.starmap(trainTransitionsEmissions, mp_items):
                            segments, trainedParams, newModels, Z = result
                            i += 1
                            Zs.append(Z)
                            segmentations.append(segments)

                            assert not np.isinf(Z), f'Z is infinit!: {Z}'
                            for j, param in enumerate(trainedParams):
                                paramCollector[param].append(trainedParams[param])

                            for kmer in newModels:
                                meanCollector[kmer].append(newModels[kmer][0])
                                stdevCollector[kmer].append(newModels[kmer][1])
                        print(f"Zs: {Zs}")

                        # update parameters
                        param_writer.write(f'{e},{batch_num},{i},') # log
                        for param in transitionParams:
                            transitionParams[param] = np.mean(paramCollector[param])
                            param_writer.write(f'{transitionParams[param]},') # log

                        for kmer in meanCollector:
                            # skip unseen models
                            if not len(meanCollector[kmer]):
                                continue
                            # TODO what is better here, median oder mean?
                            kmerModels[kmer] = [np.median(meanCollector[kmer]), np.sqrt(np.median(np.square(stdevCollector[kmer])))]
                        
                        trainedModels = baseName + f"_{e}_{batch_num}.model"
                        writeKmerModels(trainedModels, kmerModels, INITMODEL)

                        # rerun with new parameters to compare Zs
                        for j in range(len(mp_items)):
                            mp_items[j][2] = transitionParams
                            mp_items[j][4] = trainedModels
                        # print(mp_items)
                        Zdiffs = []
                        for j, result in enumerate(p.starmap(calcZ, mp_items)):
                            Z = result
                            # print('Z:', Z, end='\t')
                            # if not np.isinf(Zs[j]):
                            Zdiffs.append(Z - Zs[j])

                        print(f"Z changes: {Zdiffs}")
                        deltaZ = np.mean(Zdiffs)
                        param_writer.write(f'{deltaZ}\n') # log
                        param_writer.flush() # log
                        # print(f"Z change: {deltaZ}")

                        # End training if delta Z is 0
                        # if deltaZ <= 0:
                        #     return

                        # initialize new batch
                        paramCollector = {param : [] for param in paramCollector}
                        meanCollector = {kmer : [] for kmer in kmerModels}
                        stdevCollector = {kmer : [] for kmer in kmerModels}

                        if not same_batch:
                            mp_items = []
                            training_readids = []
                        
        
        param_writer.close()

def main() -> None:
    args = parse()
    outdir = args.out
    stdev = args.stdev
    if not exists(outdir):
        makedirs(outdir)
    param_file = join(outdir, 'params.txt')
    trainedModels = join(outdir, 'trained.model')
    polya = readPolyAStartEnd(args.polya)
    # polyAEnd= readPolyAEnd(args.polya)
    # omvKmerModels = getTrainableModels(readKmerModels(args.model_path))
    global INITMODEL
    INITMODEL = readKmerModels(args.model_path)
    assert args.minSegLen>=1, "Please choose a minimal segment length greater than 0"
    # due to algorithm min len is 1 by default, minSegLen stores number of positions to add to that length
    train(args.raw, args.fastx, polya, args.batch_size, args.epochs, param_file, args.mode, args.same_batch, args.random_batch, INITMODEL.copy(), trainedModels, stdev, args.minSegLen - 1)
    plotParameters(param_file, outdir)

if __name__ == '__main__':
    main()