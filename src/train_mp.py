#!/usr/bin/env python
# author: Jannes Spangenberg
# e-mail: jannes.spangenberg@uni-jena.de
# github: https://github.com/JannesSP
# website: https://jannessp.github.io

from argparse import ArgumentDefaultsHelpFormatter, ArgumentParser, Namespace
import numpy as np
from os.path import exists, join, dirname
from os import makedirs, name
from fileio import getFiles, loadFastx, trainSegmentation, calcZ, readPolyAEnd, plotParameters
from read5 import read
import multiprocessing as mp
from hampel import hampel

CPP_SCRIPT = join(dirname(__file__), 'segmentation')
if name == 'nt': # check for windows
    CPP_SCRIPT+='.exe'
print(f"training with {CPP_SCRIPT} script")

def parse() -> Namespace:
    parser = ArgumentParser(
        formatter_class=ArgumentDefaultsHelpFormatter
    )
    parser.add_argument('--raw', type=str, required=True, help='Raw ONT training data')
    parser.add_argument('--fastx', type=str, required=True, help='Basecalls of ONT training data')
    parser.add_argument('--polya', type=str, required=True, help='Poly A table from nanopolish polya containing the transcript starts')
    parser.add_argument('--out', type=str, required=True, help='Outpath to write files')
    parser.add_argument('--batch_size', type=int, default=16, help='Number of reads to train before updating')
    parser.add_argument('--epochs', type=int, default=8, help='Number of training epochs')
    return parser.parse_args()

def train(rawdatapath : str, fastxpath : str, polya : dict, batch_size : int, epochs :int, param_file : str) -> None:
    files = getFiles(rawdatapath, True)
    print(f'ONT Files: {len(files)}')
    basecalls = loadFastx(fastxpath)
    print(f'Training segmentation parameters with {len(basecalls)} reads')
    
    param_writer = open(param_file, 'w')
    
    # init
    params = {
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
    paramCollector = {param : 0 for param in params}
    param_writer.write("epoch,batch,read,")
    for param in params:
        param_writer.write(param+',')
    param_writer.write("Zchange\n")
    # pipe = openCPPScriptParamsTrain(CPP_SCRIPT, params)
    i = 0
    batch_num = 0
    # TODO correct for Z == -inf, should this happen? maybe a bug?
    failedInBatch = 0

    with mp.Pool(batch_size) as p:

        for e in range(epochs):
            for file in files:

                r5 = read(file)
                mp_items = []

                for readid in r5:
                    if not readid in basecalls: # or (readid not in polyAIndex)
                        continue
                    # skip reads with undetermined transcript start
                    # TODO improvement room here
                    if not readid in polya:
                        continue

                    if polya[readid] == -1:
                        polya[readid] = 0

                    signal = hampel(r5.getpASignal(readid)[polya[readid]:], 20, 2.).filtered_data
                    mp_items.append([signal, basecalls[readid][::-1], params, CPP_SCRIPT])

                    if len(mp_items) == batch_size:
                        batch_num += 1
                        Zs = []

                        for result in p.starmap(trainSegmentation, mp_items):
                            trainedParams, Z = result
                            i += 1
                            Zs.append(Z)

                            if not np.isinf(Z):
                                param_writer.write(f'{e},{batch_num},{i},') # log
                                for j, param in enumerate(trainedParams):
                                    param_writer.write(f'{params[param]},') # log
                                    paramCollector[param] += trainedParams[param]
                            else:
                                failedInBatch += 1

                        # update parameters
                        for param in params:
                            params[param] = paramCollector[param] / (batch_size - failedInBatch)

                        # TODO rerun with new parameters to compare Zs
                        for j in range(len(mp_items)):
                            mp_items[j][2] = params
                        Zdiffs = []
                        for j, result in enumerate(p.starmap(calcZ, mp_items)):
                            Z = result
                            if not np.isinf(Zs[j]):
                                Zdiffs.append(Z - Zs[j])

                        param_writer.write(f'{sum(Zdiffs)/len(Zdiffs)}\n') # log
                        param_writer.flush() # log
                        print(f"Training epoch: {e}, reads: {i}, batch: {batch_num}, failed: {failedInBatch}\n{params}\nZ change: {sum(Zdiffs)/len(Zdiffs)}")

                        # initialize new batch
                        paramCollector = {param : 0 for param in paramCollector}
                        failedInBatch = 0
                        mp_items = []
        
        param_writer.close()

def main() -> None:
    args = parse()
    outdir = args.out
    if not exists(outdir):
        makedirs(outdir)
    param_file = join(outdir, 'params.txt')
    polya=readPolyAEnd(args.polya)
    train(args.raw, args.fastx, polya, args.batch_size, args.epochs, param_file)
    plotParameters(param_file, outdir)

if __name__ == '__main__':
    main()