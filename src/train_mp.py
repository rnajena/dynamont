#!/usr/bin/env python
# author: Jannes Spangenberg
# e-mail: jannes.spangenberg@uni-jena.de
# github: https://github.com/JannesSP
# website: https://jannessp.github.io

from argparse import ArgumentDefaultsHelpFormatter, ArgumentParser, Namespace
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from os.path import exists, join, dirname
from os import makedirs, name
from fileio import getFiles, loadFastx, openCPPScriptParamsTrain, stopFeeding
from read5 import read
import multiprocessing as mp

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
    parser.add_argument('--batch_size', type=int, default=1, help='Number of reads to train before updating')
    parser.add_argument('--epochs', type=int, default=1, help='Number of training epochs')
    return parser.parse_args()

# https://stackoverflow.com/questions/32570029/input-to-c-executable-python-subprocess
def trainSegmentation(signal : np.ndarray, read : str, params : dict) -> np.ndarray:
    '''
    Parse & feed signal & read to the C++ segmentation script.

    Parameters
    ----------
    signal : np.ndarray
    read : str
    stream
        Open stdin stream of the C++ segmentation algorithm

    Returns
    -------
    params : dict
        {str : float}
    Z : float
    '''
    pipe = openCPPScriptParamsTrain(CPP_SCRIPT, params)
    # prepare cookie for segmentation
    cookie = f"{str(signal.tolist()).replace(' ', '').replace('[', '').replace(']', '')} {read}\n"
    c = open(join(dirname(__file__), 'last_cookie.txt'), 'w')
    c.write(cookie)
    c.close()
    # transfer data to bytes - needed in Python 3
    cookie = bytes(cookie, 'UTF-8')
    # feed cookie to segmentation
    pipe.stdin.write(cookie)
    pipe.stdin.flush()
    output = pipe.stdout.readline().strip().decode('UTF-8')
    try:
        params = {param.split(":")[0] : float(param.split(":")[1]) for param in output.split(";")}
    except Exception as e:
        print(output)
        print(pipe.stdout.readlines())
        print(e)
        exit(1)

    Z = float(pipe.stdout.readline().strip().decode('UTF-8').split(':')[1])
    stopFeeding(pipe)

    return params, Z

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
    paramCollector["Z"]=0
    paramCollector["ZpB"]=0
    param_writer.write("epoch,batch,read,")
    for param in params:
        param_writer.write(param+',')
    param_writer.write("Z,ZpB\n")
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

                    mp_items.append((r5.getpASignal(readid)[polya[readid]:], basecalls[readid][::-1], params))

                    if len(mp_items) == batch_size:
                        batch_num += 1

                        for result in p.starmap(trainSegmentation, mp_items):
                            trainedParams, Z = result
                            i += 1
                            if not np.isinf(Z):
                                param_writer.write(f'{e},{batch_num},{i},') # log
                                for j, param in enumerate(trainedParams):
                                    param_writer.write(f'{params[param]},') # log
                                    paramCollector[param] += trainedParams[param]
                                paramCollector["Z"]+=Z
                                paramCollector["ZpB"]+=Z/len(mp_items[j][0])
                                param_writer.write(f'{Z},{Z/len(mp_items[j][0])}\n') # log
                                param_writer.flush() # log
                            else:
                                failedInBatch += 1

                        # update parameters
                        for param in params:
                            params[param] = paramCollector[param] / (batch_size - failedInBatch)
                        print(f"Training epoch: {e}, reads: {i}, batch: {batch_num}, failed: {failedInBatch}\n{params}\nbatch avg Z: {paramCollector['Z'] / (batch_size - failedInBatch)}, batch avg ZpB: {paramCollector['ZpB'] / (batch_size - failedInBatch)}")

                        # initialize new batch
                        paramCollector = {param : 0 for param in paramCollector}
                        failedInBatch = 0
                        mp_items = []
        
        param_writer.close()

def plot(param_file : str, outdir : str) -> None:
    df = pd.read_csv(param_file, sep=',')
    for column in df:
        if column in ['epoch', 'batch']:
            continue
        sns.set_theme()
        sns.lineplot(data=df, x="batch", y=column, hue='epoch')
        plt.title(f"{column} parameter change during training")
        plt.savefig(join(outdir, f"{column}.pdf"))
        plt.cla()
        plt.close()

def readPolyA(file : str) -> dict:
    df = pd.read_csv(file, usecols=['readname', 'transcript_start'], sep='\t')
    df = df.astype({'readname' : str, 'transcript_start' : int})
    # df.set_index('readname', inplace=True)
    return pd.Series(df.transcript_start.values, index=df.readname).to_dict()

def main() -> None:
    args = parse()
    outdir = args.out
    if not exists(outdir):
        makedirs(outdir)
    param_file = join(outdir, 'params.txt')
    polya=readPolyA(args.polya)
    train(args.raw, args.fastx, polya, args.batch_size, args.epochs, param_file)
    plot(param_file, outdir)

if __name__ == '__main__':
    main()