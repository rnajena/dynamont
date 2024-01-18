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
from fileio import getFiles, loadFastqs, openCPPScriptParamsTrain, stopFeeding
from read5 import read
from subprocess import Popen

def parse() -> Namespace:
    parser = ArgumentParser(
        formatter_class=ArgumentDefaultsHelpFormatter
    )
    parser.add_argument('--raw', type=str, required=True, help='Raw ONT training data')
    parser.add_argument('--fastq', type=str, required=True, help='Basecalls of ONT training data')
    parser.add_argument('--out', type=str, required=True, help='Outpath to write files')
    parser.add_argument('--batch_size', type=int, default=1, help='Number of reads to train before updating')
    parser.add_argument('--epochs', type=int, default=1, help='Number of training epochs')
    return parser.parse_args()

# https://stackoverflow.com/questions/32570029/input-to-c-executable-python-subprocess
def trainSegmentation(signal : np.ndarray, read : str, pipe : Popen) -> np.ndarray:
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
    # prepare cookie for segmentation
    cookie = f"{str(signal.tolist()).replace(' ', '').replace('[', '').replace(']', '')} {read}\n"
    c = open(join(dirname(__file__), 'last_cookie.txt'), 'w')
    c.write(cookie)
    c.close()
    # transfer data to bytes - needed in Python 3
    cookie = bytes(cookie, 'UTF-8')
    # feed cookie to segmentation
    # try:
    pipe.stdin.write(cookie)
    # except BrokenPipeError as e:
    #     print(e)
    #     print(cookie)
    #     exit(1)
    pipe.stdin.flush()
    output = pipe.stdout.readline().strip().decode('UTF-8')
    # print("output:", output)
    params = {param.split(":")[0] : float(param.split(":")[1]) for param in output.split(";")}
    # print("params:", params)
    Z = float(pipe.stdout.readline().strip().decode('UTF-8'))
    # print("Z:", Z)
    
    # print("P matrix")
    # print(pipe.stdout.readline())
    # print("Exit")
    # exit(1)

    return params, Z

def train(rawdatapath : str, fastqpath : str, outdir : str, batch_size : int, epochs :int, param_file : str) -> None:
    files = getFiles(rawdatapath, True)
    print(f'ONT Files: {len(files)}')
    basecalls = loadFastqs(fastqpath)
    print(f'Training segmentation parameters with {len(basecalls)} reads')
    
    param_writer = open(param_file, 'w')

    CPP_SCRIPT = join(dirname(__file__), 'segment_affine_deletion')
    if name == 'nt': # check for windows
        CPP_SCRIPT+='.exe'
    
    # init
    params = {
        "p":0.5,
        "m1":0.5,
        "e1":1.0,
        "s1":0.16,
        "m2":0.16,
        "d1":0.16,
        "e2":0.16,
        "e3":0.16,
        "i1":0.16,
        "m3":0.5,
        "i2":0.5,
        "m4":0.5,
        "d2":0.5,
        "s2":1.0
    }
    paramCollector = {param : 0 for param in params}
    paramCollector["Z"]=0
    param_writer.write("epoch,batch,")
    for param in params:
        param_writer.write(param+',')
    param_writer.write("Z\n")
    pipe = openCPPScriptParamsTrain(CPP_SCRIPT, params)
    i = -1

    for e in range(epochs):
        for file in files:
            r5 = read(file)
            for readid in r5:
                if not readid in basecalls: # or (readid not in polyAIndex)
                    continue
                i+=1
                print('i:', i)
                if i % batch_size == 0 and i > 0:
                    print(f"Training epoch: {e}, batch: {i // batch_size}")
                    
                    print(paramCollector)
                    # log old parameters
                    param_writer.write(f'{e},{i // batch_size},')
                    for param in params:
                        param_writer.write(f'{params[param]},')
                    param_writer.write(f'{paramCollector["Z"]/batch_size}\n')
                    param_writer.flush()

                    # update parameters
                    for param in params:
                        params[param] = paramCollector[param] / batch_size
                    
                    # initialize new batch
                    paramCollector = {param : 0 for param in paramCollector}
                    stopFeeding(pipe)
                    pipe = openCPPScriptParamsTrain(CPP_SCRIPT, params)

                trainedParams, Z = trainSegmentation(r5.getpASignal(readid), basecalls[readid][::-1], pipe)
                print("Z:", Z)
                for param in trainedParams:
                    paramCollector[param] += trainedParams[param]
                paramCollector["Z"]+=Z
    
    # stopFeeding(pipe)
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

def main() -> None:
    args = parse()
    outdir = args.out
    if not exists(outdir):
        makedirs(outdir)
    param_file = join(outdir, 'params.txt')
    train(args.raw, args.fastq, args.out, args.batch_size, args.epochs, param_file)
    plot(param_file, outdir)

if __name__ == '__main__':
    main()