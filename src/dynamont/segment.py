#!/usr/bin/env python
# author: Jannes Spangenberg
# e-mail: jannes.spangenberg@uni-jena.de
# github: https://github.com/JannesSP
# website: https://jannessp.github.io

from argparse import ArgumentDefaultsHelpFormatter, ArgumentParser, Namespace
from os.path import exists, join, dirname
from os import makedirs, name

import read5.AbstractFileReader
import read5.Pod5Reader
from FileIO import getFiles, loadFastx, feedSegmentationAsynchronous
import read5
import multiprocessing as mp

TERM_STRING = "$"

def parse() -> Namespace:
    parser = ArgumentParser(
        formatter_class=ArgumentDefaultsHelpFormatter
    )
    parser.add_argument('--raw', type=str, required=True, help='Raw ONT training data')
    parser.add_argument('--fastx', type=str, required=True, help='Basecalls of ONT training data')
    # parser.add_argument('--polya', type=str, required=True, help='Poly A table from nanopolish polya containing the transcript starts')
    parser.add_argument('-o', '--outfile', type=str, required=True, help='Outpath to write files')
    parser.add_argument('-p', '--processes', type=int, default=None, help='Number of processes to use for segmentation')
    parser.add_argument('--model', type=str, default=None)
    parser.add_argument('-m', '--mode', choices=['basic', 'banded', 'resquiggle'], required=True)
    parser.add_argument('--pore',  type=str, required=True, choices=["rna_r9", "dna_r9", "rna_rp4", "dna_r10_260bps", "dna_r10_400bps"], default="rna_r9", help='Pore generation used to sequence the data')
    return parser.parse_args()

def listener(q : mp.Queue, outfile : str):
    '''listens for messages on the q, writes to file. '''
    with open(outfile, 'w') as f:
        f.write('readid,start,end,basepos,base,motif,state,polish\n')
        i = 0
        while 1:
            m = q.get()
            i+=1
            if i%100==0:
                print(f"Segmented {i} reads", end='\r')
            if m == 'kill':
                break
            f.write(m)
            f.flush()

# def segment(rawdatapath : str, fastxpath : str, polya : dict, batch_size : int, mode : str, outfile : str, modelpath : str) -> None:
def segment(rawdatapath : str, fastxpath : str, processes : int, mode : str, outfile : str, modelpath : str, pore : str) -> None:
    files = getFiles(rawdatapath, True)
    print(f'ONT Files: {len(files)}')
    basecalls = loadFastx(fastxpath)
    print(f'{len(basecalls)} reads found')
    
    if mode == 'basic':
        CPP_SCRIPT = join(dirname(__file__), 'dynamont_NT')
    elif mode == 'banded':
        CPP_SCRIPT = join(dirname(__file__), 'dynamont_NT_banded')
    elif mode == 'resquiggle':
        CPP_SCRIPT = join(dirname(__file__), 'dynamont_NTK')

    if name == 'nt': # check for windows
        CPP_SCRIPT+='.exe'

    if processes is None:
        processes = 2 # mp.cpu_count()-1
    print(f"Using {processes} processes in segmentation.")
    pool = mp.Pool(processes)
    queue = mp.Manager().Queue()
    watcher = pool.apply_async(listener, (queue, outfile))
    
    jobs = []
    for file in files:
        r5 = read5.read(file)
        for readid in r5.getReads():
            if not readid in basecalls: # or (readid not in polyAIndex)
                continue
            jobs.append(pool.apply_async(feedSegmentationAsynchronous, (CPP_SCRIPT, {'m': modelpath, 'r' : pore}, getSignal(r5, readid), basecalls[readid][::-1], readid, queue)))
    
    # wait for all jobs to finish
    for job in jobs:
        job.get()
        
    # tell queue to terminate
    queue.put("kill")

    # # terminate pipes
    # jobs = []
    # for _ in range(processes*2):
    #     jobs.append(pool.apply_async(feedSegmentationAsynchronous, (CPP_SCRIPT, {'m': modelpath}, TERM_STRING, TERM_STRING, "End", queue)))

    # wait for all processes to finish
    for job in jobs:
        job.get()
    watcher.get()
    pool.close()
    pool.join()
    print("Done segmenting signals")

def getSignal(r5 : read5.AbstractFileReader.AbstractFileReader, readid : str):
    return r5.getZNormSignal(readid, "mean")

def main() -> None:
    args = parse()
    outfile = args.outfile
    if not exists(dirname(outfile)):
        makedirs(dirname(outfile))
    # polya=readPolyAStartEnd(args.polya)
    segment(args.raw, args.fastx, args.processes, args.mode, outfile, args.model, args.pore)

if __name__ == '__main__':
    main()