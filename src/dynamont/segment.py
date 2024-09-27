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
    parser.add_argument('--outfile', type=str, required=True, help='Outpath to write files')
    parser.add_argument('--batch_size', type=int, default=None, help='Number of reads to train before updating')
    parser.add_argument('--model_path', type=str, default=None)
    parser.add_argument('--mode', choices=['basic', 'indel', '3d'], required=True)
    parser.add_argument('--unnormalised', action="store_true", help="Use unnormalised signal. Make sure the model matches this mode!")
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
def segment(rawdatapath : str, fastxpath : str, batch_size : int, mode : str, outfile : str, modelpath : str, unnormalised : bool) -> None:
    files = getFiles(rawdatapath, True)
    print(f'ONT Files: {len(files)}')
    basecalls = loadFastx(fastxpath)
    print(f'{len(basecalls)} reads found')
    
    if mode == '3d':
        CPP_SCRIPT = join(dirname(__file__), 'segmentation_3d_sparsed')
        if name == 'nt': # check for windows
            CPP_SCRIPT+='.exe'
    elif mode == 'basic':
        CPP_SCRIPT = join(dirname(__file__), 'segmentation_basic_sparsed')
        if name == 'nt': # check for windows
            CPP_SCRIPT+='.exe'
    elif mode == 'indel':
        CPP_SCRIPT = join(dirname(__file__), 'segmentation_indel')
        if name == 'nt': # check for windows
            CPP_SCRIPT+='.exe'

    if batch_size is None:
        batch_size = 2 # mp.cpu_count()-1
    print(f"Using {batch_size} processes in segmentation.")
    pool = mp.Pool(batch_size)
    queue = mp.Manager().Queue()
    watcher = pool.apply_async(listener, (queue, outfile))
    
    jobs = []
    for file in files:
        r5 = read5.read(file)
        for readid in r5.getReads():
            if not readid in basecalls: # or (readid not in polyAIndex)
                continue
            jobs.append(pool.apply_async(feedSegmentationAsynchronous, (CPP_SCRIPT, {'m': modelpath}, getSignal(r5, readid, unnormalised), basecalls[readid][::-1], readid, queue)))
    
    # wait for all jobs to finish
    for job in jobs:
        job.get()
        
    # tell queue to terminate
    queue.put("kill")

    # # terminate pipes
    # jobs = []
    # for _ in range(batch_size*2):
    #     jobs.append(pool.apply_async(feedSegmentationAsynchronous, (CPP_SCRIPT, {'m': modelpath}, TERM_STRING, TERM_STRING, "End", queue)))

    # wait for all processes to finish
    for job in jobs:
        job.get()
    watcher.get()
    pool.close()
    pool.join()
    print("Done segmenting signals")

def getSignal(r5 : read5.AbstractFileReader.AbstractFileReader, readid : str, unnormalised : bool):
    if unnormalised:
        return r5.getpASignal(readid)
    return r5.getZNormSignal(readid, "median")

def main() -> None:
    args = parse()
    outfile = args.outfile
    if not exists(dirname(outfile)):
        makedirs(dirname(outfile))
    # polya=readPolyAStartEnd(args.polya)
    # segment(args.raw, args.fastx, polya, args.batch_size, args.mode, join(outdir, "dynamont.csv"), args.model_path)
    segment(args.raw, args.fastx, args.batch_size, args.mode, outfile, args.model_path, args.unnormalised)

if __name__ == '__main__':
    main()