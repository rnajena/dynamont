#!/usr/bin/env python
# author: Jannes Spangenberg
# e-mail: jannes.spangenberg@uni-jena.de
# github: https://github.com/JannesSP
# website: https://jannessp.github.io

# TODO make this asynchroneous

from argparse import ArgumentDefaultsHelpFormatter, ArgumentParser, Namespace
from os.path import exists, join, dirname
from os import makedirs, name
from FileIO import getFiles, loadFastx, openCPPScript, feedSegmentationAsynchronous
from read5 import read
import multiprocessing as mp

def parse() -> Namespace:
    parser = ArgumentParser(
        formatter_class=ArgumentDefaultsHelpFormatter
    )
    parser.add_argument('--raw', type=str, required=True, help='Raw ONT training data')
    parser.add_argument('--fastx', type=str, required=True, help='Basecalls of ONT training data')
    # parser.add_argument('--polya', type=str, required=True, help='Poly A table from nanopolish polya containing the transcript starts')
    parser.add_argument('--out', type=str, required=True, help='Outpath to write files')
    parser.add_argument('--batch_size', type=int, default=16, help='Number of reads to train before updating')
    parser.add_argument('--model_path', type=str, default=None)
    parser.add_argument('--mode', choices=['basic', 'indel'], required=True)
    return parser.parse_args()

def listener(q : mp.Queue, outfile : str):
    '''listens for messages on the q, writes to file. '''
    with open(outfile, 'w') as f:
        f.write('readid,start,end,basepos,base,motif,state\n')
        i = 0
        while 1:
            m = q.get()
            i+=1
            if i%1000==0:
                print(f"Segmented {i} reads", end='\r')
            if m == 'kill':
                break
            f.write(m)
            f.flush()

# def segment(rawdatapath : str, fastxpath : str, polya : dict, batch_size : int, mode : str, outfile : str, modelpath : str) -> None:
def segment(rawdatapath : str, fastxpath : str, batch_size : int, mode : str, outfile : str, modelpath : str) -> None:
    files = getFiles(rawdatapath, True)
    print(f'ONT Files: {len(files)}')
    basecalls = loadFastx(fastxpath)
    print(f'{len(basecalls)} reads found')
    
    if mode == 'indel':
        CPP_SCRIPT = join(dirname(__file__), 'segmentation_indel')
        if name == 'nt': # check for windows
            CPP_SCRIPT+='.exe'
    elif mode == 'basic':
        CPP_SCRIPT = join(dirname(__file__), 'segmentation_basic')
        if name == 'nt': # check for windows
            CPP_SCRIPT+='.exe'
    
    i = 0
    pool = mp.Pool(batch_size)
    pipes = [openCPPScript(CPP_SCRIPT) for _ in range(batch_size - 1)]
    queue = mp.Manager().Queue()
    pool.apply_async(listener, (queue, outfile))

    jobs = []
    for file in files:
        r5 = read(file)
        for readid in r5:
            if not readid in basecalls: # or (readid not in polyAIndex)
                continue
            job = pool.apply_async(feedSegmentationAsynchronous, (r5.getpASignal(readid, "mean"), basecalls[readid][::-1], readid, pipes[i%len(pipes)], queue))
            jobs.append(job)
            i+=1
            
    for job in jobs:
        job.get()
    
    queue.put("kill")
    pool.close()
    pool.join()

def main() -> None:
    args = parse()
    outdir = args.out
    if not exists(outdir):
        makedirs(outdir)
    # polya=readPolyAStartEnd(args.polya)
    # segment(args.raw, args.fastx, polya, args.batch_size, args.mode, join(outdir, "dynamont.csv"), args.model_path)
    segment(args.raw, args.fastx, args.batch_size, args.mode, join(outdir, "dynamont.csv"), args.model_path)

if __name__ == '__main__':
    main()