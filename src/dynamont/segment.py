#!/usr/bin/env python
# author: Jannes Spangenberg
# e-mail: jannes.spangenberg@uni-jena.de
# github: https://github.com/JannesSP
# website: https://jannessp.github.io

from argparse import ArgumentDefaultsHelpFormatter, ArgumentParser, Namespace
from os.path import exists, join, dirname
from os import makedirs, name

import read5.AbstractFileReader
from FileIO import feedSegmentationAsynchronous
import read5
import multiprocessing as mp
import pysam

TERM_STRING = "$"

def parse() -> Namespace:
    parser = ArgumentParser(
        formatter_class=ArgumentDefaultsHelpFormatter
    )
    parser.add_argument('-r', '--raw',   type=str, required=True, metavar="PATH", help='Path to raw ONT data. [POD5|FAST5|SLOW5]')
    parser.add_argument('-b', '--basecalls', type=str, required=True, metavar="BAM", help='Basecalls of ONT training data as .bam file')
    parser.add_argument('-o', '--outfile', type=str, required=True, help='Outpath to write files')
    parser.add_argument('--mode',  type=str, required=True, choices=['basic', 'banded', 'resquiggle'], help='Segmentation algorithm used for segmentation')
    parser.add_argument('--processes', type=int, default=None, help='Number of processes to use for segmentation')
    parser.add_argument('--model_path', type=str, help='Which kmer model to use for segmentation')
    parser.add_argument('-p', '--pore',  type=str, required=True, choices=["rna_r9", "dna_r9", "rna_rp4", "dna_r10_260bps", "dna_r10_400bps"], default="rna_r9", help='Pore generation used to sequence the data')
    parser.add_argument('-q', '--qscore', type=int, default=0, help='Minmal allowed quality score')
    return parser.parse_args()

def listener(q : mp.Queue, outfile : str):
    '''listens for messages on the q, writes to file. '''
    with open(outfile, 'w') as f:
        f.write('readid,start,end,basepos,base,motif,state,posterior_probability,polish\n')
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

def segment(dataPath : str, basecalls : str, processes : int, mode : str, outfile : str, modelpath : str, pore : str, minQual : float = None) -> None:

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
    qualSkipped = 0
    noMatchingReadid = 0
    oldFile = None
    jobs = []

    with pysam.AlignmentFile(basecalls, "r" if basecalls.endswith('.sam') else "rb", check_sq=False) as samfile:
        for basecalled_read in samfile.fetch(until_eof=True):
            
            # skip low qual reads if activated
            qs = basecalled_read.get_tag("qs")
            if minQual and qs < minQual:
                qualSkipped+=1
                continue

            # init read, sometimes a read got split by the basecaller and got a new id
            readid = basecalled_read.get_tag("pi") if basecalled_read.has_tag("pi") else basecalled_read.query_name
            seq = basecalled_read.query_sequence
            readid = basecalled_read.query_name
            ts = basecalled_read.get_tag("ts")
            sp = basecalled_read.get_tag("sp") # split start of the signal
            ns = basecalled_read.get_tag("ns") # numbers of samples used in basecalling
            rawFile = join(dataPath, basecalled_read.get_tag("fn"))

            if oldFile != rawFile:
                oldFile = rawFile
                r5 = read5.read(rawFile)

            try:
                signal = r5.getZNormSignal(readid, "median")[sp+ts:sp+ns]
            except:
                noMatchingReadid+=1
                continue

            jobs.append(pool.apply_async(feedSegmentationAsynchronous, (CPP_SCRIPT, {'m': modelpath, 'r' : pore}, signal, seq[::-1], readid, queue)))
    
    # wait for all jobs to finish
    for job in jobs:
        job.get()
        
    # tell queue to terminate
    queue.put("kill")

    # wait for all processes to finish
    for job in jobs:
        job.get()
    watcher.get()
    pool.close()
    pool.join()
    
    print("Done segmenting signals")
    print(f"Skipped reads: low quality: {qualSkipped}\treadid not found (possible split read): {noMatchingReadid}")

def main() -> None:
    args = parse()

    outfile = args.outfile
    if not exists(dirname(outfile)):
        makedirs(dirname(outfile))

    segment(args.raw, args.basecalls, args.processes, args.mode, outfile, args.model_path, args.pore, args.qscore)

if __name__ == '__main__':
    main()