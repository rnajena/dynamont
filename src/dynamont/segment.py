#!/usr/bin/env python
# author: Jannes Spangenberg
# e-mail: jannes.spangenberg@uni-jena.de
# github: https://github.com/JannesSP
# website: https://jannessp.github.io

from argparse import ArgumentDefaultsHelpFormatter, ArgumentParser, Namespace
from os.path import exists, join, dirname
from os import makedirs, name
import read5.AbstractFileReader
from FileIO import feedSegmentationAsynchronous, hampel_filter
import read5
import multiprocessing as mp
import pysam
from numpy import float32

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
    parser.add_argument('--model_path', type=str, required=True, help='Which kmer model to use for segmentation')
    parser.add_argument('-p', '--pore',  type=str, required=True, choices=["rna_r9", "dna_r9", "rna_rp4", "dna_r10_260bps", "dna_r10_400bps"], default="rna_r9", help='Pore generation used to sequence the data')
    parser.add_argument('-q', '--qscore', type=int, default=6, help='Minimal allowed quality score')
    return parser.parse_args()

def listener(q : mp.Queue, outfile : str):
    '''listens for messages on the q, writes to file. '''
    with open(outfile, 'w') as f:
        f.write('readid,signalid,start,end,basepos,base,motif,state,posterior_probability,polish\n')

        i = 0

        while 1:
            m = q.get()

            if m == 'kill':
                break

            i+=1

            # if i%100==0:
            print(f"Segmented {i} reads", end='\r')

            f.write(m)

    print(f'Segmented {i} reads\n')

def asyncSegmentation(q : mp.Queue, script : str, modelpath : str, pore : str, rawFile : str, shift : float, scale : float, start : int, end : int, read : str, readid : str, signalid : str) -> None:
    
    r5 = read5.read(rawFile)
    signal = r5.getpASignal(readid)[start:end]
    signal = (signal - shift) / scale
    signal = hampel_filter(signal)
    
    feedSegmentationAsynchronous(
                script,
                {'m': modelpath, 'r' : pore},
                signal,
                read[::-1],
                start,
                readid,
                signalid,
                q
                )

    r5.close()

def segment(dataPath : str, basecalls : str, processes : int, CPP_SCRIPT : str, outfile : str, modelpath : str, pore : str, minQual : float = None) -> None:

    if processes is None:
        processes = 2 # mp.cpu_count()-1
    print(f"Using {processes} processes in segmentation.")
    pool = mp.Pool(processes)
    queue = mp.Manager().Queue()
    watcher = pool.apply_async(listener, (queue, outfile))
    qualSkipped = 0
    # noMatchingReadid = 0
    jobs = [None for _ in range(processes)]

    with pysam.AlignmentFile(basecalls, "r" if basecalls.endswith('.sam') else "rb", check_sq=False) as samfile:
        for bi, basecalled_read in enumerate(samfile.fetch(until_eof=True)):
            
            # skip low qual reads if activated
            qs = basecalled_read.get_tag("qs")
            if minQual and qs < minQual:
                qualSkipped+=1
                continue

            # init read
            readid = basecalled_read.query_name
            # if read got split by basecaller, another readid is assign, pi holds the read id from the pod5 file
            signalid = basecalled_read.get_tag("pi") if basecalled_read.has_tag("pi") else readid
            seq = basecalled_read.query_sequence # change direction from 5' - 3' to 3' - 5'
            ts = basecalled_read.get_tag("ts")
            ns = basecalled_read.get_tag("ns") # numbers of samples used in basecalling for this readid
            sp = basecalled_read.get_tag("sp") if basecalled_read.has_tag("sp") else 0 # if split read get start offset of the signal
            rawFile = join(dataPath, basecalled_read.get_tag("fn"))
            shift = basecalled_read.get_tag("sm")
            scale = basecalled_read.get_tag("sd")

            jobs[bi % processes] = pool.apply_async(
                asyncSegmentation, (
                    queue,
                    CPP_SCRIPT,
                    modelpath,
                    pore,
                    rawFile,
                    shift,
                    scale,
                    sp+ts,
                    sp+ns,
                    seq,
                    readid,
                    signalid
                    )
                )

    # wait for last job batch to finish
    for job in jobs:
        job.get()
        
    # tell queue to terminate
    queue.put("kill")
    watcher.get()
    pool.close()
    pool.join()
    
    print(f"Skipped reads: low quality: {qualSkipped}")

def main() -> None:
    args = parse()

    outfile = args.outfile
    if not exists(dirname(outfile)):
        makedirs(dirname(outfile))

    match args.mode:
        case "basic":
            CPP_SCRIPT = join(dirname(__file__), 'dynamont_NT')
        case "banded":
            CPP_SCRIPT = join(dirname(__file__), 'dynamont_NT_banded')
        case "resquiggle":
            CPP_SCRIPT = join(dirname(__file__), 'dynamont_NTK')

    if name == 'nt': # check for windows
        CPP_SCRIPT+='.exe'

    segment(args.raw, args.basecalls, args.processes, CPP_SCRIPT, outfile, args.model_path, args.pore, args.qscore)

if __name__ == '__main__':
    main()