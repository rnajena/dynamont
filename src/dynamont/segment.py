#!/usr/bin/env python
# author: Jannes Spangenberg
# e-mail: jannes.spangenberg@uni-jena.de
# github: https://github.com/JannesSP
# website: https://jannessp.github.io

from argparse import ArgumentDefaultsHelpFormatter, ArgumentParser, Namespace
from os.path import exists, join, dirname
from os import makedirs, name
import read5.AbstractFileReader
from FileIO import feedSegmentationAsynchronous, hampelFilter
import read5
import multiprocessing as mp
import pysam

def parse() -> Namespace:
    """
    Parse command line arguments for segmentation.

    Returns:
        Namespace: Containing the specified command line arguments
    """
    parser = ArgumentParser(
        formatter_class=ArgumentDefaultsHelpFormatter
    )
    parser.add_argument('-r', '--raw',   type=str, required=True, metavar="PATH", help='Path to raw ONT data. [POD5|FAST5|SLOW5]')
    parser.add_argument('-b', '--basecalls', type=str, required=True, metavar="BAM", help='Basecalls of ONT training data as .bam file')
    parser.add_argument('-o', '--outfile', type=str, required=True, help='Outpath to write files')
    parser.add_argument('--mode',  type=str, required=True, choices=['basic', 'banded', 'resquiggle'], help='Segmentation algorithm used for segmentation')
    parser.add_argument('--processes', type=int, default=mp.cpu_count()-1, help='Number of processes to use for segmentation')
    parser.add_argument('--model_path', type=str, required=True, help='Which kmer model to use for segmentation')
    parser.add_argument('-p', '--pore',  type=str, required=True, choices=["rna_r9", "dna_r9", "rna_rp4", "dna_r10_260bps", "dna_r10_400bps"], help='Pore generation used to sequence the data')
    parser.add_argument('-q', '--qscore', type=float, default=0.0, help='Minimal allowed quality score')
    return parser.parse_args()

def listener(q : mp.Queue, outfile : str) -> None:
    """
    Listens to a queue and writes segmentation results to a file.

    Takes a mp.Queue and a file path to write the results to. The function
    will run until it receives 'kill' from the queue. If any error occurs it
    will be printed and counted.

    Parameters
    ----------
    q : mp.Queue
        The queue to listen to
    outfile : str
        The path to write the results to
    """
    with open(outfile, 'w') as f:
        f.write('readid,signalid,start,end,basepos,base,motif,state,posterior_probability,polish\n')
        i = 0
        e = 0
        while True:
            m = q.get()
            
            if m == 'kill':
                break
            elif m.startswith('error'):
                print(m)
                e+=1
            else:
                i+=1
                f.write(m)
                f.flush()
            
            print(f"Reads segmented: {i}", f"Errors: {e}", end='\r')
    
    print(f"Reads segmented: {i}", f"Errors: {e}")
    

def asyncSegmentation(q : mp.Queue, script : str, modelpath : str, pore : str, rawFile : str, shift : float, scale : float, start : int, end : int, read : str, readid : str, signalid : str) -> None:
    """
    Asynchronously segments a raw signal using a C++ script and places the results in a queue.

    Parameters
    ----------
    q : mp.Queue
        Queue to store segmentation results or errors.
    script : str
        Path to the C++ segmentation script.
    modelpath : str
        Path to the kmer model file used for segmentation.
    pore : str
        Pore generation used, affects signal processing direction.
    rawFile : str
        Path to the file containing raw signal data.
    shift : float
        Signal shift value for normalization.
    scale : float
        Signal scale value for normalization.
    start : int
        Start index for the signal segment.
    end : int
        End index for the signal segment.
    read : str
        Nucleotide sequence in 5' -> 3' direction.
    readid : str
        Identifier for the read within the basecall file.
    signalid : str
        Signal identifier within the basecall file.

    Returns
    -------
    None
    """
    r5 = read5.read(rawFile)
    signal = r5.getpASignal(signalid)[start:end]
    r5.close()
    signal = (signal - shift) / scale
    signal = hampelFilter(signal)
    if "rna" in pore:
        read = read[::-1] # change direction from 5' - 3' to 3' - 5'
    
    feedSegmentationAsynchronous(
                script,
                {'m': modelpath, 'r' : pore},
                signal,
                read,
                start,
                readid,
                signalid,
                q
                )

def segment(dataPath : str, basecalls : str, processes : int, SCRIPT : str, outfile : str, modelpath : str, pore : str, minQual : float = 0) -> None:
    """
    Segment a set of reads using a C++ script in parallel.

    Parameters
    ----------
    dataPath : str
        Path to the directory containing raw ONT data.
    basecalls : str
        Path to the basecalls file in BAM format.
    processes : int
        Number of processes to use for segmentation.
    SCRIPT : str
        Path to the C++ script to use for segmentation.
    outfile : str
        Path to write the segmentation results to.
    modelpath : str
        Path to the kmer model file used for segmentation.
    pore : str
        Pore generation used, affects signal processing direction.
    minQual : float, optional
        If set, reads with a quality score below this threshold will be skipped.

    Returns
    -------
    None
    """
    processes = max(2, processes)
    print(f"Using {processes} processes in segmentation.")
    pool = mp.Pool(processes)
    queue = mp.Manager().Queue()
    pool.apply_async(listener, (queue, outfile))
    qualSkipped = 0
    jobs = [None for _ in range(processes)]

    with pysam.AlignmentFile(basecalls, "r" if basecalls.endswith('.sam') else "rb", check_sq=False) as samfile:
        print("Starting segmentation.")
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
            seq = basecalled_read.query_sequence
            ns = basecalled_read.get_tag("ns") # ns:i: 	the number of samples in the signal prior to trimming
            ts = basecalled_read.get_tag("ts") # ts:i: 	the number of samples trimmed from the start of the signal
            sp = basecalled_read.get_tag("sp") if basecalled_read.has_tag("sp") else 0 # if split read get start offset of the signal
            rawFile = join(dataPath, basecalled_read.get_tag("fn"))
            shift = basecalled_read.get_tag("sm")
            scale = basecalled_read.get_tag("sd")

            jobs[bi % processes] = pool.apply_async(
                asyncSegmentation, (
                    queue,
                    SCRIPT,
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
    # queue.put("kill")
    
    # Close the pool and wait for processes to finish
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
            SCRIPT = join(dirname(__file__), 'dynamont_NT')
        case "banded":
            SCRIPT = join(dirname(__file__), 'dynamont_NT_banded')
        case "resquiggle":
            SCRIPT = join(dirname(__file__), 'dynamont_NTK')

    if name == 'nt': # check for windows
        SCRIPT+='.exe'

    segment(args.raw, args.basecalls, args.processes, SCRIPT, outfile, args.model_path, args.pore, args.qscore)

if __name__ == '__main__':
    main()