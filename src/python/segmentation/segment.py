#!/usr/bin/env python
# author: Jannes Spangenberg
# e-mail: jannes.spangenberg@uni-jena.de
# github: https://github.com/JannesSP
# website: https://jannessp.github.io

import read5_ont.AbstractFileReader
import read5_ont
import multiprocessing as mp
import pysam
import psutil
from argparse import ArgumentDefaultsHelpFormatter, ArgumentParser, Namespace
from os.path import exists, join, dirname
from os import makedirs, name, getpid
from src.python.segmentation.FileIO import feedSegmentationAsynchronous, hampelFilter

def get_memory_usage():
    """
    Returns the memory usage of the current Python process in MB.
    """
    process = psutil.Process(getpid())  # Get the current process
    memory_info = process.memory_info()   # Get memory usage details
    return memory_info.rss / 1024**2      # Convert from bytes to MB

def parse() -> Namespace:
    """
    Parse command line arguments for segmentation.

    Returns:
        Namespace: Containing the specified command line arguments
    """
    parser = ArgumentParser(
        formatter_class=ArgumentDefaultsHelpFormatter
    )
    parser.add_argument('-r', '--raw',   type=str, required=True, metavar="PATH", help='Path to raw ONT data. [POD5|FAST5]')
    parser.add_argument('-b', '--basecalls', type=str, required=True, metavar="BAM", help='Basecalls of ONT training data as .bam file')
    parser.add_argument('-o', '--outfile', type=str, required=True, help='Outpath to write files')
    parser.add_argument('--mode',  type=str, required=True, choices=['basic', 'resquiggle'], help='Segmentation algorithm used for segmentation')
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
    print(f"[{'Segmented':>9} {'Errors':>8} {'Queued':>8} {'Writer memory usage (MB)':>20}]")
    # print(f"[written,\terrors,\tin queue,\tmemory]")
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
            
            print(f"[{i:>9} {e:>8} {q.qsize():>8} {get_memory_usage():>20.0f}]", end='\r')
            # print(f"[{i},\t{e},\t{q.qsize()},\t{get_memory_usage():.2f} MB]\t", end='\r')
    
    print(f"\nReads segmented: {i}", f"Errors: {e}")
    

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
    r5 = read5_ont.read(rawFile)
    if pore in ["dna_r9", "rna_r9"]:
        # for r9 pores, shift and scale are stored for pA signal in bam
        signal = r5.getpASignal(signalid)[start:end]
    else:
        # for new pores, shift and scale is directly applied to stored integer signal (DACs)
        # this way the conversion from DACs to pA is skipped
        signal = r5.getSignal(signalid)[start:end]
    r5.close()
    signal = (signal - shift) / scale
    hampelFilter(signal)
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
    
    # directly free memory
    del r5
    del signal
    del read
    del script
    del readid
    del signalid

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
    processes = max(1, processes - 1)
    print(f"Starting segmentation with {processes} processes.")
    writer = mp.Pool(1)
    pool = mp.Pool(processes)
    queue = mp.Manager().Queue()
    writer.apply_async(listener, (queue, outfile))
    qualSkipped = 0

    with pysam.AlignmentFile(basecalls, "r" if basecalls.endswith('.sam') else "rb", check_sq=False) as samfile:
        for basecalled_read in samfile.fetch(until_eof=True):
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

            pool.apply_async(
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

    # wait for jobs to finish
    pool.close()
    pool.join()
        
    # tell queue to terminate
    queue.put("kill")
    
    # Close the pool and wait for processes to finish
    writer.close()
    writer.join()
    print(f"Skipped reads: low quality: {qualSkipped}")

def main() -> None:
    args = parse()

    outfile = args.outfile
    if not exists(dirname(outfile)):
        makedirs(dirname(outfile))

    match args.mode:
        case "basic":
            SCRIPT = 'dynamont-NT'
        case "resquiggle":
            SCRIPT = 'dynamont-NTC'

    if name == 'nt': # check for windows
        SCRIPT+='.exe'

    segment(args.raw, args.basecalls, args.processes, SCRIPT, outfile, args.model_path, args.pore, args.qscore)

if __name__ == '__main__':
    main()