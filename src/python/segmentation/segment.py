#!/usr/bin/env python
# author: Jannes Spangenberg
# github: https://github.com/JannesSP
# website: https://jannessp.github.io

import multiprocessing as mp
import pysam
import read5_ont
import signal
import sys
import zstandard as zstd
from argparse import ArgumentDefaultsHelpFormatter, ArgumentParser, Namespace
from collections import OrderedDict
from multiprocessing.queues import Queue
from os.path import exists, join, dirname, splitext, basename, isdir
from os import makedirs, name
from python.segmentation.FileIO import feedSegmentationAsynchronous, hampelFilter, getModel
from python._version import __version__

# Globals that exist separately inside every worker process
SCRIPT = None
MODELPATH = None
PORE = None
QUEUE = None
OPTIONS = None
IS_RNA = None
RAW_CACHE = None
RAW_CACHE_SIZE = None
KMERSIZE = None

def init_worker(script, modelpath, pore, queue, is_rna, kmer_size):
    global SCRIPT, MODELPATH, PORE, QUEUE, OPTIONS, IS_RNA, RAW_CACHE, KMERSIZE, RAW_CACHE_SIZE

    # Let the parent process handle Ctrl+C
    signal.signal(signal.SIGINT, signal.SIG_IGN)

    SCRIPT = script
    MODELPATH = modelpath
    PORE = pore
    QUEUE = queue
    IS_RNA = is_rna
    KMERSIZE = kmer_size
    RAW_CACHE = OrderedDict()
    RAW_CACHE_SIZE = 4

    OPTIONS = {
        'm': modelpath,
        'r': pore,
        't': 4
    }

def parse() -> Namespace:
    """
    Parse command line arguments for segmentation.

    Returns:
        Namespace: Containing the specified command line arguments
    """
    parser = ArgumentParser(
        formatter_class=ArgumentDefaultsHelpFormatter,
        prog='dynamont-resquiggle',
    )
    parser.add_argument('-r', '--raw',   type=str, required=True, metavar="PATH", help='Path to raw ONT data. [POD5|FAST5]')
    parser.add_argument('-b', '--basecalls', type=str, required=True, metavar="BAM", help='Basecalls of ONT training data as .bam file')
    parser.add_argument('-o', '--outfile', type=str, required=True, help='Path to output file. Will be zstd level 3 compressed. If directory is given, will write to dynamont.csv in that directory.')
    parser.add_argument('--mode',  type=str, required=True, choices=['basic', 'resquiggle'], help='Segmentation algorithm used for segmentation')
    parser.add_argument('--processes', type=int, default=mp.cpu_count()-1, help='Number of processes to use for segmentation')
    parser.add_argument('-p', '--pore',  type=str, required=True, choices=["rna002", "rna004", "dna_r10_260bps", "dna_r10_400bps"], help='Pore generation used to sequence the data') # "dna_r9", 
    parser.add_argument('--model_path', type=str, help='Which kmer model to use for segmentation')
    parser.add_argument('-q', '--qscore', type=float, default=0.0, help='Minimal allowed quality score')
    parser.add_argument('--version', action='version', version=f'%(prog)s {__version__}')
    return parser.parse_args()

def listener(queue: Queue, outfile: str) -> None:
    """
    Listens to a queue and writes segmentation results to a zstd compressed file.
    """
    print(f"{'Segmented':>9} | {'Errors':>8}", file=sys.stderr)
    errfile = splitext(splitext(outfile)[0])[0] + ".errors"
    compressor = zstd.ZstdCompressor(level=3)
    num_reads = 0
    num_err = 0

    with open(outfile, "wb") as raw:
        with compressor.stream_writer(raw) as output:
            # write header
            output.write(b'readid,signalid,start,end,basepos,base,motif,state,posterior_probability,polish\n')

            while True:
                result = queue.get() # this is already in bytes

                if result == "kill":
                    break

                elif isinstance(result, str): # and result.startswith("error"):
                    with open(errfile, "a") as err:
                        err.write(f"{result}\n")
                    num_err += 1

                else:
                    num_reads += 1
                    output.write(result)

                if num_reads%10 == 0:
                    print(f"{num_reads:>9} | {num_err:>8}", end="\r", file=sys.stderr)
    print(f"\nReads segmented: {num_reads}", f"Errors: {num_err}", file=sys.stderr)

def get_raw(path):
    global RAW_CACHE

    if path in RAW_CACHE:
        RAW_CACHE.move_to_end(path)
        return RAW_CACHE[path]

    if len(RAW_CACHE) >= RAW_CACHE_SIZE:
        _, old = RAW_CACHE.popitem(last=False)
        old.close()

    RAW_CACHE[path] = read5_ont.read(path)

    return RAW_CACHE[path]

def asyncSegmentation(args) -> None:
    """
    Asynchronously segments a raw signal using a C++ script and places the results in a queue.

    Parameters
    ----------
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
    rawFile, shift, scale, start, end, read, readid, signalid = args
    r5 = get_raw(rawFile)

    #! I do not know anymore in which version, but in some 0.9.x dorado version, the shift and scale values were taken from the raw DACS values instead of the pA signal
    if shift > 400:
        signal = r5.getSignal(signalid)[start:end]
    else:
        signal = r5.getpASignal(signalid)[start:end]

    # shift = np.median(signal[start:end])
    # scale = np.median(np.abs(signal[start:end] - shift))

    # standardize signal
    signal -= shift
    signal /= scale
    hampelFilter(signal)

    if IS_RNA:
        read = read[::-1] # change direction from 5' - 3' to 3' - 5'
        if not read.startswith("AAAAAAAAA"):
            read = "AAAAAAAAA" + read
    feedSegmentationAsynchronous(
                SCRIPT,
                OPTIONS,
                signal,
                read,
                start,
                readid,
                signalid,
                KMERSIZE,
                QUEUE,
                IS_RNA
                )

def generate_jobs(dataPath : str, basecalls : str, minQual : float = 0):
    """
    Generate segmentation jobs from a BAM/SAM file.

    Yields tuples matching asyncSegmentation() arguments.

    Parameters
    ----------
    dataPath : str
        Path containing raw ONT files.
    basecalls : str
        Path to BAM/SAM basecalls.
    minQual : float
        Minimum allowed quality score.
    queue : Queue, optional
        Queue used for reporting missing raw files.

    Yields
    ------
    tuple
        (
            rawFile,
            shift,
            scale,
            start,
            end,
            sequence,
            readid,
            signalid
        )
    """
    qualSkipped = 0

    with pysam.AlignmentFile(basecalls, "rb", check_sq=False) as samfile:
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
            rawFile = join(dataPath, basecalled_read.get_tag("fn")) if basecalled_read.has_tag("fn") else join(dataPath, basecalled_read.get_tag("f5"))
            # if not exists(rawFile): 
            #     queue.put(f"error: 6, no raw file found\t{rawFile}\t{readid}\t{signalid}")
            #     continue

            #! normalize whole signal
            shift = basecalled_read.get_tag("sm")
            scale = basecalled_read.get_tag("sd")
            
            yield (
                rawFile,
                shift,
                scale,
                sp+ts,
                sp+ns,
                seq,
                readid,
                signalid
            )

    print(f"Skipped reads due to low quality: {qualSkipped}", file=sys.stderr)
    

def segment(dataPath : str, basecalls : str, processes : int, script : str, outfile : str, modelpath : str, pore : str, minQual : float = 0) -> None:
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
    script : str
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
    processes = max(1, processes // 2) # half of the processes will spawn the cpp segmentation process
    print(f"Starting segmentation with {processes} processes.", file=sys.stderr)
    queue = mp.Queue(maxsize=1000) # workers naturally backpressure
    # queue = mp.SimpleQueue()
    
    writer = mp.Process(target=listener, args=(queue, outfile))
    writer.start()

    pool = mp.Pool(
        processes,
        initializer=init_worker,
        initargs=(
            script,
            modelpath,
            pore,
            queue,
            "rna" in pore,
            5 if pore in ["dna_r9", "rna002"] else 9
        )
    )

    try:
        for _ in pool.imap_unordered(
            asyncSegmentation,
            generate_jobs(dataPath, basecalls, minQual),
            chunksize=16
        ):
            pass

    except KeyboardInterrupt:
        print("\nKeyboardInterrupt detected! Terminating...", file=sys.stderr)

        pool.terminate()
        writer.terminate()

        pool.join()
        writer.join()

        queue.close()
        queue.join_thread()

        return

    except Exception:
        pool.terminate()
        writer.terminate()

        pool.join()
        writer.join()

        raise

    else:
        pool.close()
        pool.join()

        queue.put("kill")
        writer.join()

    finally:
        queue.close()
        queue.join_thread()

def main() -> None:
    args = parse()

    outfile = args.outfile

    if isdir(outfile):
        outfile = join(outfile, "dynamont.csv.zst")
    elif not outfile.endswith(".zst"):
        outfile += ".zst"
    parent = dirname(outfile)
    if parent and not exists(parent):
        makedirs(parent)

    match args.mode:
        case "basic":
            # SCRIPT = 'dynamont-NT'
            CPP_SCRIPT = 'dynamont-NT-banded'
        case "resquiggle":
            CPP_SCRIPT = 'dynamont-NTC'

    if name == 'nt': # check for windows
        CPP_SCRIPT+='.exe'

    if args.model_path:
        model_path = args.model_path
        assert exists(model_path), "Model path does not exist"
    else:
        model_path = getModel(args.pore)
        assert exists(model_path), f"Default model not found for pore: {args.pore}, {model_path}"
    print(f"Loaded model: {basename(model_path)}", file=sys.stderr)

    segment(args.raw, args.basecalls, args.processes, CPP_SCRIPT, outfile, model_path, args.pore, args.qscore)

if __name__ == '__main__':
    mp.set_start_method("fork")
    main()