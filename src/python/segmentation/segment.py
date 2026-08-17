#!/usr/bin/env python
# author: Jannes Spangenberg
# github: https://github.com/JannesSP
# website: https://jannessp.github.io

import multiprocessing as mp
import pysam
import signal
import sys
import zstandard as zstd
from argparse import ArgumentDefaultsHelpFormatter, ArgumentParser, Namespace
from collections import OrderedDict
from multiprocessing.queues import Queue
from os.path import exists, join, dirname, splitext, basename, isdir
from os import makedirs
from python import Aligner, __version__
from python.pod5_io import get_signal, open_pod5
from python.segmentation.utils import (
    segmentation_to_string,
    hampel,
    get_model,
)
from tqdm import tqdm
import numpy as np

# Globals that exist separately inside every worker process
QUEUE = None
IS_RNA = None
RAW_CACHE = None
RAW_CACHE_SIZE = None
KMERSIZE = None
ALIGNER = None

def init_worker(model_path, pore, queue, is_rna, kmer_size, mode, threads):
    global QUEUE, IS_RNA, RAW_CACHE, KMERSIZE, RAW_CACHE_SIZE, ALIGNER

    # Let the parent process handle Ctrl+C
    signal.signal(signal.SIGINT, signal.SIG_IGN)

    QUEUE = queue
    IS_RNA = is_rna
    KMERSIZE = kmer_size
    RAW_CACHE = OrderedDict()
    RAW_CACHE_SIZE = 3 # the pod5 files should more or less be ordered
    ALIGNER = Aligner(model_path, pore, mode=mode, threads=threads, band=400)

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
    parser.add_argument('-r', '--raw',   type=str, required=True, metavar="PATH", help='Path to raw ONT data. [POD5]')
    parser.add_argument('-b', '--basecalls', type=str, required=True, metavar="BAM", help='Basecalls of ONT training data as .bam file')
    parser.add_argument('-o', '--outfile', type=str, required=True, help='Path to output file. Will be zstd level 3 compressed. If directory is given, will write to dynamont.csv in that directory.')
    parser.add_argument('--mode',  type=str, required=True, choices=['basic', 'resquiggle'], help='Segmentation algorithm used for segmentation')
    parser.add_argument('--processes', type=int, default=1, help='Number of segmentation processes to use in parallel for segmentation (at least 3).\nPlease be sure that the number of processes times the number of threads does not exceed the number of available CPU cores!')
    parser.add_argument('-p', '--pore',  type=str, required=True, choices=["rna002", "rna004", "dna_r10_260bps", "dna_r10_400bps"], help='Pore generation used to sequence the data') # "dna_r9", 
    parser.add_argument('--model_path', type=str, help='Which kmer model to use for segmentation')
    parser.add_argument('-q', '--qscore', type=float, default=0.0, help='Minimal allowed quality score')
    parser.add_argument('--version', action='version', version=f'%(prog)s {__version__}')
    return parser.parse_args()

def listener(queue: Queue, outfile: str) -> None:
    """
    Listens to a queue and writes segmentation results to a zstd compressed file.
    """
    errfile = splitext(splitext(outfile)[0])[0] + ".errors"
    compressor = zstd.ZstdCompressor(level=3)
    num_err = 0

    with open(outfile, "wb") as raw:
        with compressor.stream_writer(raw) as output:
            # write header
            output.write(b'readid,signalid,start,end,basepos,base,motif,state,posterior_probability,polish\n')

            with tqdm(
                desc="Segmented reads",
                unit=" reads",
                dynamic_ncols=False,
                mininterval=10,     # update at most every 5 seconds
                file=sys.stderr,
            ) as pbar:

                while True:
                    result = queue.get()

                    if result == "kill":
                        break

                    elif isinstance(result, str): # and result.startswith("error"):
                        with open(errfile, "a") as err:
                            err.write(result + '\n')
                        num_err += 1
                        pbar.set_postfix(errors=num_err)
                        pbar.update(1)

                    else:
                        output.write(result)
                        pbar.update(1)

    print("Done segmenting reads.", file=sys.stderr)

def close_raw_cache():
    global RAW_CACHE

    if RAW_CACHE is None:
        return

    while RAW_CACHE:
        _, reader = RAW_CACHE.popitem(last=False)
        try:
            reader.close()
        except Exception:
            pass


def get_raw(path):
    global RAW_CACHE

    if path in RAW_CACHE:
        RAW_CACHE.move_to_end(path)
        return RAW_CACHE[path]

    if len(RAW_CACHE) >= RAW_CACHE_SIZE:
        _, old = RAW_CACHE.popitem(last=False)
        try:
            old.close()
        except Exception:
            pass

    RAW_CACHE[path] = open_pod5(path)

    return RAW_CACHE[path]

def _segment_one(args) -> None:
    """Segment one read with the worker's persistent native Aligner."""
    rawFile, shift, scale, start, end, read, readid, signalid = args
    r5 = get_raw(rawFile)

    signal = np.array(
        get_signal(r5, signalid, calibrated=shift <= 400)[start:end],
        dtype=np.float64,
        copy=True,
    )
    signal -= shift
    signal /= scale
    hampel(signal)

    if IS_RNA:
        read = read[::-1]
        if not read.startswith("AAAAAAAAA"):
            read = "AAAAAAAAA" + read # artificially add polyA tail to the read if it is not present, sometimes dorado cuts out the polyA tail, but the signal still contains it, which leads to a misalignment of the read and signal

    try:
        result = ALIGNER.align(signal, read, calc_probabilities=True)
        QUEUE.put(segmentation_to_string(
            result,
            readid,
            signalid,
            start,
            len(signal) + start,
            read,
            KMERSIZE,
            IS_RNA,
        ))
    except Exception as error:
        QUEUE.put(
            f"error: native, {error}\tT: {len(signal)}\tN: {len(read)}"
            f"\tRid: {readid}\tSid: {signalid}"
        )

def segmentation(args) -> None:
    """Run one segmentation job and report failures to the listener."""
    try:
        _segment_one(args)
    except Exception as error:
        _, _, _, _, _, read, readid, signalid = args
        QUEUE.put(
            f"error: worker, {error}\tN: {len(read)}"
            f"\tRid: {readid}\tSid: {signalid}"
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
    

def segment(data_path : str, basecalls : str, processes : int, outfile : str, model_path : str, pore : str, mode : str, minq : float = 0) -> None:
    """
    Segment a set of reads using a C++ script in parallel.

    Parameters
    ----------
    data_path : str
        Path to the directory containing raw ONT data.
    basecalls : str
        Path to the basecalls file in BAM format.
    processes : int
        Number of processes to use for segmentation.
    threads : int
        Number of threads to use for the native aligner.
    outfile : str
        Path to write the segmentation results to.
    model_path : str
        Path to the kmer model file used for segmentation.
        Number of processes to use for segmentation.
    outfile : str
        Path to write the segmentation results to.
    model_path : str
        Path to the kmer model file used for segmentation.
    pore : str
        Pore generation used, affects signal processing direction.
    minq : float, optional
        If set, reads with a quality score below this threshold will be skipped.

    Returns
    -------
    None
    """
    processes = max(1, processes)
    print(f"Starting segmentation with {processes} processes.", file=sys.stderr)
    # Keep the result queue unbounded so a failed writer cannot block workers.
    queue = mp.Queue()
    
    writer = mp.Process(target=listener, args=(queue, outfile))
    writer.start()

    threads = 1 # currently open mp is disabled
    processes = (processes - 2) // threads
    print(f"{threads} threads and {processes} segmentation processes", file=sys.stderr)

    pool = mp.Pool(
        processes,
        initializer=init_worker,
        initargs=(
            model_path,
            pore,
            queue,
            "rna" in pore,
            5 if pore in ["dna_r9", "rna002"] else 9,
            mode,
            threads, # number of threads used by openmp for loop parallelization
        )
    )

    # The try block runs the worker pool and consumes every generated job.
    try:
        for _ in pool.imap_unordered(
            segmentation,
            generate_jobs(data_path, basecalls, minq),
            chunksize=128
        ):
            pass

        print("Done with segmentation.", file=sys.stderr, flush=True)
    except KeyboardInterrupt:
        # Reached when the user presses Ctrl+C in the parent process.
        print("\nKeyboardInterrupt detected! Terminating...", file=sys.stderr)

        pool.terminate()
        writer.terminate()

        pool.join()
        writer.join()

        queue.close()
        queue.join_thread()

        return

    except Exception as error:
        # Reached when job generation or pool iteration raises an exception.
        print(f"Pool exception: {error!r}", file=sys.stderr, flush=True)
        pool.terminate()
        writer.terminate()

        pool.join()
        writer.join()

        raise

    else:
        # Reached only when the try block completes without an exception.
        print("Closing worker pool.", file=sys.stderr, flush=True)
        pool.close()
        pool.join()
        if writer.exitcode not in (None, 0):
            raise RuntimeError(
                f"Output writer failed with exit code {writer.exitcode}"
            )
        queue.put("kill")
        writer.join()

    finally:
        # Always reached, after normal completion, interruption, or exception.
        print("Closing queue and raw-file cache.", file=sys.stderr, flush=True)
        close_raw_cache()
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

    if args.model_path:
        model_path = args.model_path
        assert exists(model_path), "Model path does not exist"
    else:
        model_path = get_model(args.pore)
        assert exists(model_path), f"Default model not found for pore: {args.pore}, {model_path}"
    print(f"Loaded model: {basename(model_path)}", file=sys.stderr)

    segment(args.raw, args.basecalls, args.processes, outfile, model_path, args.pore, args.mode, args.qscore)

if __name__ == '__main__':
    mp.set_start_method("fork")
    main()