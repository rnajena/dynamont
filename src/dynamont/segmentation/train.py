#!/usr/bin/env python
# author: Jannes Spangenberg
# github: https://github.com/JannesSP
# website: https://jannessp.github.io

import numpy as np
import multiprocessing as mp
import pysam
import sys
from argparse import ArgumentDefaultsHelpFormatter, ArgumentParser, Namespace
from os.path import exists, join, dirname, basename
from os import makedirs
from datetime import datetime
from collections import deque
from dynamont.segmentation.utils import calcZ, plt_parameters, train_transition_emission, read_kmer_model, write_kmer_model, hampel, cnt_nts_ratios, get_model
from dynamont import __version__
from dynamont.pod5_io import get_signal, open_pod5

class ManagedList:
    def __init__(self, values, max_size=100):
        """Initialize the ManagedList with a maximum size."""
        self.values = deque(values, maxlen=max_size)  # deque automatically limits size

    def add(self, value):
        """Add a new value to the list, removing the oldest if necessary."""
        self.values.append(value)

    def get_list(self):
        """Returns the current list of values."""
        return list(self.values)

    def __repr__(self):
        """String representation of the ManagedList."""
        return f"ManagedList({list(self.values)})"
    
    def mean(self):
        """Returns the mean of the values in the managed list."""
        if len(self.values) == 0:
            return None  # Handle the case where the list is empty
        return np.mean(self.values)

    def median(self):
        """Returns the median of the values in the managed list."""
        if len(self.values) == 0:
            return None  # Handle the case where the list is empty
        return np.median(self.values)

def parse() -> Namespace:
    parser = ArgumentParser(
        formatter_class=ArgumentDefaultsHelpFormatter,
        prog='dynamont-train',
    )
    # required
    parser.add_argument('-r', '--raw',   type=str, required=True, metavar="PATH", help='Path to raw ONT data. [POD5|FAST5]')
    parser.add_argument('-b', '--basecalls', type=str, required=True, metavar="BAM", help='Basecalls of ONT training data as .bam file')
    parser.add_argument('-o', '--outdir',   type=str, required=True, metavar="PATH", help='Outpath to write files')
    parser.add_argument('-p', '--pore',  type=str, required=True, choices=["rna002", "rna004", "dna_r10_260bps", "dna_r10_400bps"], help='Pore generation used to sequence the data') # "dna_r9"
    # parser.add_argument('--mode',  type=str, required=True, choices=['basic', 'resquiggle'], help='Segmentation algorithm used for segmentation')
    # optional
    parser.add_argument('--model_path', type=str, help='Which initial kmer models to use for training')
    parser.add_argument('--batch_size', type=int, default=24, help='Number of reads to train before updating')
    parser.add_argument('--max_batches', type=int, default=None, help='Numbers of batches to train each epoch')
    parser.add_argument('-e', '--epochs', type=int, default=1, help='Number of training epochs')
    parser.add_argument('-q', '--qscore', type=float, default=10.0, help='Minimal allowed quality score')
    parser.add_argument('--version', action='version', version=f'%(prog)s {__version__}')
    return parser.parse_args()

def train(data_path : str, basecalls : str, batch_size : int, epochs :int, param_file : str, mode : str, model_path : str, max_batches : int, pore : str, minq : float = None) -> None:

    model = read_kmer_model(model_path)
    trained_model = join(dirname(param_file), 'trained_0_0.model')
    write_kmer_model(trained_model, model)
    
    param_writer = open(param_file, 'w')

    if mode == 'basic':
        ALIGNER_MODE = 'basic'
        transition_params = {
            'e1': 1.0,
            'm1': 0.03,
            'e2': 0.97
            }
    elif mode == 'resquiggle':
        ALIGNER_MODE = 'resquiggle'
        transition_params = {
            'a1': 0.012252440188168037,
            'a2': 0.246584724985145,
            'p1': 0.04477093133243305,
            'p2': 0.007687811003133089,
            'p3': 0.4469623669791557,
            's1': 0.05321209670114726,
            's2': 0.0007555035568187239,
            's3': 0.21999557711272136,
            'e1': 1.0,
            'e2': 0.9467879033992115,
            'e3': 0.9552290685034269,
            'e4': 0.9792321612614708,
            'i1': 7.208408117990252e-05,
            'i2': 0.08645733058947891
            }
    else:
        print(f'Mode {mode} not implemented', file=sys.stderr)
        sys.exit(1)

    # collect trained parameters to get an ensemble training in the end to reduce outlier trainings
    param_collector = {kmer:(ManagedList([model[kmer][0]]), ManagedList([model[kmer][1]])) for kmer in model}
    param_collector = param_collector | {param : ManagedList([transition_params[param]]) for param in transition_params}
    kmers_seen = set()
    param_writer.write("epoch,batch,read,")
    for param in transition_params:
        param_writer.write(param+',')
    param_writer.write("Zchange\n")
    i = 0
    qskips = 0
    readid_mismatches = 0

    with mp.Pool(batch_size) as pool:

        for e in range(epochs):
            mp_items = []
            trainIDs = []
            cbatch = 0
            oldFile = None

            with pysam.AlignmentFile(basecalls, "r" if basecalls.endswith('.sam') else "rb", check_sq=False) as samfile:
                for read in samfile.fetch(until_eof=True):
                    
                    # skip low qual reads if activated
                    qual = read.get_tag("qs")
                    if minq and qual < minq:
                        qskips+=1
                        continue

                    # init read
                    seq = read.query_sequence
                    counts = cnt_nts_ratios(seq)

                    # saw weird signals in rp4 data, a very homogenous signal jumping between two values
                    # does not look like a normal read, more like an artifact
                    # produces reads with high quality but the read consists only of repeats
                    if any(counts[base] >= 0.6 for base in counts.keys()):
                        continue

                    # init read, sometimes a read got split by the basecaller and got a new id
                    readid = read.get_tag("pi") if read.has_tag("pi") else read.query_name
                    sp = read.get_tag("sp") if read.has_tag("sp") else 0 # if split read get start offset of the signal
                    ns = read.get_tag("ns") # ns:i: 	the number of samples in the signal prior to trimming
                    ts = read.get_tag("ts") # ts:i: 	the number of samples trimmed from the start of the signal
                    ont_file = join(data_path, read.get_tag("fn"))
                    start = sp+ts
                    end = sp + ns
                    shift = read.get_tag("sm")
                    scale = read.get_tag("sd")

                    if oldFile != ont_file:
                        oldFile = ont_file
                        r5 = open_pod5(ont_file)

                    # fill batch
                    if len(mp_items) < batch_size:
                        # saw more consistency for short reads when using the mean
                        try:
                            signal = get_signal(r5, readid, calibrated=shift <= 400)[start:end]
                        except:
                            readid_mismatches+=1
                            continue

                        signal -= shift
                        signal /= scale
                        hampel(signal, 7, 5.) # small window and high variance allowed: just to filter outliers that result from sensor errors, rest of the original signal should be kept
                        if "rna" in pore:
                            seq = seq[::-1]
                            if not seq.startswith("AAAAAAAAA"):
                                seq = "AAAAAAAAA" + seq
                        mp_items.append([signal, seq, transition_params | {'r' : pore, 't' : 4}, ALIGNER_MODE, trained_model, readid])
                        trainIDs.append(readid)

                    if len(mp_items) == batch_size:
                        print("============================", file=sys.stderr)
                        print(f"{datetime.now().strftime('%Y-%m-%d_%H-%M-%S')}: Training epoch: {e}, reads: {i}, batch: {cbatch}\n{transition_params}", file=sys.stderr)
                        print("Training with read:", trainIDs, file=sys.stderr)
                        cbatch += 1
                        preZ = []

                        for result in pool.starmap(train_transition_emission, mp_items):
                            if isinstance(result, str):
                                print(f"No segmentation calculated for {result} in {e}: {trained_model}.", file=sys.stderr)
                                continue

                            trainedParams, newModels, Z = result

                            i += 1
                            preZ.append(Z)

                            for param in trainedParams:
                                param_collector[param].add(trainedParams[param])

                            #! skip weird trainings
                            if ('AAAAAAAAA' in newModels and newModels['AAAAAAAAA'][0] < 0.5) or ('AAAAA' in newModels and newModels['AAAAA'][0] < 0.5):
                                continue

                            for kmer in newModels:
                                kmers_seen.add(kmer)
                                param_collector[kmer][0].add(newModels[kmer][0])
                                param_collector[kmer][1].add(newModels[kmer][1])

                        print(f"Zs: {preZ}", file=sys.stderr)

                        # update parameters
                        param_writer.write(f'{e},{cbatch},{i},') # log
                        for param in transition_params:
                            try:
                                transition_params[param] = param_collector[param].mean()
                            except:
                                print(param, param_collector[param].get_list(), file=sys.stderr)
                                sys.exit(1)
                            param_writer.write(f'{transition_params[param]},') # log

                        for kmer in kmers_seen:
                            model[kmer] = [param_collector[kmer][0].mean(), param_collector[kmer][1].mean()]

                        trained_model = join(dirname(trained_model), f"trained_{e}_{cbatch}.model")
                        write_kmer_model(trained_model, model)
                        param_writer.flush()

                        # rerun with new parameters to compare Zs
                        for j in range(len(mp_items)):
                            mp_items[j][2] = transition_params | {'r' : pore, 't' : 4}
                            mp_items[j][4] = trained_model
                        postZ = []
                        for result in pool.starmap(calcZ, mp_items):
                            if isinstance(result, str):
                                print(f"No segmentation calculated for {result} in {e} calcZ.", file=sys.stderr)
                                continue
                            postZ.append(result)

                        dZ = np.array(postZ) - np.array(preZ)

                        print(f"Z changes: {dZ}", file=sys.stderr)
                        deltaZ = np.mean(dZ) if len(dZ) else 0.0
                        param_writer.write(f'{deltaZ}\n') # log
                        param_writer.flush() # log
                        # initialize new batch
                        kmers_seen = set()
                        mp_items = []
                        trainIDs = []

                        if max_batches is not None and cbatch >= max_batches:
                            break

        param_writer.close()
        print("Done training", file=sys.stderr)
        print(f"Skipped reads due to low quality: {qskips}", file=sys.stderr)

def main() -> None:
    args = parse()

    # init outdir
    outdir = args.outdir + f'_{datetime.now().strftime("%Y-%m-%d_%H-%M-%S")}'
    if not exists(outdir):
        makedirs(outdir)

    paramFile = join(outdir, 'params.csv')

    if args.model_path:
        model_path = args.model_path
        assert exists(model_path), "Model path does not exist"
    else:
        model_path = get_model(args.pore)
        assert exists(model_path), f"Default model not found for pore: {args.pore}, {model_path}"
    print(f"Loaded model: {basename(model_path)}", file=sys.stderr)

    train(args.raw, args.basecalls, args.batch_size, args.epochs, paramFile, "basic", model_path, args.max_batches, args.pore, args.qscore)
    plt_parameters(paramFile, outdir)

if __name__ == '__main__':
    main()