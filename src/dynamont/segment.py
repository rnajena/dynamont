#!/usr/bin/env python
# author: Jannes Spangenberg
# e-mail: jannes.spangenberg@uni-jena.de
# github: https://github.com/JannesSP
# website: https://jannessp.github.io

from argparse import ArgumentDefaultsHelpFormatter, ArgumentParser, Namespace
from os.path import exists, join, dirname
from os import makedirs, name
from FileIO import getFiles, loadFastx, genSegmentTable
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

# def segment(rawdatapath : str, fastxpath : str, polya : dict, batch_size : int, mode : str, outfile : str, modelpath : str) -> None:
def segment(rawdatapath : str, fastxpath : str, batch_size : int, mode : str, outfile : str, modelpath : str) -> None:
    files = getFiles(rawdatapath, True)
    print(f'ONT Files: {len(files)}')
    basecalls = loadFastx(fastxpath)
    print(f'{len(basecalls)} reads found')

    out = open(outfile, 'w')
    
    if mode == 'indel':
        CPP_SCRIPT = join(dirname(__file__), 'segmentation_indel')
        if name == 'nt': # check for windows
            CPP_SCRIPT+='.exe'
    elif mode == 'basic':
        CPP_SCRIPT = join(dirname(__file__), 'segmentation_basic')
        if name == 'nt': # check for windows
            CPP_SCRIPT+='.exe'
    i = 0

    out.write('readid,start,end,basepos,base,motif,state\n')

    with mp.Pool(batch_size) as p:
        for file in files:

            r5 = read(file)
            mp_items = []

            for readid in r5:
                if not readid in basecalls: # or (readid not in polyAIndex)
                    continue

                # skip reads with undetermined transcript start, only take really good reads for training
                # if not readid in polya or polya[readid][1] == -1 or polya[readid][1] - polya[readid][0] < 30:
                #     signal = r5.getPolyAStandardizedSignal(readid, polya[readid][0], polya[readid][1])[polya[readid][1]:]
                # else:

                signal = r5.getpASignal(readid)
                
                # if polya[readid][1] - polya[readid][0] < 30:
                #     continue

                # signal = r5.getPolyAStandardizedSignal(readid, polya[readid][0], polya[readid][1])[polya[readid][1]:]
                mp_items.append([signal, basecalls[readid][::-1], CPP_SCRIPT, readid, modelpath])

                if len(mp_items)%batch_size == 0:
                    for result in p.starmap(genSegmentTable, mp_items):
                        out.write(result)
                        i += 1
                    print(f"Segmented {i} reads", end='\r')
                    # initialize new batch
                    mp_items = []

            # last unfilled batch
            if len(mp_items):
                for result in p.starmap(genSegmentTable, mp_items):
                    out.write(result)
                    i += 1
                print(f"Segmented {i} reads")
                # initialize new batch
                mp_items = []

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