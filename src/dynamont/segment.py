#!/usr/bin/env python
# author: Jannes Spangenberg
# e-mail: jannes.spangenberg@uni-jena.de
# github: https://github.com/JannesSP
# website: https://jannessp.github.io

from argparse import ArgumentDefaultsHelpFormatter, ArgumentParser, Namespace
from os.path import exists, join, dirname
from os import makedirs, name
from FileIO import getFiles, loadFastx, readPolyAEnd, genSegmentTable
from read5 import read
import multiprocessing as mp
from hampel import hampel

def parse() -> Namespace:
    parser = ArgumentParser(
        formatter_class=ArgumentDefaultsHelpFormatter
    )
    parser.add_argument('--raw', type=str, required=True, help='Raw ONT training data')
    parser.add_argument('--fastx', type=str, required=True, help='Basecalls of ONT training data')
    parser.add_argument('--polya', type=str, required=True, help='Poly A table from nanopolish polya containing the transcript starts')
    parser.add_argument('--out', type=str, required=True, help='Outpath to write files')
    parser.add_argument('--batch_size', type=int, default=16, help='Number of reads to train before updating')
    parser.add_argument('--mode', choices=['basic', 'indel'], required=True)
    return parser.parse_args()

def train(rawdatapath : str, fastxpath : str, polya : dict, batch_size : int, mode : str, outfile : str) -> None:
    files = getFiles(rawdatapath, True)
    print(f'ONT Files: {len(files)}')
    basecalls = loadFastx(fastxpath)
    print(f'Training segmentation parameters with {len(basecalls)} reads')

    out = open(outfile, 'w')
    
    if mode == 'indel':
        CPP_SCRIPT = join(dirname(__file__), 'segmentation_indel')
        if name == 'nt': # check for windows
            CPP_SCRIPT+='.exe'
        # init
        # params = {
        #     "e1":1.,
        #     "m2":.03333,
        #     "d1":.00001,
        #     "e2":.96664,
        #     "e3":.00001,
        #     "i1":.00001,
        #     "m3":.99,
        #     "i2":.01,
        #     "m4":.99,
        #     "d2":.01,
        #     # "s1":0.16,
        #     # "s2":1.0
        # }
    elif mode == 'basic':
        CPP_SCRIPT = join(dirname(__file__), 'segmentation_basic')
        if name == 'nt': # check for windows
            CPP_SCRIPT+='.exe'
        # params = {
        #     "e1":1.,
        #     "m2":.33,
        #     "e2":.33,
        #     "e3":.33,
        # }

    i = 0
    out.write('readid,start,end,basepos,base,state\n')

    with mp.Pool(batch_size) as p:
        for file in files:

            r5 = read(file)
            mp_items = []

            for readid in r5:
                if not readid in basecalls: # or (readid not in polyAIndex)
                    continue

                # currently depending on nanopolish polya TODO improvement room here
                if polya[readid] == -1:
                    polya[readid] = 0

                signal = hampel(r5.getpASignal(readid)[polya[readid]:], 20, 2.).filtered_data
                mp_items.append([signal, basecalls[readid][::-1], CPP_SCRIPT, readid])

                if len(mp_items)%batch_size == 0:
                    for result in p.starmap(genSegmentTable, mp_items):
                        out.write(result)
                        i += 1
                    print(f"Segmented {i} reads")
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
    polya=readPolyAEnd(args.polya)
    train(args.raw, args.fastx, polya, args.batch_size, args.mode)

if __name__ == '__main__':
    main()