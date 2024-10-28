#!/usr/bin/env python
# author: Jannes Spangenberg
# e-mail: jannes.spangenberg@uni-jena.de
# github: https://github.com/JannesSP
# website: https://jannessp.github.io

from argparse import ArgumentDefaultsHelpFormatter, ArgumentParser, Namespace
import os
import h5py
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from collections import defaultdict

def parse() -> Namespace:
    parser = ArgumentParser(
        formatter_class=ArgumentDefaultsHelpFormatter
    )
    parser.add_argument("--changepoints", type=str, required=True, metavar="HDF5", help="HDF5 file with ground truth change points")
    parser.add_argument("--dynamont", type=str, required=True, metavar="CSV", help="Dynamont segmentation output in CSV format")
    parser.add_argument("--dorado", type=str, required=True, metavar="TSV", help="Dorado segmentation output in TSV format")
    parser.add_argument("--f5cresquiggle", type=str, required=True, metavar="TSV", help="f5c resquiggle segmentation output in TSV format")
    parser.add_argument("--f5ceventalign", type=str, required=True, metavar="TSV", help="f5c eventalign segmentation output in TSV format\nSummary file must exists in the same path with extension .sum")
    parser.add_argument("--tombo", type=str, required=True, metavar="PATH", help="Parent directory of single fast5s processed by Tombo")
    parser.add_argument("--output", type=str, required=True, metavar="CSV", help="Output CSV containing pandas data frame")
    return parser.parse_args()

def readDynamont(file: str) -> dict:
    readMap = defaultdict(set)
    with open(file, 'r') as f:
        next(f) # skip header
        # readid,signalid,start,end,basepos,base,motif,state,posterior_probability,polish
        for line in f:
            _, signalid, start, end, *_ = line.strip().split(',')
            readMap[signalid].update([int(start), int(end)])
    return readMap

def readDorado(file: str) -> dict:
    readMap = defaultdict(set)
    with open(file, 'r') as f:
        next(f) # skip header
        # readid    position    base    motif   start   end
        for line in f:
            readid, _, _, _, start, end = line.strip().split('\t')
            readMap[readid].update([int(start), int(end)])
    return readMap

def readF5CEventalign(file: str, summary : str) -> dict:
    readIdMap = {}
    with open(summary, 'r') as f:
        next(f) # skip header
        # read_index    read_name   fast5_path  model_name  strand  num_events  num_steps   num_skips   num_stays   total_duration  shift   scale   drift   var
        for line in f:
            readid, read_name, *rest = line.strip().split('\t')
            readIdMap[readid] = read_name

    readMap = defaultdict(set)
    with open(file, 'r') as f:
        next(f) # skip header
        # contig    position    reference_kmer  read_index  strand  event_index event_level_mean    event_stdv  event_length    model_kmer  model_mean  model_stdv  standardized_level  start_idx   end_idx
        for line in f:
            _, _, _, readid, _, _, _, _, _, _, _, _, _, start, end = line.strip().split('\t')
            readMap[readIdMap[readid]].update([int(start), int(end)])
    return readMap

def readF5CResquiggle(file: str) -> dict:
    readMap = defaultdict(set)
    with open(file, 'r') as f:
        next(f) # skip header
        # contig    position    reference_kmer  read_index  strand  event_index event_level_mean    event_stdv  event_length    model_kmer  model_mean  model_stdv  standardized_level  start_idx   end_idx
        for line in f:
            readid, _, start, end = line.strip().split('\t')
            readMap[readid].update([int(start), int(end)])
    return readMap

def getFast5s(directory: str) -> list:
    fast5_files = []
    for root, _, files in os.walk(directory):
        for file in files:
            if file.endswith('.fast5'):
                fast5_files.append(os.path.join(root, file))
    return fast5_files

def readTombo(directory: str) -> dict:
    fast5s = getFast5s(directory)
    readMap = defaultdict(set)
    for fast5 in fast5s:
        with h5py.File(fast5, 'r') as h5:
            readid = os.path.basename(fast5).split('.')[0]
            events : np.ndarray = h5['Analyses/RawGenomeCorrected_000/BaseCalled_template/Events'][:]
            readMap[readid].update(events['start'].tolist())
    return readMap

def readChangepoints(file : str) -> dict:
    readMap = defaultdict(set)
    with h5py.File(file, 'r') as h5:
        for readid in h5.keys():
            changepoints = h5[readid][:]
            readMap[readid].update(changepoints)
    return readMap

def toNumpy(readMap: dict) -> None:
    for readid in readMap:
        readMap[readid] = np.array(list(readMap[readid]))
        readMap[readid].sort()

def evaluate(groundTruth: dict, prediction: dict, maxDistance : int) -> dict:
    found = 0
    for gt_value in groundTruth:
        min_distance = min(np.abs(gt_value - prediction))
        found += min_distance < maxDistance
    return found

def plot(df: pd.DataFrame, outfile : str) -> None:
    sns.set_theme()
    sns.barplot(df, x='maxDistance', y='foundEdges', hue='tool', log_scale=True)
    plt.savefig(outfile)

def main() -> None:
    args = parse()

    if not os.path.exists(args.output):
        groundTruths = readChangepoints(args.changepoints)

        toolsResult = {
            'dynamont' : readDynamont(args.dynamont),
            'f5cEventalign' : readF5CEventalign(args.f5ceventalign, os.path.splitext(args.f5ceventalign)[0] + '.sum'),
            'f5cResquiggle' : readF5CResquiggle(args.f5cresquiggle),
            'dorado' : readDorado(args.dorado),
            'tombo' : readTombo(args.tombo),
        }

        for tool in toolsResult:
            toNumpy(toolsResult[tool])

        df = pd.DataFrame()
        for tool, tool_results in toolsResult.items():
            print(f"Evaluating {tool}...")

            for maxDistance in range(0, 1000, 1):
                print(f'Distance: {maxDistance}')
                totalEdges = 0 # total edges in ground truth
                foundEdges = 0 # found edges by tool
                totalReads = 0 # total reads in ground truth
                segmentedReads = 0 # segmented reads by tool
            
                for i, (readid, gt_set) in enumerate(groundTruths.items()):
                    if (i+1) % 1000 == 0:
                        print(f'{i+1}/{len(groundTruths)}', end='\r')
                    totalEdges += len(gt_set)
                    totalReads += 1

                    if readid in tool_results and tool_results[readid]:
                        segmentedReads += 1
                        foundEdges += evaluate(gt_set, tool_results, maxDistance)

                new_entry = pd.DataFrame({
                    'tool': [tool],
                    'maxDistance': [maxDistance],
                    'foundEdges': [foundEdges],
                    'totalEdges': [totalEdges],
                    'foundEdgesRatio': [foundEdges / totalEdges if totalEdges > 0 else 0],
                    'segmentedReads': [segmentedReads],
                    'totalReads': [totalReads],
                    'segmentedReadsRatio': [segmentedReads / totalReads if totalReads > 0 else 0],
                })
                df = pd.concat((df, new_entry), ignore_index=True)

            print(f'Done: {i}/{len(groundTruths)}')

        df.to_csv(args.output, index=False)    

    else:
        df = pd.read_csv(args.output, index=False)

    print('Plotting...')
    plot(df, os.path.splitext(args.output)[0] + '.svg')

if __name__ == '__main__':
    main()