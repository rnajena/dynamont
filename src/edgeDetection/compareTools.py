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
import multiprocessing as mp

# TODO: add support for DNA
# ! - add rna/dna parameter
# ! - add different data processing in readTombo

def parse() -> Namespace:
    """
    Parse command line arguments for tool comparison.

    Returns:
        Namespace: Containing the specified command line arguments
    """
    parser = ArgumentParser(
        formatter_class=ArgumentDefaultsHelpFormatter
    )
    parser.add_argument("--changepoints", type=str, required=True, metavar="HDF5", help="HDF5 file with ground truth change points")
    parser.add_argument("--dorado", type=str, required=True, metavar="TSV", help="Dorado segmentation output in TSV format")
    parser.add_argument("--f5cresquiggle", type=str, required=True, metavar="TSV", help="f5c resquiggle segmentation output in TSV format")
    parser.add_argument("--f5ceventalign", type=str, required=True, metavar="TSV", help="f5c eventalign segmentation output in TSV format\nSummary file must exists in the same path with extension .sum")
    parser.add_argument("--tombo", type=str, required=True, metavar="PATH", help="Parent directory of single fast5s processed by Tombo")
    parser.add_argument("--basecalls", type=str, required=True, metavar="BAM", help="Basecalls of ONT training data as .bam file")
    # parser.add_argument("--pod5", type=str, required=True, metavar="POD5", help="raw signals")
    parser.add_argument("--output", type=str, required=True, metavar="CSV", help="Output CSV containing pandas data frame")
    parser.add_argument("--dynamont", type=str, nargs='+', required=True, metavar="CSV", help="Dynamont segmentation output in CSV format")
    parser.add_argument("--labels", type=str, nargs='+', required=True, metavar="CSV", help="Dynamont segmentation output in CSV format")
    parser.add_argument("-p", "--processes", type=int, default=7, metavar="INT", help="Number of processes to use for parallel processing (default: all available CPUs)")
    return parser.parse_args()

def readDynamont(file: str) -> dict:
    """
    Reads a CSV file containing Dynamont segmentation output and extracts
    signal IDs along with their corresponding start and end positions.

    Parameters
    ----------
    file : str
        Path to the CSV file containing Dynamont segmentation output.

    Returns
    -------
    dict
        A dictionary mapping each signal ID to a set of start and end positions.
    """
    print("Reading Dynamont output from " + file)
    readMap = defaultdict(set)
    with open(file, 'r') as f:
        next(f) # skip header
        # readid,signalid,start,end,basepos,base,motif,state,posterior_probability,polish
        for line in f:
            try:
                _, signalid, start, end, *_ = line.strip().split(',')
            except ValueError: # empty line in file
                pass
            readMap[signalid].update([int(start), int(end)])
    return readMap

def readDorado(file: str) -> dict:
    """
    Reads a TSV file containing Dorado segmentation output and extracts
    signal IDs along with their corresponding start and end positions.

    Parameters
    ----------
    file : str
        Path to the TSV file containing Dorado segmentation output.

    Returns
    -------
    dict
        A dictionary mapping each signal ID to a set of start and end positions.
    """
    print("Reading Dorado output from " + file)
    readMap = defaultdict(set)
    with open(file, 'r') as f:
        next(f) # skip header
        # readid    position    base    motif   start   end
        for line in f:
            readid, _, _, _, start, end = line.strip().split('\t')
            readMap[readid].update([int(start), int(end)])
    return readMap

def readF5CEventalign(file: str, summary : str) -> dict:
    """
    Reads a TSV file containing f5c eventalign segmentation output and extracts
    read IDs along with their corresponding start and end positions.

    Parameters
    ----------
    file : str
        Path to the TSV file containing f5c eventalign segmentation output.
    summary : str
        Path to the .sum file containing the read names corresponding to the read indices in the segmentation output.

    Returns
    -------
    dict
        A dictionary mapping each read ID to a set of start and end positions.
    """
    print("Reading f5c eventalign output from " + file)
    readIdMap = {}
    with open(summary, 'r') as f:
        next(f) # skip header
        # read_index    read_name   fast5_path  model_name  strand  num_events  num_steps   num_skips   num_stays   total_duration  shift   scale   drift   var
        for line in f:
            readid, readName, *rest = line.strip().split('\t')
            readIdMap[readid] = readName

    readMap = defaultdict(set)
    with open(file, 'r') as f:
        next(f) # skip header
        # contig    position    reference_kmer  read_index  strand  event_index event_level_mean    event_stdv  event_length    model_kmer  model_mean  model_stdv  standardized_level  start_idx   end_idx
        for line in f:
            _, _, _, readid, _, _, _, _, _, _, _, _, _, start, end = line.strip().split('\t')
            readMap[readIdMap[readid]].update([int(start), int(end)])
    return readMap

def readF5CResquiggle(file: str) -> dict:
    """
    Parses a TSV file to extract read IDs along with their corresponding start and end positions.

    Parameters
    ----------
    file : str
        Path to the TSV file containing resquiggle segmentation output.

    Returns
    -------
    dict
        A dictionary mapping each read ID to a set of start and end positions.
    """
    print("Reading f5c resquiggle output from " + file)
    readMap = defaultdict(set)
    with open(file, 'r') as f:
        next(f) # skip header
        # contig    position    reference_kmer  read_index  strand  event_index event_level_mean    event_stdv  event_length    model_kmer  model_mean  model_stdv  standardized_level  start_idx   end_idx
        for line in f:
            # print(line)
            readid, _, start, end = line.strip().split('\t')
            try:
                readMap[readid].update([int(start), int(end)])
            except ValueError: # no segmentation for that position, e.g. 00b77311-c22c-472e-ad39-59f5f6aea4ec    467     .       .
                pass
    return readMap

def getFast5s(directory: str) -> list:
    """
    Recursively walks through a directory and its subdirectories to find all FAST5 files.

    Parameters
    ----------
    directory : str
        The root directory to start searching from.

    Returns
    -------
    list
        A list of paths to all FAST5 files found in the directory and its subdirectories.
    """
    print("Reading FAST5 files from " + directory)
    fast5Files = []
    for root, _, files in os.walk(directory):
        # check if root is absolute path
        if not root.startswith("/"):
            root = os.path.join(os.getcwd(), root)
        for file in files:
            if file.endswith('.fast5'):
                filepath = os.path.join(root, file)
                if not os.path.exists(filepath):
                    print("File not found, path may be broken!: " + filepath)
                    exit(1)
                fast5Files.append(filepath)
    return fast5Files

def readTombo(directory: str) -> dict:
    """
    Recursively walks through a directory and its subdirectories to find all FAST5 files, then parses the 'RawGenomeCorrected_000/BaseCalled_template/Events' group in each FAST5 file to extract read IDs along with their corresponding start positions.

    Parameters
    ----------
    directory : str
        The root directory to start searching from.

    Returns
    -------
    dict
        A dictionary mapping each read ID to a set of start positions.
    """
    print("Reading tombo output from " + directory)
    fast5s = getFast5s(directory)
    readMap = defaultdict(set)
    for fast5 in fast5s:
        with h5py.File(fast5, 'r') as h5:
            readid = os.path.basename(fast5).split('.')[0]
            try:
                if h5['Analyses/RawGenomeCorrected_000/BaseCalled_template'].attrs['status'] == 'Alignment not produced':
                    continue
                # https://nanoporetech.github.io/tombo/resquiggle.html
                # Minor RNA note: RNA reads pass through the pore in the 3’ to 5’ direction during sequencing. As such, the raw signal and albacore events are stored in the reverse direction from DNA reads. Tombo events for RNA data are stored in the opposite direction (corresponding to the genome sequence direction, not sequencing time direction) for several practical reasons. Thus if events are to be compared to the raw signal, the raw signal must be reversed. Tombo RNA models are stored in the same direction and thus may be considered inverted as compared to some other RNA HMM signal level models.
                signalLen = len(h5['Raw/Reads'][list(h5['Raw/Reads'].keys())[0]]['Signal'][:])
                starts = h5['Analyses/RawGenomeCorrected_000/BaseCalled_template/Events'][:]['start'] + h5['Analyses/RawGenomeCorrected_000/BaseCalled_template/Events'].attrs['read_start_rel_to_raw']
                ends = starts + h5['Analyses/RawGenomeCorrected_000/BaseCalled_template/Events'][:]['length']
                borders = np.unique((starts, ends))
                # reverse indices
                borders = signalLen - borders - 1
            except KeyError:
                continue
            readMap[readid].update(borders.tolist())
    return readMap

def readChangepoints(file : str) -> dict:
    """
    Reads a HDF5 file containing changepoints and extracts read IDs along with their corresponding changepoints.

    Parameters
    ----------
    file : str
        Path to the HDF5 file containing changepoints.

    Returns
    -------
    dict
        A dictionary mapping each read ID to a set of changepoints.
    """
    print("Reading changepoints from " + file)
    readMap = defaultdict(set)
    with h5py.File(file, 'r') as h5:
        for readid in h5.keys():
            # print(h5[readid].keys())
            changepoints = h5[readid + "/waveletEdge"][:]
            readMap[readid].update(changepoints)
    return readMap

def toNumpy(readMap: dict) -> None:
    """
    Converts all values in the given readMap from a set to a sorted numpy array.

    Parameters
    ----------
    readMap : dict
        A dictionary mapping read IDs to a set of positions.
    """
    for readid in readMap:
        readMap[readid] = np.array(list(readMap[readid]))
        readMap[readid].sort()

def evaluate(groundTruth: np.ndarray, prediction: np.ndarray, maxDistance: int) -> int:
    """
    Masterclass solution to find the number of ground truth values
    with a closest prediction within a specified maxDistance using two pointers.
    Runtime: O(n + m)
    Space: O(1)

    Parameters
    ----------
    groundTruth : np.ndarray
        Sorted array of ground truth change point positions.
    prediction : np.ndarray
        Sorted array of predicted change point positions.
    maxDistance : int
        Maximum allowable distance for a match.

    Returns
    -------
    int
        Count of ground truth values with a nearby prediction within maxDistance.
    """
    gt_idx = pred_idx = found = before = after = 0
    
    # Two-pointer technique
    while gt_idx < len(groundTruth) and pred_idx < len(prediction):
        gt_val = groundTruth[gt_idx]
        pred_val = prediction[pred_idx]
        
        # Check the distance
        if abs(gt_val - pred_val) <= maxDistance:
            found += 1
            # Check if the prediction is before or after the ground truth
            if pred_val < gt_val:
                before += 1
            else:
                after += 1
            gt_idx += 1  # Move to next ground truth point as it's already matched
        elif pred_val < gt_val:
            pred_idx += 1  # Move prediction pointer to the right
        else:
            gt_idx += 1  # Move ground truth pointer to the right

    return found, before, after

def plot(df: pd.DataFrame, outfile : str) -> None:
    """
    Plots a bar plot of the found edges for each tool at each maxDistance.

    Parameters
    ----------
    df : pd.DataFrame
        A DataFrame containing the found edges for each tool at each maxDistance.
    outfile : str
        The path to save the plot to.
    """
    sns.set_theme()

    plt.rcParams.update({'axes.labelsize': 15, 'xtick.labelsize': 15, 'ytick.labelsize': 15, 'legend.fontsize': 15, 'legend.title_fontsize': 15})

    plt.figure(figsize=(12,12), dpi=300)
    sns.lineplot(df, x='maxDistance', y='foundEdges', hue='tool')
    plt.xlabel("Distance Treshold")
    plt.ylabel("Found Edges Ratio")
    plt.savefig(outfile + "foundEdges.svg", dpi=300)
    plt.xlim((0, 10))
    plt.savefig(outfile + "foundEdges_0-10.svg", dpi=300)
    plt.close()
    
    plt.figure(figsize=(12,12), dpi=300)
    sns.lineplot(df, x='maxDistance', y='foundEdgesRatio', hue='tool')
    plt.xlabel("Distance Treshold")
    plt.ylabel("Found Edges Ratio")
    plt.savefig(outfile + "foundEdgesRatio.svg", dpi=300)
    plt.xlim((0, 10))
    plt.savefig(outfile + "foundEdgesRatio_0-10.svg", dpi=300)
    plt.close()

    plt.figure(figsize=(12,12), dpi=300)
    sns.barplot(df.loc[df["maxDistance"] == 0], x='tool', y='segmentedReadsRatio')
    # Annotate each bar with the height
    for index, (_, row) in enumerate(df.loc[df["maxDistance"] == 0].iterrows()):
        plt.text(x=index, y=row['segmentedReadsRatio'],  # Add a small offset for better visibility
             s=f"{row['segmentedReadsRatio']:.3f}",        # Format the text with 2 decimals
             ha='center',                     # Center the text horizontally
             va='bottom')                     # Align the text to the top of the bar
    plt.ylim((0.8, 1.01))
    plt.yticks(np.arange(0.8, 1.01, 0.05))
    plt.xlabel("Distance Treshold")
    plt.ylabel("Found Edges Ratio")
    plt.savefig(outfile + "segmentedReadsRatio.svg", dpi=300)
    plt.close()

def generateControl(bamFile : str) -> tuple:
    """
    Generates a control set for a given bam file.

    For each read in the bam file, two types of control sets are generated:
    - A set of randomly drawn borders without replacement.
    - A set of equidistant borders.

    The start and end positions of each read are stored in a dictionary for later use.

    Returns
    -------
    tuple
        A tuple containing two dictionaries: randomBorders and uniformBorders.
        Each dictionary contains a list of borders for each read in the bam file.
    """
    print("Generating control set from " + bamFile)
    import pysam
    basecalledRegions = {}
    with pysam.AlignmentFile(bamFile, "rb", check_sq=False) as bamfile:
        for read in bamfile.fetch(until_eof=True):
            readid = read.get_tag("pi") if read.has_tag("pi") else read.query_name
            sp = read.get_tag("sp") if read.has_tag("sp") else 0 # if split read get start offset of the signal
            ts = read.get_tag("ts") # ts:i: 	the number of samples trimmed from the start of the signal
            ns = read.get_tag("ns") # ns:i: 	the number of samples in the signal prior to trimming
            nts = len(read.query_sequence)
            
            if readid not in basecalledRegions:
                basecalledRegions[readid] = []

            # store information: number of nucleotides, start position in signal, end position in signal
            basecalledRegions[readid].append([nts, sp+ts, sp+ns])

    randomBorders = {} # randomly drawn borders without replacement
    uniformBorders = {} # equidistant borders
    for readid in basecalledRegions:
        regions = np.array([])
        uniformBorders[readid] = []
        nts = 0
        for region in basecalledRegions[readid]:
            nts += region[0]
            regions = np.append(regions, np.arange(region[1], region[2], 1))
            uniformBorders[readid].extend(np.linspace(region[1], region[2] - 1, region[0]).tolist())

        randomBorders[readid] = np.random.choice(regions, size=nts, replace=False)
        uniformBorders[readid] = np.round(np.array(uniformBorders[readid])).astype(int)
        
    return randomBorders, uniformBorders

def main() -> None:
    args = parse()
    print(args)

    pool = mp.Pool(args.processes)

    assert len(args.dynamont) == len(args.labels), "Number of dynamont results must match the number of labels"

    if not os.path.exists(args.output):
        gtReturn = pool.apply_async(readChangepoints, args=(args.changepoints,))
        # dynamontReturn = pool.apply_async(readDynamont, args=(args.dynamont,))
        dynamontReturn = {label : pool.apply_async(readDynamont, args=(result,)) for label, result in zip(args.labels, args.dynamont)}
        f5ceventalignReturn = pool.apply_async(readF5CEventalign, args=(args.f5ceventalign, os.path.splitext(args.f5ceventalign)[0] + '.sum'))
        f5cresquiggleReturn = pool.apply_async(readF5CResquiggle, args=(args.f5cresquiggle,))
        doradoReturn = pool.apply_async(readDorado, args=(args.dorado,))
        tomboReturn = pool.apply_async(readTombo, args=(args.tombo,))
        controlReturn = pool.apply_async(generateControl, args=(args.basecalls,))

        # pool.close()
        # pool.join()

        print("Done Reading,\tStart Evaluating...")

        groundTruths = gtReturn.get()  #readChangepoints(args.changepoints)
        # generating controls
        randomBorders, uniformBorders = controlReturn.get()

        toolsResult = {
            f"dynamont_{label}" : dynamontReturn[label].get() for label in dynamontReturn
        } | {
            'f5cEventalign' : f5ceventalignReturn.get(), # readF5CEventalign(args.f5ceventalign, os.path.splitext(args.f5ceventalign)[0] + '.sum'),
            'f5cResquiggle' : f5cresquiggleReturn.get(), # readF5CResquiggle(args.f5cresquiggle),
            'dorado' : doradoReturn.get(), # readDorado(args.dorado),
            'tombo' : tomboReturn.get(), # readTombo(args.tombo),
            'random' : randomBorders,
            'uniform' : uniformBorders,
        }

        for tool in toolsResult:
            toNumpy(toolsResult[tool])
        toNumpy(groundTruths)

        # df = pd.DataFrame()
        outfile = open(args.output, 'w')
        outfile.write("tool,maxDistance,foundEdges,foundBefore,foundAfter,totalEdges,foundEdgesRatio,segmentedReads,totalReads,segmentedReadsRatio\n")
        outfile.flush()
        for tool, result in toolsResult.items():
            print(f"Evaluating {tool}...")

            for maxDistance in range(0, 50, 1):
                print(f'Distance: {maxDistance}')
                totalEdges = 0 # total edges in ground truth
                foundEdges = 0 # found edges by tool
                foundBefore = 0
                foundAfter = 0
                totalReads = 0 # total reads in ground truth
                segmentedReads = 0 # segmented reads by tool
            
                jobs = []
                for i, (readid, gtSet) in enumerate(groundTruths.items()):
                    totalEdges += len(gtSet)
                    totalReads += 1

                    if readid in result and result[readid].size:
                        segmentedReads += 1
                        jobs.append(pool.apply_async(evaluate, args=(gtSet, result[readid], maxDistance)))

                for i, job in enumerate(jobs):
                    if (i+1) % 1000 == 0:
                        print(f'{i+1}/{len(groundTruths)}', end='\r')
                    result = job.get()
                    foundEdges += result[0]
                    foundBefore += result[1]
                    foundAfter += result[2]
                    
                outfile.write(f"{tool},{maxDistance},{foundEdges},{foundBefore},{foundAfter},{totalEdges},{foundEdges / totalEdges if totalEdges > 0 else 0},{segmentedReads},{totalReads},{segmentedReads / totalReads if totalReads > 0 else 0}\n")
                outfile.flush()

            print(f'Done: {i+1}/{len(groundTruths)}')
        
        outfile.close()
    
    pool.close()
    pool.join()

    print("Reading file " + args.output)
    df = pd.read_csv(args.output)
    print('Plotting...')
    plot(df, os.path.splitext(args.output)[0] + "_")

if __name__ == '__main__':
    main()