#!/usr/bin/env python
# author: Jannes Spangenberg
# e-mail: jannes.spangenberg@uni-jena.de
# github: https://github.com/JannesSP
# website: https://jannessp.github.io

import os
import h5py
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import multiprocessing as mp
from argparse import ArgumentDefaultsHelpFormatter, ArgumentParser, Namespace
from collections import defaultdict
import itertools
# from venn import venn
import read5_ont
from scipy.stats import median_abs_deviation as mad

# TODO: add support for DNA

def parse() -> Namespace:
    """
    Parse command line arguments for tool comparison.

    Returns:
        Namespace: Containing the specified command line arguments
    """
    parser = ArgumentParser(
        formatter_class=ArgumentDefaultsHelpFormatter
    )
    # parser.add_argument("--changepoints", type=str, required=True, metavar="HDF5", help="HDF5 file with ground truth change points")
    parser.add_argument("--dorado", type=str, required=True, metavar="TSV", help="Dorado segmentation output in TSV format")
    parser.add_argument("--f5cresquiggle", type=str, required=True, metavar="TSV", help="f5c resquiggle segmentation output in TSV format")
    parser.add_argument("--f5ceventalign", type=str, required=True, metavar="TSV", help="f5c eventalign segmentation output in TSV format\nSummary file must exists in the same path with extension .sum")
    parser.add_argument("--uncalled", type=str, metavar="TSV", required=True, help="Uncalled4 segmentation output in csv format")
    parser.add_argument("--basecalls", type=str, required=True, metavar="BAM", help="Basecalls of ONT training data as .bam file")
    parser.add_argument("--pod5", type=str, required=True, metavar="POD5", help="raw signals")
    parser.add_argument("--output", type=str, required=True, metavar="CSV", help="Output CSV containing pandas data frame")
    parser.add_argument("--dynamont", type=str, nargs='+', required=True, metavar="CSV", help="Dynamont segmentation output in CSV format")
    parser.add_argument("--labels", type=str, nargs='+', required=True, metavar="CSV", help="Dynamont segmentation output in CSV format")
    parser.add_argument("--tombo", type=str, metavar="PATH", help="Parent directory of single fast5s processed by Tombo")
    # parser.add_argument("--distance", type=int, default=50, metavar="INT", help="Maximum distance to search for")
    parser.add_argument("-p", "--processes", type=int, default=7, metavar="INT", help="Number of processes to use for parallel processing (default: all available CPUs)")
    parser.add_argument("-w", "--window", type=int, default=0, metavar="INT", help="Window size to use for comparison")
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
            readid, signalid, _, _, _, start, end = line.strip().split('\t')
            readMap[signalid].update([int(start), int(end)])
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

def readUncalled4(file : str) -> dict:
    print("Reading uncalled 4 output from " + file)
    readMap = defaultdict(set)
    with open(file, 'r') as f:
        next(f) # skip header
        # seq.name        seq.pos seq.strand      aln.id  seq.kmer        aln.read_id     dtw.start       dtw.length      dtw.current     dtw.current_sd
        for line in f:
            _, _, _, _, _, readid, start, length, *_ = line.strip().split('\t')
            if start == "*":
                continue
            readMap[readid].update([int(start), int(start) + int(length)])
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
        readMap[readid] = np.array(list(readMap[readid]), dtype=int)
        readMap[readid].sort()

# def evaluate(groundTruth: np.ndarray, prediction: np.ndarray, maxDistance: int) -> int:
#     """
#     Masterclass solution to find the number of ground truth values
#     with a closest prediction within a specified maxDistance using two pointers.
#     Runtime: O(n + m)
#     Space: O(1)

#     Parameters
#     ----------
#     groundTruth : np.ndarray
#         Sorted array of ground truth change point positions.
#     prediction : np.ndarray
#         Sorted array of predicted change point positions.
#     maxDistance : int
#         Maximum allowable distance for a match.

#     Returns
#     -------
#     int
#         Count of ground truth values with a nearby prediction within maxDistance.
#     """
#     gt_idx = pred_idx = found = 0
    
#     # Two-pointer technique
#     while gt_idx < len(groundTruth) and pred_idx < len(prediction):
#         gt_val = groundTruth[gt_idx]
#         pred_val = prediction[pred_idx]
        
#         # Check the distance
#         if abs(gt_val - pred_val) <= maxDistance:
#             found += 1
#             gt_idx += 1  # Move to next ground truth point as it's already matched
#         elif pred_val < gt_val:
#             pred_idx += 1  # Move prediction pointer to the right
#         else:
#             gt_idx += 1  # Move ground truth pointer to the right

#     return found

def evaluate(gt : np.ndarray, pred : np.ndarray, maxDistance : int) -> pd.DataFrame:
    """
    gt : sorted np.ndarray
    pred : sorted np.ndarray
    maxDistance : int
    """
    
    found = np.array([{} for _ in range(len(gt))])

    # collect all predictions within maxDistance for each ground truth
    pred_start = 0
    for gt_idx in range(len(gt)):
        gt_val = int(gt[gt_idx])

        for pred_idx in range(pred_start, len(pred)):
            pred_val = int(pred[pred_idx])

            if gt_val - pred_val > maxDistance:
                pred_idx+=1
                continue
            elif gt_val - pred_val < -maxDistance:
                break

            found[gt_idx][pred_idx] = gt_val - pred_val

        # start next loop at min position of previous gt value
        pred_start = min(found[gt_idx].keys()) if found[gt_idx] else 0

    result = np.zeros(2 * maxDistance + 1, dtype=int)
    for gt_idx in range(len(gt)):

        # no prediction for this gt value
        if not found[gt_idx]:
            continue

        cur_gt = found[gt_idx]

        # get closest prediction
        pred_idx = min(cur_gt, key=lambda k: abs(cur_gt[k]))

        prev_gt = found[gt_idx - 1] if gt_idx > 0 else {}
        # remove current prediction from current gt if prediction is closer to previous ground truth value
        # should only happen once
        if prev_gt and cur_gt and pred_idx in prev_gt and abs(prev_gt[pred_idx]) <= abs(cur_gt[pred_idx]):
            del cur_gt[pred_idx]

            if cur_gt:
                pred_idx = min(cur_gt, key=lambda k: abs(cur_gt[k]))

        next_gt = found[gt_idx + 1] if gt_idx < len(gt) - 1 else {}
        if next_gt and cur_gt and pred_idx in next_gt and abs(next_gt[pred_idx]) < abs(cur_gt[pred_idx]):
        # remove current prediction from current gt if prediction is closer to succeeding ground truth value
        # should only happen once
            del cur_gt[pred_idx]

            if cur_gt:
                pred_idx = min(cur_gt, key=lambda k: abs(cur_gt[k]))

        # no prediction for this gt value
        if not cur_gt:
            continue

        minDist = int(cur_gt[pred_idx])

        if minDist >= 0:
            # Map the difference to the appropriate index
            result[: maxDistance - minDist + 1] += 1
        if minDist <= 0:
            # Map the difference to the appropriate index
            result[maxDistance - minDist :] += 1
        if minDist == 0:
            result[maxDistance] -= 1

    return result

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

    # outfile.write("Tool,Distance Threshold,Found Change Points,Total Change Points,Found Change Points Ratio,Segmented Reads,Total Reads,Segmented Reads Ratio\n")

    # Get the unique order of Tools based on their appearance in the original DataFrame
    tool_order = df['Tool'].unique()

    plt.figure(figsize=(8,8), dpi=300)
    sns.lineplot(df, x='Distance Threshold', y='Found Change Points', hue='Tool', style='Tool')
    plt.title("Found Change Points vs Directional Distance Threshold", fontsize=18)
    plt.savefig(outfile + "foundEdgesDirectional.svg", dpi=300)
    plt.savefig(outfile + "foundEdgesDirectional.pdf", dpi=300)
    plt.xlim((-10, 10))
    plt.savefig(outfile + "foundEdgesDirectional_-10-10.svg", dpi=300)
    plt.savefig(outfile + "foundEdgesDirectional_-10-10.pdf", dpi=300)
    plt.close()
    
    plt.figure(figsize=(8,8), dpi=300)
    sns.lineplot(df, x='Distance Threshold', y='Found Change Points Ratio', hue='Tool', style='Tool')
    plt.title("Found Change Point Ratios vs Directional Distance Threshold", fontsize=18)
    plt.savefig(outfile + "foundEdgesRatioDirectional.svg", dpi=300)
    plt.savefig(outfile + "foundEdgesRatioDirectional.pdf", dpi=300)
    plt.xlim((-10, 10))
    plt.savefig(outfile + "foundEdgesRatioDirectional_-10-10.svg", dpi=300)
    plt.savefig(outfile + "foundEdgesRatioDirectional_-10-10.pdf", dpi=300)
    plt.close()

    # get counts for undirectional absolute threshold
    df["Absolute Distance"] = df["Distance Threshold"].abs()
    cumulativeDistance = df.groupby(["Absolute Distance", "Tool"], as_index=False)["Found Change Points"].sum()

    zero_counts = cumulativeDistance[cumulativeDistance['Absolute Distance'] == 0].set_index('Tool')['Found Change Points']

    cumulativeDistance['Found Change Points'] = cumulativeDistance.apply(
        lambda row: row['Found Change Points'] - zero_counts[row['Tool']] if row['Absolute Distance'] != 0 else zero_counts[row['Tool']],
        axis=1
    )
    cumulativeDistance["Found Change Points Ratio"] = cumulativeDistance["Found Change Points"] / df["Total Change Points"].values[0]

    # Set Tool column as categorical with the original order
    cumulativeDistance['Tool'] = pd.Categorical(cumulativeDistance['Tool'], categories=tool_order, ordered=True)

    # Sort grouped DataFrame by Tool to maintain the original order
    cumulativeDistance = cumulativeDistance.sort_values(by=['Tool', 'Absolute Distance']).reset_index(drop=True)

    cumulativeDistance.to_csv(outfile + "cumulativeFoundEdges.csv", index=False)

    plt.figure(figsize=(8,8), dpi=300)
    sns.lineplot(cumulativeDistance, x='Absolute Distance', y='Found Change Points Ratio', hue='Tool', style='Tool')
    plt.title("Found Change Point Ratios vs Absolute Distance Threshold", fontsize=18)
    plt.savefig(outfile + "foundEdgesRatio.svg", dpi=300)
    plt.savefig(outfile + "foundEdgesRatio.pdf", dpi=300)
    plt.xlim((0, 10))
    plt.savefig(outfile + "foundEdgesRatio_0-10.svg", dpi=300)
    plt.savefig(outfile + "foundEdgesRatio_0-10.pdf", dpi=300)
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

def getResults(args : Namespace, pool) -> dict:
    """
    Read all changepoint files and generate control set from a given bam file.
    
    Parameters
    ----------
    args : Namespace
        A Namespace object containing all the arguments from the command line.
    pool : multiprocessing.Pool
        A multiprocessing pool to parallelize the reading of changepoint files.
    
    Returns
    -------
    tuple[dict, np.ndarray]
        A tuple containing two elements: a dictionary containing all changepoint results and a numpy array containing the ground truth changepoints.
    """
    # gtReturn = pool.apply_async(readChangepoints, args=(args.changepoints,))
    # dynamontReturn = pool.apply_async(readDynamont, args=(args.dynamont[0],))
    dynamontReturn = {label : pool.apply_async(readDynamont, args=(result,)) for label, result in zip(args.labels, args.dynamont)}
    f5ceventalignReturn = pool.apply_async(readF5CEventalign, args=(args.f5ceventalign, os.path.splitext(args.f5ceventalign)[0] + '.sum'))
    f5cresquiggleReturn = pool.apply_async(readF5CResquiggle, args=(args.f5cresquiggle,))
    uncalled4return = pool.apply_async(readUncalled4, args=(args.uncalled,))
    doradoReturn = pool.apply_async(readDorado, args=(args.dorado,))
    controlReturn = pool.apply_async(generateControl, args=(args.basecalls,))

    # groundTruths = gtReturn.get()  #readChangepoints(args.changepoints)
    randomBorders, uniformBorders = controlReturn.get() # generating controls
    
    if args.tombo:
        tomboReturn = pool.apply_async(readTombo, args=(args.tombo,))
    
    toolsResult = {
        f"Dynamont {' '.join(label.split('_'))}" : dynamontReturn[label].get() for label in dynamontReturn
    } | {
        'f5c Eventalign' : f5ceventalignReturn.get(), # readF5CEventalign(args.f5ceventalign, os.path.splitext(args.f5ceventalign)[0] + '.sum'),
        'f5c Resquiggle' : f5cresquiggleReturn.get(), # readF5CResquiggle(args.f5cresquiggle),
        'Dorado' : doradoReturn.get(), # readDorado(args.dorado),
        'Uncalled4' : uncalled4return.get(), # readUncalled4(args.uncalled),
        'Control Random' : randomBorders,
        'Control Uniform' : uniformBorders,
    }

    if args.tombo:
        toolsResult = toolsResult | {
            'Tombo' : tomboReturn.get(), # readTombo(args.tombo),
        }

    for tool in toolsResult:
        toNumpy(toolsResult[tool])
    # toNumpy(groundTruths)

    return toolsResult #, groundTruths

def evaluateAgainstGroundTruth(groundTruths : np.ndarray, toolsResult : dict, maxDistance : int, output : str, pool) -> pd.DataFrame:

    print("Done\nStart Evaluating...")
    # df = pd.DataFrame()
    outfile = open(output, 'w')
    outfile.write("Tool,Distance Threshold,Found Change Points,Total Change Points,Found Change Points Ratio,Segmented Reads,Total Reads,Segmented Reads Ratio\n")
    outfile.flush()
    for tool, result in toolsResult.items():
        print(f"Evaluating {tool} on distance {maxDistance}...")

        totalEdges = 0 # total edges in ground truth
        foundEdges = np.zeros(2*maxDistance + 1, dtype=int)
        totalReads = 0 # total reads in ground truth
        segmentedReads = 0 # segmented reads by tool
    
        jobs = []
        for readid, gtSet in groundTruths.items():
            totalEdges += len(gtSet)
            totalReads += 1

            if readid in result and result[readid].size:
                predictions = result[readid]
                segmentedReads += 1
                jobs.append(pool.apply_async(evaluate, args=(gtSet, predictions, maxDistance)))

        [job.wait() for job in jobs]

        for i, job in enumerate(jobs):
            # if (i+1) % 1000 == 0:
            print(f'{i+1}/{len(groundTruths)}', end='\r')
            for j, n in enumerate(job.get()):
                foundEdges[j] += n
        
        for distance in range(2*maxDistance + 1):
            outfile.write(f"{tool},{distance - maxDistance},{foundEdges[distance]},{totalEdges},{foundEdges[distance] / totalEdges if totalEdges > 0 else 0},{segmentedReads},{totalReads},{segmentedReads / totalReads if totalReads > 0 else 0}\n")
        outfile.flush()
    
    outfile.close()

def plotSegmentSizeDistribution(toolsResult : dict, output : str) -> None:
    for tool in toolsResult:
        print("Plotting segment size distribution for " + tool)

        distance_list = [np.diff(toolsResult[tool][readid]) for readid in toolsResult[tool]]
        distances = np.concatenate(distance_list) if distance_list else np.array([])
        
        distances_clipped = np.clip(distances, 0, 200)

        sns.set_theme()
        sns.histplot(distances_clipped, bins=np.append(np.linspace(0, 200, 21), np.inf))
        plt.yscale("log")
        plt.xlim((0, 210))
        plt.xlabel("Segment Size")
        plt.title(tool + " Segmentsize Histogram, Overflow Bin for Values > 1000")
        plt.savefig(os.path.splitext(output)[0] + "_" + tool.replace(' ', '_') + "_segmentsize.pdf", dpi=300)
        plt.savefig(os.path.splitext(output)[0] + "_" + tool.replace(' ', '_') + "_segmentsize.svg", dpi=300)
        plt.close()
        sns.reset_defaults()

def createUpsetPlot(df, output):
    """
    Create an UpSet plot to visualize shared segment borders across tools.

    Parameters:
    - segmentCounts (dict): Dictionary where keys are tuples of tools and values are counts of shared segment borders.
    - toolNames (list): List of all tool names.
    - output (str): Output file prefix for saving the plot.
    """
    from upsetplot import UpSet

    try:
        # Create UpSet plot
        plt.figure(figsize=(10, 6))
        upset = UpSet(
            df,
            show_counts=True,
            # sort_by="cardinality",
            # max_degree=3,
            min_subset_size=500_000,
            shading_color=0.15,
            )
        upset.plot()

        # Customize count labels (shorten & rotate)
        ax = plt.gca()
        for text in ax.texts:
            try:
                num = int(text.get_text().replace(",", ""))  # Convert text to integer
                if num >= 1_000_000:  # Format large numbers
                    new_text = f"{num / 1_000_000:.1f}M"
                elif num >= 1_000:
                    new_text = f"{num / 1_000:.1f}K"
                else:
                    new_text = str(num)
                text.set_text(new_text)
                text.set_rotation(45)  # Rotate numbers
                text.set_fontsize(10)  # Adjust font size
            except ValueError:
                pass  # Skip non-numeric labels

        plt.title("Shared Segment Borders Across Tools (Min. 500K Segments)", fontsize=14)
        # plt.title("Shared Segment Borders Across Tools", fontsize=14)
        plt.tight_layout()

        # Save plot
        plt.savefig(f"{output}_upset_500k.pdf", dpi=300)
        plt.savefig(f"{output}_upset_500k.svg", dpi=300)
        plt.close()
    except IndexError:
        pass

    #! plot complete upset plot

    # Create UpSet plot
    plt.figure(figsize=(10, 6))
    upset = UpSet(
        df,
        show_counts=True,
        shading_color=0.15,
        )
    upset.plot()

    # Customize count labels (shorten & rotate)
    ax = plt.gca()
    for text in ax.texts:
        try:
            num = int(text.get_text().replace(",", ""))  # Convert text to integer
            if num >= 1_000_000:  # Format large numbers
                new_text = f"{num / 1_000_000:.1f}M"
            elif num >= 1_000:
                new_text = f"{num / 1_000:.1f}K"
            else:
                new_text = str(num)
            text.set_text(new_text)
            text.set_rotation(45)  # Rotate numbers
            text.set_fontsize(10)  # Adjust font size
        except ValueError:
            pass  # Skip non-numeric labels

    plt.title("Shared Segment Borders Across Tools", fontsize=14)
    # plt.title("Shared Segment Borders Across Tools", fontsize=14)
    plt.tight_layout()

    # Save plot
    plt.savefig(f"{output}_upset.pdf", dpi=300)
    plt.savefig(f"{output}_upset.svg", dpi=300)
    plt.close()

def processReadSegments(args):
    """
    Worker function to compute segment border overlaps per read.
    """
    readid, toolsResult, toolNames, pod5, windowSize = args
    signalLength = len(read5_ont.read(pod5).getSignal(readid))
    
    # Initialize an array to count segment borders at each position
    segmentPresence = np.zeros((len(toolNames), signalLength + 1), dtype=bool)

    # Mark positions where each tool has segment borders
    for i, tool in enumerate(toolNames):
        if readid in toolsResult[tool]:
            segmentPositions = toolsResult[tool][readid]
            segmentPresence[i, segmentPositions] = True  # Mark the entire window as True

    # Compute shared segment counts across tool combinations
    segmentCounts = {}
    for i in range(1, len(toolNames) + 1):
        for toolSubset in itertools.combinations(range(len(toolNames)), i):
            subsetPresence = segmentPresence[list(toolSubset)]  # Extract the subset tools

            # Instead of `np.all(subsetPresence, axis=0)`, check within window range
            sharedBorders = 0
            for pos in np.where(np.any(subsetPresence, axis=0))[0]:  # Only check marked positions
                # Check if all tools in the subset have a segment within the window
                if np.any(np.all(subsetPresence[:, max(0, pos - windowSize): min(signalLength, pos + windowSize + 1)], axis=0)):
                    sharedBorders += 1

            if sharedBorders > 0:
                segmentCounts[tuple(toolNames[j] for j in toolSubset)] = sharedBorders

    return segmentCounts

def processReadScores(args):
    """
    Worker function to compute segmentation scores per read.
    """
    readid, toolsResult, toolNames, pod5, window = args
    r5 = read5_ont.read(pod5)
    signal = r5.getSignal(readid)

    toolScores = {tool: [] for tool in toolNames}

    # Compute scores per tool
    for tool in toolNames:
        
        segPos = np.array(toolsResult[tool][readid])
        for i in range(1, len(segPos) - 1):
            start = max(0, segPos[i] - window)
            end = min(segPos[i] + window, len(signal))
            
            window1 = signal[start : segPos[i]]
            window2 = signal[segPos[i] : end]  # Next region
            medShift = abs(np.median(window2) - np.median(window1))
            madShift = abs(mad(window2) - mad(window1))
            s1 = medShift + madShift
            # s1 *= (1 - overSegPen) * (1 - underSegPen) # add weight to scores which penalize over- and undersegmentation

            segment = signal[segPos[i] : segPos[i+1]]
            if len(segment) >= 10:
                trim = max(int(0.1 * len(segment)), 1)  # 10% of segment length trimmed, but at least 1
                s2 = mad(segment[trim:-trim])
            else:
                s2 = np.nan

            toolScores[tool].append([s1, s2])

        toolScores[tool] = np.array(toolScores[tool])

    return toolScores

def scoreTools(toolsResult: dict, pod5: str, output: str, pool, window : int) -> None:
    
    if not os.path.exists(output + "_score.csv"):
        toolNames = list(toolsResult.keys())

        # Find reads segmented by all tools
        print("Finding reads segmented by all tools...")
        all_reads = list(set.intersection(*[set(toolsResult[tool].keys()) for tool in toolNames]))

        # Parallel processing of reads
        print("Start multiprocessing...")
        read_chunks = [(read, toolsResult, toolNames, pod5, window) for read in all_reads]
        toolScoresList = pool.map(processReadScores, read_chunks)

        # Aggregate scores
        print("Aggregating scores...")
        scoreSum = []
        for toolReadScores in toolScoresList:
            for tool, scores in toolReadScores.items():
                # if scores:  # Ensure non-empty scores
                scoreSum.extend(
                    [{"Tool": tool, "Score": s, "Segment Quality": "Contrast"} for s in scores[:, 0]] +
                    [{"Tool": tool, "Score": s, "Segment Quality": "Homogeneity"} for s in scores[:, 1]]
                )

        print("Computing statistics...")

        # Save results
        df = pd.DataFrame(scoreSum)
        df.to_csv(output + "_score.csv", index=False)

        print("Scoring complete. Results saved to", output + "_score.csv")

    else:

        df = pd.read_csv(output + "_score.csv", index_col = None)


    # Compute median and mean scores per tool
    median_contrasts = (
        df[df["Segment Quality"] == "Contrast"]
        .groupby("Tool")["Score"]
        .median()
        .sort_values(ascending=False)
    )

    mean_contrasts = (
        df[df["Segment Quality"] == "Contrast"]
        .groupby("Tool")["Score"]
        .mean()
        .reindex(median_contrasts.index)
    )

    median_homogeneity = (
        df[df["Segment Quality"] == "Homogeneity"]
        .groupby("Tool")["Score"]
        .median()
        .reindex(median_contrasts.index)
    )

    mean_homogeneity = (
        df[df["Segment Quality"] == "Homogeneity"]
        .groupby("Tool")["Score"]
        .mean()
        .reindex(median_contrasts.index)
    )


    print("Plotting scores")
    plt.figure(figsize=(12, 6))
    sns.set_theme(style="whitegrid")

    ax = sns.boxplot(
        data=df,
        x="Tool",
        y="Score",
        order=median_contrasts.index,  # Use sorted order
        hue="Segment Quality",
        legend=True,
        showfliers=True,
        flierprops={"marker": "o", "markersize": 3, "alpha": 0.3},  # Outlier marker
        notch=True,
        showcaps=False,
        width=0.6,  # Box width
        linewidth=1.5,  # Box border thickness
        )
    
    # Extract median values from the boxplot
    for i, tool in enumerate(median_contrasts.index):
        mc = median_contrasts[tool]  # Extract precomputed median
        mh = median_homogeneity[tool]

        ax.text(
            i - 0.15, mc, f"{mc:.2f}",
            ha="center", va="bottom", fontsize=10, fontweight="bold", color="black"
        )
        ax.text(
            i + 0.15, mh, f"{mh:.2f}",
            ha="center", va="bottom", fontsize=10, fontweight="bold", color="black"
        )

    plt.xticks(rotation=30, fontsize=12)  # Rotate tool names
    plt.yticks(fontsize=12)  # Adjust y-axis font size
    plt.title("Segmentation Tool Score Distribution", fontsize=16, fontweight="bold")
    plt.ylabel("Score", fontsize=14)
    plt.xlabel("Tool", fontsize=14)
    plt.grid(True, linestyle="--", alpha=0.5)
    ax.legend(loc="upper right", bbox_to_anchor=(1, 1))

    plt.tight_layout()

    ax.set_yscale("log")
    plt.savefig(output + "_score_log.pdf", dpi=300)
    plt.savefig(output + "_score_log.svg", dpi=300)
    plt.savefig(output + "_score_log.png", dpi=300)

    ax.set_yscale("linear")
    plt.ylim((0, 150))
    plt.savefig(output + "_score.pdf", dpi=300)
    plt.savefig(output + "_score.svg", dpi=300)
    plt.savefig(output + "_score.png", dpi=300)

    plt.close()
    sns.reset_defaults()

    print("Saving scores")
    # Create DataFrame
    score_summary = pd.DataFrame({
        "Tool": list(median_contrasts.index) * 2,  # Repeat tools for both segment qualities
        "Median Score": list(median_contrasts) + list(median_homogeneity),
        "Mean Score": list(mean_contrasts) + list(mean_homogeneity),
        "Segment Quality": ["Contrast"] * len(median_contrasts) + ["Homogeneity"] * len(median_homogeneity)
    })
    score_summary.to_csv(output + "_score.txt", index=False, sep="\t")

def compareTools(toolsResult: dict, pod5: str, output: str, pool, window : int) -> None:
    """
    Compare segment borders across tools and visualize them with a Venn diagram.

    Parameters:
    - toolsResult (dict): Dictionary containing segment borders per tool.
    - pod5 (str): Path to the POD5 file.
    - output (str): Path prefix for saving the Venn diagram.
    - pool (mp.Pool): Multiprocessing pool object.
    """

    # print("tools result", toolsResult.keys())

    if not os.path.exists(output + "_upset.csv"):
        toolNames = list(toolsResult.keys())

        # Find reads segmented by all tools
        print("Finding reads segmented by all tools...")
        all_reads = list(set.intersection(*[set(toolsResult[tool].keys()) for tool in toolNames]))

        # Parallel processing of reads
        print("Start multiprocessing...")
        read_chunks = [(read, toolsResult, toolNames, pod5, window) for read in all_reads]
        segmentCounts_list = pool.map(processReadSegments, read_chunks)

        # Merge segment counts from all workers
        print("Aggregating counts...")
        segmentCounts = {}
        for local_counts in segmentCounts_list:
            for key, value in local_counts.items():
                segmentCounts[key] = segmentCounts.get(key, 0) + value

        from upsetplot import from_memberships

        df = from_memberships([list(k) for k in segmentCounts.keys()], list(segmentCounts.values()))
        df.to_csv(output + "_upset.csv")
    else:
        print("Reading data from existing file...")
        df = pd.read_csv(output + "_upset.csv", index_col = None)
        df.rename(columns={df.columns[-1]: "count"}, inplace=True)
        toolNames = list(df.columns)[:-1]
        df[toolNames] = df[toolNames].astype(bool)
        df = df.set_index(toolNames).squeeze()

    print("Creating UpSet plot...")
    # if df.empty:
    #     print("No shared segment borders found for all tools.")
    # else:
    createUpsetPlot(df, output)
    

    # Plot Venn diagram
    # plt.figure(figsize=(10, 8))
    # venn(subsets=venn_labels, set_labels=toolNames)
    # plt.title("Shared Segment Coordinates Across Tools")
    # plt.savefig(output + "_venn.pdf", dpi=300)
    # plt.savefig(output + "_venn.svg", dpi=300)
    # plt.close()

def getReadStatsDict(basecalls: str) -> dict:
    """
    Extracts read lengths from a BAM file and returns them as a dictionary.

    Args:
        basecalls (str): Path to the BAM file containing read alignments.

    Returns:
        dict: A dictionary mapping read identifiers to their corresponding lengths.
              The read identifier is determined by the 'pi' tag if present, otherwise
              the query name is used.
    """

    import pysam
    readStats = {}
    with pysam.AlignmentFile(basecalls, "rb", check_sq=False) as bamfile:
        for read in bamfile.fetch(until_eof=True):
            readid = read.get_tag("pi") if read.has_tag("pi") else read.query_name
            readStats[readid] = (len(read.query_sequence), read.get_tag("qs"))
    return readStats

def plotSegmentationRate(toolsResult : dict, outfile : str) -> None:

    numReads = len(toolsResult["Dorado"])

    df = pd.DataFrame(
        {
            "Tool" : list(toolsResult.keys()),
            "Segmented Reads Ratio" : [len(value) / numReads for value in toolsResult.values()],
        }
    )

    plt.figure(figsize=(10,5), dpi=300)
    # sns.barplot(df.loc[df["Distance Threshold"] == 0], x='Tool', y='Segmented Reads Ratio', hue='Tool')
    sns.barplot(df, x='Tool', y='Segmented Reads Ratio', hue='Tool')
    # Annotate each bar with the height
    for index, (_, row) in enumerate(df.iterrows()):
        plt.text(x=index, y=row['Segmented Reads Ratio'],  # Add a small offset for better visibility
             s=f"{row['Segmented Reads Ratio']:.3f}",        # Format the text with 2 decimals
             ha='center',                     # Center the text horizontally
             va='bottom')                     # Align the text to the top of the bar
    # plt.ylim((0.8, 1.01))
    # plt.yticks(np.arange(0.8, 1.01, 0.05))
    plt.xticks(rotation = 65)
    plt.title("Ratio Of Segmented Reads Per Tool", fontsize=18)
    plt.tight_layout()
    plt.savefig(outfile + "_segmentedReadsRatio.svg", dpi=300)
    plt.savefig(outfile + "_segmentedReadsRatio.pdf", dpi=300)
    plt.close()
    sns.reset_defaults()

    df.to_csv(outfile + "_segmentedReadsRatio.csv", index=False)

def calcN50(readLength : np.ndarray) -> int:
    """
    Compute the N50 value from a list of read lengths.

    Parameters:
        read_lengths (list or np.array): List of sequencing read lengths.

    Returns:
        int: N50 value.
    """
    sortedLengths = np.sort(readLength)[::-1]
    totalLength = np.sum(sortedLengths)
    halfLength = totalLength / 2
    cumulLength = 0
    for length in sortedLengths:
        cumulLength += length
        if cumulLength >= halfLength:
            return length
    return 0

def plotReadStats(readLens : dict, output : str) -> None:

    data = np.array(list(readLens.values()))

    df = pd.DataFrame(
        {
            "Readlength" : data[:, 0],
            "Quality Score" : data[:, 1]
        }
    )

    sns.set_theme()
    sns.histplot(df, x = "Readlength", kde=True)
    plt.title("Read Length Distribution")
    plt.yscale('log')
    plt.tight_layout()
    plt.savefig(output + "_readLengths.svg", dpi=300)
    plt.savefig(output + "_readLengths.pdf", dpi=300)
    plt.close()
    sns.reset_defaults()

    sns.set_theme()
    sns.histplot(df, x = "Quality Score", kde=True)
    plt.title("Read Length Distribution")
    plt.yscale('log')
    plt.tight_layout()
    plt.savefig(output + "_readQuality.svg", dpi=300)
    plt.savefig(output + "_readQuality.pdf", dpi=300)
    plt.close()
    sns.reset_defaults()

    # Compute statistics
    n50 = calcN50(df["Readlength"].to_numpy())
    mean_qual = np.mean(df["Readlength"].to_numpy())
    median_qual = np.median(df["Readlength"].to_numpy())
    max_qual = np.max(df["Readlength"].to_numpy())
    min_qual = np.min(df["Readlength"].to_numpy())

    # Create a DataFrame
    df_stats = pd.DataFrame({
        "Metric": ["N50", "Mean", "Median", "Max", "Min"],
        "Value": [n50, mean_qual, median_qual, max_qual, min_qual]
    })

    # Save to CSV
    df_stats.to_csv(output + "_readLengths.csv", index=False, sep='\t')

    mean_qual = np.mean(df["Quality Score"].to_numpy())
    median_qual = np.median(df["Quality Score"].to_numpy())
    max_qual = np.max(df["Quality Score"].to_numpy())
    min_qual = np.min(df["Quality Score"].to_numpy())

    # Create a DataFrame
    df_stats = pd.DataFrame({
        "Metric": ["Mean", "Median", "Max", "Min"],
        "Value": [mean_qual, median_qual, max_qual, min_qual]
    })

    # Save to CSV
    df_stats.to_csv(output + "_readQuality.csv", index=False, sep='\t')

# def compute_log_likelihood(signal: np.ndarray, segmentation: list, kmer_means: dict, kmer_stds: dict) -> float:
#     """
#     Compute the log-likelihood score of the segmentation given the k-mer model.

#     Parameters:
#     - signal (np.ndarray): The observed signal values.
#     - segmentation (list): List of tuples (start, end, kmer), where each tuple defines a segment.
#     - kmer_means (dict): Dictionary mapping k-mers to their expected mean signal values.
#     - kmer_stds (dict): Dictionary mapping k-mers to their expected standard deviation.

#     Returns:
#     - float: Total log-likelihood score for the segmentation.
#     """
#     total_log_likelihood = 0.0
    
#     for start, end, kmer in segmentation:
#         segment_signal = signal[start:end]
        
#         if kmer not in kmer_means or kmer not in kmer_stds:
#             raise ValueError(f"K-mer {kmer} not found in the model.")
        
#         mean_k = kmer_means[kmer]
#         std_k = kmer_stds[kmer]
        
#         if std_k == 0:
#             std_k = 1e-6  # Prevent division by zero
        
#         # Compute Gaussian log-likelihood
#         log_likelihood = np.sum(
#             -0.5 * np.log(2 * np.pi * std_k**2) - ((segment_signal - mean_k) ** 2) / (2 * std_k**2)
#         )
        
#         total_log_likelihood += log_likelihood
    
#     return total_log_likelihood

def main() -> None:
    args = parse()
    print(args)

    pool = mp.Pool(args.processes)

    assert len(args.dynamont) == len(args.labels), "Number of dynamont results must match the number of labels"
    
    toolsResult = getResults(args, pool)
    
    plotSegmentationRate(toolsResult, os.path.splitext(args.output)[0])

    compareTools(toolsResult, args.pod5, os.path.splitext(args.output)[0], pool, args.window)

    readStats = getReadStatsDict(args.basecalls)

    plotReadStats(readStats, os.path.splitext(args.output)[0])

    scoreTools(toolsResult, args.pod5, os.path.splitext(args.output)[0], pool, 6)

    pool.close()
    pool.join()

if __name__ == '__main__':
    main()