#!/usr/bin/env python
# author: Jannes Spangenberg
# e-mail: jannes.spangenberg@uni-jena.de
# github: https://github.com/JannesSP
# website: https://jannessp.github.io

from argparse import ArgumentDefaultsHelpFormatter, ArgumentParser, Namespace
import matplotlib.pyplot as plt
import read5_ont as r5
# from Bio import SeqIO
from os.path import splitext, join
import numpy as np
from matplotlib.patches import Rectangle

def parse() -> Namespace:
    parser = ArgumentParser(
        formatter_class=ArgumentDefaultsHelpFormatter
    )
    parser.add_argument("--dorado", type=str, required=True, metavar="TSV", help="Dorado segmentation output in TSV format")
    parser.add_argument("--f5cresquiggle", type=str, required=True, metavar="TSV", help="f5c resquiggle segmentation output in TSV format")
    parser.add_argument("--f5ceventalign", type=str, required=True, metavar="TSV", help="f5c eventalign segmentation output in TSV format\nSummary file must exists in the same path with extension .sum")
    parser.add_argument("--uncalled4", type=str, metavar="TSV", required=True, help="Uncalled4 segmentation output in csv format")
    parser.add_argument("--dynamont", type=str, required=True, metavar="CSV", help="Dynamont segmentation output in CSV format")
    # parser.add_argument("--tombo", type=str, metavar="PATH", help="Parent directory of single fast5s processed by Tombo")
    parser.add_argument("--basecalls", type=str, required=True, metavar="BAM", help="Basecalls of ONT training data as .bam file")
    parser.add_argument("--pod5", type=str, required=True, metavar="POD5", help="raw signals")
    parser.add_argument("-k", type=int, required=True, metavar="INT", help="kmer length")
    parser.add_argument("-o", "--out", type=str, required=True, metavar="PATH", help="Output path")
    parser.add_argument("--readid", type=str, required=True, metavar="ID", help="Read ID to plot")
    return parser.parse_args()

def readDynamont(file: str, readid : str) -> list:
    print("Reading Dynamont output from " + file)
    segments = []
    with open(file, 'r') as f:
        next(f) # skip header
        # readid,signalid,start,end,basepos,base,motif,state,posterior_probability,polish
        for line in f:
            if not line.startswith(readid):
                continue
            try:
                _, _, start, end, _, base, _, *_ = line.strip().split(',')
            except ValueError: # empty line in file
                pass
            segments.append((int(start), int(end), base))
    return segments

def readUncalled4(file : str, readid : str) -> list:
    print("Reading Uncalled4 output from " + file)
    segments = []
    with open(file, 'r') as f:
        next(f) # skip header
        # seq.name        seq.pos seq.strand      aln.id  seq.kmer        aln.read_id     dtw.start       dtw.length      dtw.current     dtw.current_sd
        for line in f:
            _, _, _, _, kmer, rid, start, length, *_ = line.strip().split('\t')
            if not rid == readid:
                continue
            if start == "*":
                continue
            segments.append((int(start), int(start) + int(length), kmer[len(kmer)//2]))
    return segments

def readF5CResquiggle(file : str, readid : str, read : str, halfKmer : int) -> list:
    print("Reading f5c Resquiggle output from " + file)
    segments = []
    with open(file, 'r') as f:
        next(f) # skip header
        # contig    position    reference_kmer  read_index  strand  event_index event_level_mean    event_stdv  event_length    model_kmer  model_mean  model_stdv  standardized_level  start_idx   end_idx
        for line in f:
            if not line.startswith(readid):
                continue
            _, kmer_idx, start, end = line.strip().split('\t')
            if start == '.' or end == '.':
                continue
            segments.append((int(start), int(end), read[int(kmer_idx) + halfKmer]))
    return segments

def readF5CEventalign(file: str, summary : str, readid : str) -> list:
    print("Reading f5c Eventalign output from " + file)
    readNum = None
    with open(summary, 'r') as f:
        # read_index    read_name   fast5_path  model_name  strand  num_events  num_steps   num_skips   num_stays   total_duration  shift   scale   drift   var
        for line in f:
            rid, readName, *_ = line.strip().split('\t')
            if not readName == readid:
                continue
            readNum = rid
            break

    segments = []
    with open(file, 'r') as f:
        next(f) # skip header
        # contig    position    reference_kmer  read_index  strand  event_index event_level_mean    event_stdv  event_length    model_kmer  model_mean  model_stdv  standardized_level  start_idx   end_idx
        for line in f:
            _, _, _, rNum, _, _, _, _, _, modelKmer, _, _, _, start, end = line.strip().split('\t')
            if not rNum == readNum:
                continue
            segments.append((int(start), int(end), modelKmer[len(modelKmer)//2]))
    return segments

def readDorado(file : str, readid : str) -> list:
    print("Reading Dorado output from " + file)
    segments = []
    with open(file, 'r') as f:
        next(f) # skip header
        # readid    position    base    motif   start   end
        for line in f:
            if not line.startswith(readid):
                continue
            rid, signalid, _, base, motif, start, end = line.strip().split('\t')
            segments.append((int(start), int(end), base))
    return segments

def readBasecalls(file : str, readid : str) -> str:
    import pysam
    with pysam.AlignmentFile(file, "rb", check_sq=False) as bamfile:
        for read in bamfile.fetch(until_eof=True):
            if not readid == read.query_name:
                continue
            return read.query_sequence

def main() -> None:
    basecolors = {
    'A':'#377eb8',
    'a':'#377eb8',
    'C':'#ff7f00',
    'c':'#ff7f00',
    'G':'#4daf4a',
    'g':'#4daf4a',
    'T':'#f781bf',
    't':'#f781bf',
    'U':'#f781bf',
    'u':'#f781bf',
    'N':'#ffffff',
    'n':'#ffffff',
    }
    
    args = parse()

    read = readBasecalls(args.basecalls, args.readid)

    tools = {
        "Dynamont" : np.array(readDynamont(args.dynamont, args.readid)),
        "Uncalled4" : np.array(readUncalled4(args.uncalled4, args.readid)),
        "f5c Resquiggle" : np.array(readF5CResquiggle(args.f5cresquiggle, args.readid, read, args.k // 2)),
        "f5c Eventalign" : np.array(readF5CEventalign(args.f5ceventalign, splitext(args.f5ceventalign)[0] + '.sum', args.readid)),
        # "Tombo" : readTombo(args.tombo, args.readid),
        "Dorado" : np.array(readDorado(args.dorado, args.readid)),
    }

    outpath = args.out
    signal = r5.read(args.pod5).getpASignal(args.readid)

    fig, ax = plt.subplots(nrows = 5, figsize=(110,15), dpi=300)
    fig.suptitle("Segmentation of the same read by different tools")
    fig.supylabel("Current (pA)")
    fig.supxlabel("Rel. Time (Sequencing Data Points)")

    for i, (tool, segments) in enumerate(tools.items()):
        segments = np.unique(segments, axis=1)
        # segments = segments[np.argsort(segments[:, 0].astype(int))]
        ax[i].plot(signal, linewidth=1.0, c='black')
        ax[i].set_xlim((0, len(signal)))
      #   ax[i].vlines(segments[:, 0], ymin=min(signal), ymax=max(signal), linewidth=0.5, linestyle='--', colors='red')
        ax[i].set_title(tool)
        # ax[i].legend()
        ax[i].set_xticks(np.arange(0, len(signal), 1000))

        for s in segments:
            start = int(s[0])
            end = int(s[1])
            base = str(s[2])
            # ax[i].text(start + (end - start)/2 - 6, -3.5, base, fontdict={'size' : 6, 'color':'black'}) # rotation=90
            ax[i].vlines([start, end], ymin=min(signal), ymax=max(signal), colors=basecolors[base], linestyles='--', linewidth=0.7)
            ax[i].add_patch(Rectangle((start, min(signal)), end - start, max(signal)-min(signal), alpha=0.4, edgecolor=basecolors[base], facecolor=basecolors[base]))

    plt.tight_layout()
    plt.savefig(join(outpath, args.readid + "_tool_segmentation.svg"), dpi=300)
    plt.savefig(join(outpath, args.readid + "_tool_segmentation.pdf"), dpi=300)
    plt.savefig(join(outpath, args.readid + "_tool_segmentation.png"), dpi=300)
    plt.close()

    # fig, ax = plt.subplots(nrows = 5, figsize=(20,10), dpi=300)
    # fig.suptitle(f"Segmentation of {args.readid} by Different Tools", fontsize=20)
    # fig.supylabel("Current (pA)", fontsize=18)
    # fig.supxlabel("Rel. Time (Sequencing Data Points)", fontsize=18)

    # for i, (tool, segments) in enumerate(tools.items()):
    #     segments = np.unique(segments, axis=1)
    #     # segments = segments[np.argsort(segments[:, 0].astype(int))]
    #     ax[i].plot(signal, linewidth=1.0, c='black')
    #   #   ax[i].vlines(segments[:, 0], ymin=min(signal), ymax=max(signal), linewidth=0.5, linestyle='--', colors='red')
    #     ax[i].set_title(tool, fontsize=18)
    #     # ax[i].legend()

    #     for s in segments:
    #         start = int(s[0])
    #         end = int(s[1])
    #         base = str(s[2])
    #         # ax[i].text(start + (end - start)/2 - 6, -3.5, base, fontdict={'size' : 6, 'color':'black'}) # rotation=90
    #         ax[i].vlines([start], ymin=min(signal), ymax=max(signal), colors=basecolors[base], linestyles='--', linewidth=0.7)
    #         ax[i].add_patch(Rectangle((start, min(signal)), end - start, max(signal)-min(signal), alpha=0.4, edgecolor=basecolors[base], facecolor=basecolors[base]))

    #     ax[i].set_xticks(np.arange(0, len(signal), 500))
    #     ax[i].set_xlim((19000, 22000))

    # plt.tight_layout()
    # plt.savefig(join(outpath, args.readid + "_tool_segmentation_region.svg"), dpi=300)

if __name__ == '__main__':
    main()