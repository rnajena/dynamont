#!/usr/bin/env python
# author: Jannes Spangenberg
# e-mail: jannes.spangenberg@uni-jena.de
# github: https://github.com/JannesSP
# website: https://jannessp.github.io

import sys
sys.path.append('/home/yi98suv/projects/dynamont/src')

from argparse import ArgumentDefaultsHelpFormatter, ArgumentParser, Namespace
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from os.path import exists, join
from os import makedirs
from dynamont.FileIO import getFiles, readKmerModels
from read5 import read
from scipy.stats import norm
from itertools import product

# STANDARDONTMODEL = pd.read_csv("/home/yi98suv/projects/dynamont/data/template_median69pA_extended.model", sep='\t', index_col = "kmer")
# STANDARDONTMODEL = readKmerModels("/data/fass5/jannes/dynamont/training_ecoli_wt_basic_norm_t24/_2024-06-28_14-31-53/trained_0_1000.model")
MODEL : pd.DataFrame = None

def parse() -> Namespace:
    parser = ArgumentParser(
        formatter_class=ArgumentDefaultsHelpFormatter
    )
    parser.add_argument('raw', type=str, help='Path to raw ONT training data')
    parser.add_argument('out', type=str, help='Output directory')
    return parser.parse_args()

def getPDFs(signal : np.ndarray, motif : str) -> np.ndarray:
    mean, stdev = MODEL[motif]
    return norm.pdf(signal, mean, stdev)

def iterateReads(rawFiles : list) -> pd.DataFrame:
    THRESHOLD = 100000000
    motifs = [''.join(motif) for motif in product('ACGT', repeat=5)]
    pdfs = np.array([])
    for file in rawFiles:
        r5 = read(file)
        for readid in r5.getReads():
            normSignal = r5.getZNormSignal(readid, "mean")
            for motif in motifs:
                pdfs = np.append(pdfs, getPDFs(normSignal, motif))
                if len(pdfs) > THRESHOLD:
                    print("Number of analysed points:", len(pdfs))
                    return pd.DataFrame({'PDF': pdfs})
                else:
                    print("Number of analysed points:", len(pdfs), end="\r")

def plot(pdfs : pd.DataFrame, outdir : str):
    print("Plotting")
    g = sns.histplot(pdfs, x="PDF", kde=True, bins=100)
    plt.savefig(join(outdir, "pdfs.svg"), dpi=300)
    plt.savefig(join(outdir, "pdfs.pdf"), dpi=300)
    g.set_yscale("log")
    plt.savefig(join(outdir, "pdfs_ylog.svg"), dpi=300)
    plt.savefig(join(outdir, "pdfs_ylog.pdf"), dpi=300)
    plt.close()

def main() -> None:
    global MODEL
    args = parse()
    # print(args)
    rawFiles = getFiles(args.raw, True)
    print(f'ONT Files: {len(rawFiles)}')
    MODEL = readKmerModels("/home/yi98suv/projects/dynamont/data/norm_models/rna_r9.4_180mv_70bps_extended_stdev1.model")
    if not exists(join(args.out, "stdev1")):
        makedirs(join(args.out, "stdev1"))
    plot(iterateReads(rawFiles), join(args.out, "stdev1"))

    MODEL = readKmerModels("/home/yi98suv/projects/dynamont/data/norm_models/rna_r9.4_180mv_70bps_extended_stdev0_5.model")
    if not exists(join(args.out, "stdev0_5")):
        makedirs(join(args.out, "stdev0_5"))
    plot(iterateReads(rawFiles), join(args.out, "stdev0_5"))

    MODEL = readKmerModels("/data/fass5/jannes/dynamont/training_ecoli_wt_basic_norm_t24/_2024-06-28_14-31-53/trained_0_1.model")
    if not exists(join(args.out, "trained")):
        makedirs(join(args.out, "trained"))
    plot(iterateReads(rawFiles), join(args.out, "trained"))

if __name__ == '__main__':
    main()