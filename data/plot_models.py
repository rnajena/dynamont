#!/usr/bin/env python
# author: Jannes Spangenberg
# e-mail: jannes.spangenberg@uni-jena.de
# github: https://github.com/JannesSP
# website: https://jannessp.github.io

import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from itertools import repeat

def getParams(file : str) -> dict:
    ret = {}
    with open(file, 'r') as r:
        r.readline()
        for line in r:
            kmer, level_mean, level_stdv, sd_mean, sd_stdv, ig_lambda, weight = line.strip().split()
            ret[kmer] = (float(level_mean), float(level_stdv))
    return ret

def createDataFrame(models : dict, n : int = 1000) -> pd.DataFrame:
    df = pd.DataFrame(columns=['model', 'value'])
    for kmer in models:
        new_entries = pd.DataFrame({
            'model' : repeat(kmer, n),
            'simulated pA' : np.random.normal(*models[kmer], n)
        })
        df = pd.concat((df, new_entries), ignore_index=True)
    return df

def plot(data : pd.DataFrame, outfile : str) -> None:
    sns.set_theme(style="ticks")
    sns.histplot(
        data,
        x='simulated pA',
        hue='model',
        multiple="stack",
        # palette="light:m_r",
        edgecolor=".3",
        linewidth=.5,
        # log_scale=(False, True),
        legend=False,
    )
    plt.title("R9 stacked 5mer model distributions (1024 models)")
    plt.grid(True, 'both', 'both', linestyle='--', color='grey', alpha=.6)
    plt.savefig(outfile)

def main() -> None:
    models = os.path.join(os.path.dirname(__file__), "template_median69pA.model")
    outfile = os.path.join(os.path.dirname(__file__), "model_distributions.pdf")
    print("Loading models")
    model_params = getParams(models)
    print("Creating plot data")
    dataframe = createDataFrame(model_params)
    print("Plotting")
    plot(dataframe, outfile)
    print("Done")

if __name__ == '__main__':
    main()