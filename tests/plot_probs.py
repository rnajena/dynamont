#!/usr/bin/env python
# author: Jannes Spangenberg
# e-mail: jannes.spangenberg@uni-jena.de
# github: https://github.com/JannesSP
# website: https://jannessp.github.io

import os
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

def readM() -> np.ndarray:
    l = []
    with open(os.path.join(os.path.dirname(__file__), 'probs.txt'), 'r') as r:
        for line in r:
            line = line.strip()[:-1].split(',')
            l.append(list(map(float, line)))
    return np.exp(l)

def main() -> None:
    matrix = readM()
    sns.heatmap(matrix)
    plt.savefig(os.path.join(os.path.dirname(__file__), 'props_heatmaps.pdf'))
    
if __name__ == '__main__':
    main()
