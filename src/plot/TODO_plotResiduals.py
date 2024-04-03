#!/usr/bin/env python
# author: Jannes Spangenberg
# e-mail: jannes.spangenberg@uni-jena.de
# github: https://github.com/JannesSP
# website: https://jannessp.github.io

'''
Use the given segmentation and calculate the segment variance and the segmentwise signal distance to the model mean.
'''

from argparse import ArgumentDefaultsHelpFormatter, ArgumentParser, Namespace
import os
import h5py
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import scipy

def parse() -> Namespace:
    parser = ArgumentParser(
        formatter_class=ArgumentDefaultsHelpFormatter
    )
    #TODO insert arguments
    return parser.parse_args()

def main() -> None:
    args = parse()
    #TODO do stuff

if __name__ == '__main__':
    main()