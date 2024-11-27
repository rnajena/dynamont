#!/usr/bin/env python
# author: Jannes Spangenberg
# e-mail: jannes.spangenberg@uni-jena.de
# github: https://github.com/JannesSP
# website: https://jannessp.github.io

from argparse import ArgumentDefaultsHelpFormatter, ArgumentParser, Namespace
from src.python.segmentation.FileIO import plotParameters

def parse() -> Namespace:
    parser = ArgumentParser(
        formatter_class=ArgumentDefaultsHelpFormatter
    )
    parser.add_argument("params")
    parser.add_argument("outdir")
    return parser.parse_args()

def main() -> None:
    args = parse()
    plotParameters(args.params, args.outdir)

if __name__ == '__main__':
    main()