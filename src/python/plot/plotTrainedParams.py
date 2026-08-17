#!/usr/bin/env python
# author: Jannes Spangenberg
# github: https://github.com/JannesSP
# website: https://jannessp.github.io

from argparse import ArgumentDefaultsHelpFormatter, ArgumentParser, Namespace
from python.segmentation.utils import plt_parameters

def parse() -> Namespace:
    parser = ArgumentParser(
        formatter_class=ArgumentDefaultsHelpFormatter
    )
    parser.add_argument("params")
    parser.add_argument("outdir")
    return parser.parse_args()

def main() -> None:
    args = parse()
    plt_parameters(args.params, args.outdir)

if __name__ == '__main__':
    main()