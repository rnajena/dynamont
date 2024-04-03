#!/usr/bin/env python
# author: Jannes Spangenberg
# e-mail: jannes.spangenberg@uni-jena.de
# github: https://github.com/JannesSP
# website: https://jannessp.github.io

from argparse import ArgumentDefaultsHelpFormatter, ArgumentParser, Namespace
import pandas as pd

def parse() -> Namespace:
    parser = ArgumentParser(
        formatter_class=ArgumentDefaultsHelpFormatter
    )
    parser.add_argument("estimatedModels", help="Contains empty kmers (with nan values)")
    parser.add_argument("fillModels", help="Fully filled model file to fill missing values")
    return parser.parse_args()

def main() -> None:
    args = parse()

    estModels = pd.read_csv(args.estimatedModels, sep='\t')
    fillModels = pd.read_csv(args.fillModels, sep='\t')
    
    print(f"Filling {len(estModels[estModels['level_mean'].isna()]['kmer'].values)} missing models")

    for kmer in estModels[estModels['level_mean'].isna()]['kmer'].values:
        estModels.loc[estModels['kmer'] == kmer, 'level_mean'] = fillModels.loc[fillModels['kmer'] == kmer.replace('U', 'T'), 'level_mean']
        estModels.loc[estModels['kmer'] == kmer, 'level_stdv'] = fillModels.loc[fillModels['kmer'] == kmer.replace('U', 'T'), 'level_stdv']

    estModels.to_csv(args.estimatedModels + '.filled', index=False, sep='\t')

if __name__ == '__main__':
    main()