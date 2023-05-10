#!/usr/bin/env python
# author: Jannes Spangenberg
# e-mail: jannes.spangenberg@uni-jena.de
# github: https://github.com/JannesSP
# website: https://jannessp.github.io

from os.path import join, dirname
from itertools import chain, combinations, product
import itertools
import pandas as pd

# TODO change this to the correct kmer size for the model before executing!
K = 5

BASES = ['A', 'C', 'G', 'T']
ALL_KMERS = set(map(lambda x: ''.join(x), product(BASES, repeat = K)))

def powerset(iterable):
    "powerset([1,2,3]) --> () (1,) (2,) (3,) (1,2) (1,3) (2,3) (1,2,3)"
    s = list(iterable)
    return list(chain.from_iterable(combinations(s, r) for r in range(len(s)+1)))[1:]

def main() -> None:
    df = pd.read_csv(join(dirname(__file__), 'template_median69pA.model'), sep='\t', header=0)
    add_on = pd.DataFrame(columns=df.columns)
    seen_patterns = set()
    N_pos_list = powerset(list(range(K)))

    # all possible positions for N
    for idx, N_set in enumerate(N_pos_list):

        print(f'{idx+1}/{len(N_pos_list)} patterns extended', end='\r')

        for kmer in ALL_KMERS:
            pattern = list(kmer)
            for pos in N_set:
                pattern[pos] = '.'
            pattern = ''.join(pattern)
            # skip seen patterns
            if pattern in seen_patterns:
                continue

            seen_patterns.add(pattern)
            subset = df[df.kmer.str.match(pattern)]
            level_mean = level_stdv = sd_mean = sd_stdv = ig_lambda = weight = 0
            for _, row in subset.iterrows():
                level_mean+=row['level_mean']
                level_stdv+=row['level_stdv']
                sd_mean+=row['sd_mean']
                sd_stdv+=row['sd_stdv']
                ig_lambda+=row['ig_lambda']
                weight+=row['weight']

            # setup new entry for pattern, simply use the mean of all value matching the pattern
            new_entry = {
                'kmer': [pattern.replace('.', 'N')],
                'level_mean': [level_mean/len(subset.index)],
                'level_stdv': [level_stdv/len(subset.index)],
                'sd_mean': [sd_mean/len(subset.index)],
                'sd_stdv': [sd_stdv/len(subset.index)],
                'ig_lambda': [ig_lambda/len(subset.index)],
                'weight': [weight/len(subset.index)],
                }
            add_on = pd.concat([add_on, pd.DataFrame.from_dict(new_entry)], ignore_index=True)

    print(f'{idx+1}/{len(N_pos_list)} patterns extended')
    new_model = pd.concat([df, add_on], ignore_index=True)
    fcols = new_model.select_dtypes('float').columns
    new_model[fcols] = new_model[fcols].apply(pd.to_numeric, downcast='float')
    new_model.to_csv(join(dirname(__file__), 'template_median69pA_extended.model'), index=False, sep='\t')

if __name__ == '__main__':
    main()