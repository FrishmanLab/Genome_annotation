#!/usr/bin/env python

import sys
import pandas as pd
import re

def gtfseq_to_fasta(gtfseq_file):
    with open(gtfseq_file, 'r') as tsv:
        with open(gtfseq_file + '.fasta', 'w') as fasta:
            for line in tsv:
                columns = line.strip().split('\t')
                header = '_'.join(columns[:9])
                sequence = columns[9]
                fasta.write(f'>{header}\n{sequence}\n')


def main(argv):
        gtfseq_file = argv[0]
        gtfseq_to_fasta(gtfseq_file)

if __name__ == "__main__":
        main(sys.argv[1:])