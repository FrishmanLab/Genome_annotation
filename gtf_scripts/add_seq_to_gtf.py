#!/usr/bin/env python

import sys
import pandas as pd
from Bio import SeqIO
def add_seqs(input_genome_path, input_gtf_path):
    column_names = ['Seqid', 'Source', 'Type', 'Start', 'End', 'Score', 'Strand', 'Frame', 'Suppl']
    df = pd.read_csv(input_gtf_path, sep='\t', header=None, names=column_names, dtype={'Start': int, 'End': int})
    df['Seq'] = ""

    genome_sequences = SeqIO.to_dict(SeqIO.parse(input_genome_path, "fasta"))

    for index, row in df.iterrows():
        seqid = row['Seqid']
        strand = row['Strand']
        start = row['Start'] - 1
        end = row['End']
        seq = str(genome_sequences[seqid].seq[start:end])
        df.at[index, 'Seq'] = seq

    df.to_csv(input_gtf_path + '+seqs.gtf', sep='\t', header=False, index=False)
    return df

def main(argv):
    input_genome = argv[0]
    input_gtf = argv[1]
    add_seqs(input_genome, input_gtf)

if __name__ == "__main__":
    main(sys.argv[1:])



