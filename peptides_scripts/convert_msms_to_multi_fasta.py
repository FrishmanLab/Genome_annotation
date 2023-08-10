#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
takes msms.txt file and convert to multifasta file
the header of each peptide will be the raw file name + scan number + length + sequence
command: python extract_from_msms_files.py /path/to/msms.txt /path/to/out_unfiltered_peptides.faa
"""
import os
import pandas as pd
import re
import sys
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO


def get_peptides_from_msms(input_file, output_path):
    # Read the TSV file
    df = pd.read_csv(input_file, sep='\t', header=0)

    # Create a list of SeqRecord objects
    records = []
    for i, row in df.iterrows():
        header = row[0] + '-' + str(row[1]) + '-' + str(row[2]) + '-' + row[3] + '-' + str(row[4]) + '-' + str(row[5])
        sequence = row[3]
        seq = Seq(sequence)
        record = SeqRecord(seq, id=header, name=header, description='')
        records.append(record)

    # Write the FASTA file
    SeqIO.write(records, output_path, 'fasta')


def main(argv):
    input_file = argv[0]
    output_path = argv[1]
    get_peptides_from_msms(input_file, output_path)


if __name__ == "__main__":
   main(sys.argv[1:])