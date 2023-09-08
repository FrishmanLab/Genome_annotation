#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Evaluation of gene prediction(s) based on reference annotation
Important: sequence ids should have the same names in both .gtf files
"""

import pandas as pd
import re
import sys
import csv


def process_headers(genome_file):
    # Input file name
    file_in = genome_file
    # Output file name
    file_out = genome_file + '.modified'

    with open(file_in, 'r') as f_in, open(file_out, 'w') as f_out:
        # Read each line of the input file
        for line in f_in:
            # Check if the line contains '>'
            if '>' in line:
                # Replace spaces with '_'
                line = line.replace(' ', '_')
            # Write the line to the output file
            f_out.write(line)



def main(argv):
    genome_file = argv[0]
    process_headers(genome_file)

if __name__ == "__main__":
    main(sys.argv[1:])