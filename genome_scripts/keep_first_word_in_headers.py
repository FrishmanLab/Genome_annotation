#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
keep only the first word in headers "most likely NC_**"
"""

import pandas as pd
import re
import sys
import csv


def keep_first_word(input_file):
    output_file = input_file + '.modified'
    with open(input_file, 'r') as in_file, open(output_file, 'w') as out_file:
        for line in in_file:
            if line.startswith(">"):
                first_word = line.strip().split()[0][1:]
                out_file.write(">" + first_word + "\n")
            else:
                out_file.write(line)



def main(argv):
    genome_file = argv[0]
    keep_first_word(genome_file)

if __name__ == "__main__":
    main(sys.argv[1:])

