#!/usr/bin/env python3

from spectrum_io.raw.thermo_raw import ThermoRaw
from pathlib import Path
import sys

def convert(raw_directory):
    path = Path(raw_directory)
    raw_files = path.glob("*.raw")
    for raw_file in raw_files:
        ThermoRaw.convert_raw_mzml(raw_file)

def main(argv):
    raw_directory = argv[0]
    convert(raw_directory)
    
if __name__ == "__main__":
    main(sys.argv[1:])