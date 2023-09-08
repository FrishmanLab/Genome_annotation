
#!/usr/bin/env python
import sys

def combine_sequnces(fasta_file):
    # combine sequences into one sequence
    with open(fasta_file, "r") as f_in, open(fasta_file + '.combined.fna', "w") as f_out:
        header = ""
        sequence = ""
        for line in f_in:
            line = line.strip()
            if line.startswith(">"):
                if header:
                    f_out.write(f"{header}\n{sequence}\n")
                header = line
                sequence = ""
            else:
                sequence += line
        f_out.write(f"{header}\n{sequence}\n")

def main(argv):
    fasta_file = argv[0]
    combine_sequnces(fasta_file)

if __name__ == "__main__":
    main(sys.argv[1:])

