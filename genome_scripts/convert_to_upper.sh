
#!/bin/bash
genome=/home/students/q.abbas/+proj-q.abbas/AT/0-raw_data/1-genomes/GCF_000001735.4_TAIR10.1_genomic_nuclear.fna

/home/students/q.abbas/Tools/gmes_linux_64/probuild --reformat_fasta \
--in $genome \
--out $genome.upper.fa \
--uppercase 1 --letters_per_line 60 --original
