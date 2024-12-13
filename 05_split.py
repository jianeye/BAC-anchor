# Usage: python xx.py input.fa -enzyme BstBI (or ClaI)

import os
import sys
import argparse

# Create command-line argument parser
parser = argparse.ArgumentParser(description="Split sequences based on enzyme cutting site.")
parser.add_argument("input_file", help="Input FASTA file")
parser.add_argument("-enzyme", required=True, choices=['BstBI', 'ClaI'], help="Specify enzyme name (BstBI or ClaI)")
args = parser.parse_args()

# Enzyme name to cutting site mapping
enzyme_map = {
    'BstBI': 'TTCGAA',
    'ClaI': 'ATCGAT'
}

# Get command-line arguments
fa = args.input_file
enzyme_name = args.enzyme
enzyme_sequence = enzyme_map[enzyme_name]

# Output file
oup = open(f"{fa}.split.fa", "w")

# Function to read FASTA file
def readFa(fa):
    with open(fa, 'r') as FA:
        seqName, seq = '', ''
        while True:
            line = FA.readline()
            line = line.strip('\n')
            if (line.startswith('>') or not line) and seqName:
                yield (seqName, seq)
            if line.startswith('>'):
                seqName = line[1:]
                seq = ''
            else:
                seq += line
            if not line:
                break

# Process FASTA file based on enzyme cutting site
for seqName, seq in readFa(fa):
    seqName1 = list(seqName)
    seqName1.append("_F")
    F = "".join(seqName1)
    
    seqName2 = list(seqName)
    seqName2.append("_R")
    R = "".join(seqName2)
    
    if enzyme_sequence in seq:
        i1 = seq.split(enzyme_sequence)[0]
        i2 = seq.split(enzyme_sequence)[1]
    else:
        i1 = seq
        i2 = ''
    
    oup.write(f">{F}\n{i1}\n>{R}\n{i2}\n")

# Close output file
oup.close()
