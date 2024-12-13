# Usage: python xx.py input.fa -enzyme BstBI

import os
import sys
import argparse

# Create command-line argument parser
parser = argparse.ArgumentParser(description="Extract the left end of the first enzyme site and the right end of the last enzyme site from each sequence.")
parser.add_argument("input_file", help="Input FASTA file")
parser.add_argument("-enzyme", required=True, choices=['BstBI', 'ClaI'], help="Specify enzyme name (BstBI or ClaI)")
args = parser.parse_args()

# Enzyme name to cutting site mapping
enzyme_map = {
    'BstBI': 'TTCGAA',
    'ClaI': 'ATCGAT'
}

# Get enzyme cutting sequence based on provided enzyme name
enzyme_sequence = enzyme_map[args.enzyme]

# Output file
output_file = f"{args.input_file}.split.fa"
with open(output_file, "w") as oup:

    # Function to read FASTA file
    def readFa(fa):
        with open(fa, 'r') as FA:
            seqName, seq = '', ''
            for line in FA:
                line = line.strip('\n')
                if line.startswith('>'):
                    if seqName:
                        yield (seqName, seq)
                    seqName = line[1:]
                    seq = ''
                else:
                    seq += line
            if seqName:
                yield (seqName, seq)

    # Process each sequence
    for seqName, seq in readFa(args.input_file):
        positions = []
        start = 0

        # Find all enzyme cutting sites
        while True:
            start = seq.find(enzyme_sequence, start)
            if start == -1:
                break
            positions.append(start)
            start += len(enzyme_sequence)

        if positions:
            # Extract the left end of the first enzyme site
            first_enzyme_start = positions[0]
            left_end = seq[:first_enzyme_start]
            
            # Extract the right end of the last enzyme site
            last_enzyme_end = positions[-1] + len(enzyme_sequence)
            right_end = seq[last_enzyme_end:]
            
            # Write to output file
            oup.write(f">{seqName}_F\n{left_end}\n>{seqName}_R\n{right_end}\n")

print(f"Results written to {output_file}")

