#python 01_filter_enzyme_site.fa.py -enzy ATCGAT  poc.filter.novel.fa  poc.filter.novel.1.fa poc.filter.novel.more1.fa poc.filter.novel.count


import os
import sys
import argparse
import pyfastx

# Create a command-line argument parser
parser = argparse.ArgumentParser(description="Count enzyme cutting sites and output related information")
parser.add_argument("input_file", help="Input FASTA file")
parser.add_argument("output_file1", help="Output file for sequences with 1 enzyme site")
parser.add_argument("output_file2", help="Output file for sequences with more than 1 enzyme site")
parser.add_argument("output_file3", help="Output file for counting the number of sequences with different enzyme site counts")
parser.add_argument("-enzy", "--enzyme", required=True, help="Specify enzyme sequence, e.g., TTCGAA, claI: ATCGAT")

# Parse command-line arguments
args = parser.parse_args()

# Read command-line arguments
enzyme_sequence = args.enzyme
input_file = args.input_file
output_file1 = args.output_file1
output_file2 = args.output_file2
output_file3 = args.output_file3

# Prepare the command-line string as a comment
command_line = ' '.join(sys.argv)

# Process FASTA file
fa = pyfastx.Fastx(input_file)

# Function to write the command-line comment and then the data
def write_with_comment(file_name, data):
    with open(file_name, 'w') as f:
        f.write(f"# Command line: {command_line}\n")
        f.write(data)

# Collect data for each output file
data_file1 = ""
data_file2 = ""
data_file3 = ""

enzyme_count = {}

for record in fa:
    name = record[0]
    seq = record[1]
    n = seq.count(enzyme_sequence)

    # Collect data based on enzyme site count
    if n == 1:
        data_file1 += f">{name}\n{seq}\n"
    elif n > 1:
        data_file2 += f">{name}\n{seq}\n"  # Change to FASTA format

    # Count enzyme site occurrences
    if n in enzyme_count:
        enzyme_count[n] += 1
    else:
        enzyme_count[n] = 1

# Output the enzyme site count statistics to data_file3
for count, num_sequences in sorted(enzyme_count.items()):
    data_file3 += f"{count} enzyme site(s): {num_sequences} sequences\n"

# Write data to output files with the command-line comment
write_with_comment(output_file1, data_file1)
write_with_comment(output_file2, data_file2)
write_with_comment(output_file3, data_file3)
