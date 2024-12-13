#!/usr/bin/env python
# Usage: python xx.py input.sam

import sys
import pysam
from collections import defaultdict

# Get command line arguments
input_sam = sys.argv[1]
uniq_map_output = sys.argv[1] + ".uniq"
multiple_map_output = sys.argv[1] + ".multiple"
unmap_output = sys.argv[1] + ".unmap"

# Initialize counters
uniq_map_count = 0
multiple_map_seq_count = 0  # To count unique sequences
unmap_count = 0

# Count occurrences of each read name
read_count = defaultdict(int)

# First pass to count read names and store the header
with pysam.AlignmentFile(input_sam, "r") as infile:
    header = infile.header
    for read in infile:
        if read.is_unmapped:
            unmap_count += 1
        else:
            read_count[read.query_name] += 1

# Create sets to keep track of which sequences have been seen
multiple_map_seq_set = set()

# Create lists to hold the full read lines
uniq_mapped_reads = []
multiple_mapped_reads = []
unmapped_reads = []

# Second pass to classify the reads
with pysam.AlignmentFile(input_sam, "r") as infile:
    for read in infile:
        if read.is_unmapped:
            unmapped_reads.append(read)
        else:
            if read_count[read.query_name] == 1:
                uniq_map_count += 1
                uniq_mapped_reads.append(read)
            else:
                # Add the sequence to the set only once
                if read.query_name not in multiple_map_seq_set:
                    multiple_map_seq_count += 1
                    multiple_map_seq_set.add(read.query_name)
                multiple_mapped_reads.append(read)

# Print the counts
print(f"Uniq mapped reads: {uniq_map_count}")
print(f"Multiple mapped sequences: {multiple_map_seq_count}")
print(f"Unmapped reads: {unmap_count}")

# Write unique mapped reads to the output file
with pysam.AlignmentFile(uniq_map_output, "w", header=header) as outfile:
    for read in uniq_mapped_reads:
        outfile.write(read)

# Write multiple mapped reads to the output file
with pysam.AlignmentFile(multiple_map_output, "w", header=header) as outfile:
    for read in multiple_mapped_reads:
        outfile.write(read)

# Write unmapped reads to the output file
with pysam.AlignmentFile(unmap_output, "w", header=header) as outfile:
    for read in unmapped_reads:
        outfile.write(read)

print(f"Unique mapped reads have been written to {uniq_map_output}")
print(f"Multiple mapped reads have been written to {multiple_map_output}")
print(f"Unmapped reads have been written to {unmap_output}")
