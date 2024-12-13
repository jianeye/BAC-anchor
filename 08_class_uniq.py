#!/usr/bin/env python3
# Usage: python xx.py input.sam output_prefix
#python 07_class_uniq.py  poc.filter.novel.1.sam.uniq pob.filter.novel.1.sam.uniq
import sys
import re

# Function to process SAM file and classify reads into four categories
def classify_reads(input_file, output_prefix):
    reads_data = {}  # Store reads by their base names (_F/_R removed)

    # Open four output files for the different categories
    with open(f"{output_prefix}_diff_chr.txt", 'w') as cat1, \
         open(f"{output_prefix}_same_chr.txt", 'w') as cat2, \
         open(f"{output_prefix}_same_hap.txt", 'w') as cat3, \
         open(f"{output_prefix}_unpair.txt", 'w') as cat4, \
         open(input_file, 'r') as infile:

        for line in infile:
            if line.startswith('@'):  # Skip header lines
                continue

            columns = line.strip().split('\t')
            if len(columns) < 11:  # Ensure the SAM file has at least 11 fields
                continue

            read_name = columns[0]  # First column: read name
            chr_info = columns[2]   # Third column: chromosome info
            pos = columns[3]        # Fourth column: mapping position

            # Only process reads with _F or _R in their names
            if '_F' not in read_name and '_R' not in read_name:
                continue

            base_read_name = re.sub(r'_F$|_R$', '', read_name)

            # If this is the first time we've seen this base read name, store the whole line
            if base_read_name not in reads_data:
                reads_data[base_read_name] = [line.strip()]  # Store the entire line
            else:
                # We found a pair, compare the current read with the stored one
                first_line = reads_data[base_read_name][0]
                first_columns = first_line.split('\t')

                first_read = first_columns[0]
                first_chr_info = first_columns[2]
                first_pos = first_columns[3]

                # Parse chromosome and haplotype
                first_chr_split = first_chr_info.split('_')
                second_chr_split = chr_info.split('_')

                if len(first_chr_split) < 2 or len(second_chr_split) < 2:
                    continue  # Skip reads with invalid chr_info format

                first_chr_num = first_chr_split[0]
                second_chr_num = second_chr_split[0]
                first_haplotype = first_chr_split[1]
                second_haplotype = second_chr_split[1]

                # Calculate the position difference
                pos_diff = abs(int(first_pos) - int(pos))

                # Output to the appropriate category file
                if first_chr_num != second_chr_num:
                    # Category 1: Different chromosomes
                    cat1.write(f"{first_read}\t{first_chr_info}\t{first_pos} : {read_name}\t{chr_info}\t{pos}\n")
                elif first_haplotype != second_haplotype:
                    # Category 2: Same chromosome, different haplotype
                    cat2.write(f"{first_read}\t{first_chr_info}\t{first_pos} : {read_name}\t{chr_info}\t{pos}\n")
                else:
                    # Category 3: Same chromosome, same haplotype
                    cat3.write(f"{first_read}\t{first_chr_info}\t{first_pos} : {read_name}\t{chr_info}\t{pos}: {pos_diff}\n")

                # Remove the processed pair from the dictionary
                del reads_data[base_read_name]

        # Write unpaired reads to the unpair file
        if reads_data:  # Only write to the file if there are unpaired reads
            for base_read_name, read_info in reads_data.items():
                for read in read_info:
                    cat4.write(f"{read}\n")  # Output the entire line for unpaired reads

# Main function to handle command-line arguments
if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python xx.py input.sam output_prefix")
        sys.exit(1)

    input_file = sys.argv[1]
    output_prefix = sys.argv[2]

    classify_reads(input_file, output_prefix)
