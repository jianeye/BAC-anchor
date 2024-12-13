import argparse

# Read fasta names into a set
def read_fasta_names(file):
    with open(file, 'r') as f:
        return set(line.strip() for line in f if line.strip() and not line.startswith('>'))

# Write fasta names to a file
def write_fasta_names(file, names):
    with open(file, 'w') as f:
        for name in sorted(names):
            f.write(f"{name}\n")
    print(f"Number of names written to {file}: {len(names)}")

def process_files(file_A, file_B, file_C, output_dir):
    # Read names from files A, B, C
    names_A = read_fasta_names(file_A)
    names_B = read_fasta_names(file_B)
    names_C = read_fasta_names(file_C)

    # First category: Intersection of B and C
    names_D = names_B.intersection(names_C)

    # Second category: Unique to B after excluding D
    names_B_unique = names_B.difference(names_D)

    # Third category: Unique to C after excluding D
    names_C_unique = names_C.difference(names_D)

    # Fourth category: Names in A but not in B, C, or D
    names_A_unique = names_A.difference(names_B.union(names_C).union(names_D))

    # Write results to files and print counts
    write_fasta_names(f'{output_dir}/D.name', names_D)
    write_fasta_names(f'{output_dir}/B_unique.name', names_B_unique)
    write_fasta_names(f'{output_dir}/C_unique.name', names_C_unique)
    write_fasta_names(f'{output_dir}/A_unique.name', names_A_unique)

    print("Processing complete, results written to", output_dir)

if __name__ == "__main__":
    # Command-line argument parsing
    parser = argparse.ArgumentParser(description="Process FASTA files and generate four result sets")
    parser.add_argument('input_A', type=str, help='Path to input file A')
    parser.add_argument('input_B', type=str, help='Path to input file B')
    parser.add_argument('input_C', type=str, help='Path to input file C')
    parser.add_argument('output_dir', type=str, help='Directory to store output files')

    # Parse arguments
    args = parser.parse_args()

    # Process the files
    process_files(args.input_A, args.input_B, args.input_C, args.output_dir)

