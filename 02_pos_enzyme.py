import argparse

def find_sequence_positions_per_fasta(file_path, sequence):
    sequences = []
    current_sequence = ''
    header = ''

    with open(file_path, 'r') as file:
        for line in file:
            line = line.strip()
            if line.startswith('>'):
                if current_sequence:
                    sequences.append((header, current_sequence))
                header = line
                current_sequence = ''
            else:
                current_sequence += line
        
        if current_sequence:
            sequences.append((header, current_sequence))
    
    return sequences

def filter_sequences(sequences, sequence):
    filtered_sequences = []
    
    for header, seq in sequences:
        positions = []
        start = 0
        while True:
            start = seq.find(sequence, start)
            if start == -1:
                break
            positions.append(start + 1)  # 1-based index
            start += len(sequence)
        
        # Check the two conditions
        if any(pos < 2 for pos in positions) or any(len(seq) - pos == 5 for pos in positions):
            filtered_sequences.append((header, seq, positions))
    
    return filtered_sequences

def main(file_path, enzyme):
    # Parse enzyme to get the sequence
    enzyme_name, sequence = enzyme.split(':')
    
    # Find sequences and filter
    sequences = find_sequence_positions_per_fasta(file_path, sequence)
    filtered_sequences = filter_sequences(sequences, sequence)

    for header, seq, positions in filtered_sequences:
        positions_str = ', '.join(map(str, positions))
        print("{}".format(header))
        print("Sequence length: {} : Positions of '{}': {}".format(len(seq), sequence, positions_str))

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Process FASTA file and find enzyme positions.")
    parser.add_argument('file_path', type=str, help='Path to the input FASTA file')
    parser.add_argument('-enzyme', type=str, required=True, help='Enzyme in the format EnzymeName:Sequence')

    args = parser.parse_args()
    main(args.file_path, args.enzyme)

