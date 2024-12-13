import sys

def filter_sam(input_file):
    with open(input_file, 'r') as f:
        lines = f.readlines()

    hap_file = open(input_file + '.hap', 'w')
    chr_file = open(input_file + '.chr', 'w')
    diff_file = open(input_file + '.diff', 'w')

    filtered_reads = {}

    for line in lines:
        if line.startswith('@'):
            continue
        
        fields = line.split('\t')
        read_name = fields[0]
        chr_name = fields[2]
        pos = fields[3]

        if chr_name == '*':
            continue

        haplo_type = chr_name  # 假设单倍型就是染色体名称

        base_read_name = read_name.rsplit('/', 1)[0]

        if base_read_name not in filtered_reads:
            filtered_reads[base_read_name] = {}
        
        filtered_reads[base_read_name][read_name] = (chr_name, pos, haplo_type)

    for base_read_name, mappings in filtered_reads.items():
        if len(mappings) == 2:  # 只有一对读名
            read_F = f"{base_read_name}/ccs_F"
            read_R = f"{base_read_name}/ccs_R"

            if read_F in mappings and read_R in mappings:
                chr_F, pos_F, haplo_F = mappings[read_F]
                chr_R, pos_R, haplo_R = mappings[read_R]

                if haplo_F == haplo_R:  # 同一单倍型
                    output_line = f"{read_F}\t{chr_F}\t{pos_F}: {read_R}\t{chr_R}\t{pos_R}\n"
                    chr_file.write(output_line)
                elif chr_F == chr_R:  # 同一染色体但不同单倍型
                    diff_file.write(f"{read_F}\t{chr_F}\t{pos_F} (不同单倍型)\n")
                    diff_file.write(f"{read_R}\t{chr_R}\t{pos_R} (不同单倍型)\n")
                else:  # 不同染色体
                    diff_file.write(f"{read_F}\t{chr_F}\t{pos_F} (不同染色体)\n")
                    diff_file.write(f"{read_R}\t{chr_R}\t{pos_R} (不同染色体)\n")
        elif len(mappings) > 2:  # 多个比对
            for read_name, (chr_name, pos, haplo_type) in mappings.items():
                diff_file.write(f"{read_name}\t{chr_name}\t{pos}\n")

    hap_file.close()
    chr_file.close()
    diff_file.close()

if __name__ == '__main__':
    if len(sys.argv) != 2:
        print("Usage: python filter_sam.py <input.sam>")
        sys.exit(1)
    
    input_file = sys.argv[1]
    filter_sam(input_file)
