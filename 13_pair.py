#!/usr/bin/env python
# Usage: python xx.py inp oup1 oup2

import sys
import os

def process_sam_file(inp_file, oup1_file, oup2_file):
    reads_dict = {}

    # 检查输入文件是否存在
    if not os.path.isfile(inp_file):
        print(f"Error: Input file '{inp_file}' does not exist.")
        sys.exit(1)

    # 读取输入文件并将reads分类
    with open(inp_file, 'r') as inp:
        for line_num, line in enumerate(inp, 1):
            if line.startswith('@'):  # 跳过SAM header行
                continue
            fields = line.strip().split('\t')
            if len(fields) < 11:
                print(f"Warning: Line {line_num} does not have enough fields. Skipping.")
                continue
            read_name = fields[0]

            # 检查read名称是否以'ccs_F'或'ccs_R'结尾
            if read_name.endswith('ccs_F'):
                prefix = read_name[:-5]  # 去掉'ccs_F'
                suffix = 'F'
            elif read_name.endswith('ccs_R'):
                prefix = read_name[:-5]  # 去掉'ccs_R'
                suffix = 'R'
            else:
                # 如果不以'ccs_F'或'ccs_R'结尾，跳过
                continue

            if prefix not in reads_dict:
                reads_dict[prefix] = {'F': [], 'R': []}

            reads_dict[prefix][suffix].append(line)

    # 将处理后的结果写入输出文件
    paired_count = 0
    unpaired_count = 0

    with open(oup1_file, 'w') as oup1, open(oup2_file, 'w') as oup2:
        for prefix, pair in reads_dict.items():
            F_reads = pair['F']
            R_reads = pair['R']

            if F_reads and R_reads:
                paired_count += 1
                # 将所有F和R的组合写入oup1，F在前，R在后
                for F in F_reads:
                    oup1.write(F)
                for R in R_reads:
                    oup1.write(R)
            elif F_reads or R_reads:
                unpaired_count += 1
                # 只有F或R，写入oup2
                for F in F_reads:
                    oup2.write(F)
                for R in R_reads:
                    oup2.write(R)

    print(f"Processing complete.")
    print(f"Paired reads written to '{oup1_file}': {paired_count} pairs.")
    print(f"Unpaired reads written to '{oup2_file}': {unpaired_count} entries.")

if __name__ == "__main__":
    if len(sys.argv) != 4:
        print("Usage: python xx.py inp oup1 oup2")
        sys.exit(1)

    inp = sys.argv[1]
    oup1 = sys.argv[2]
    oup2 = sys.argv[3]

    process_sam_file(inp, oup1, oup2)
