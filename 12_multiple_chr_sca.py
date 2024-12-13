#!/usr/bin/env python3

import argparse
from collections import defaultdict

def is_chromosome(rname):
    """
    判断 RNAME 是否对应于染色体。
    根据实际情况修改此函数，以匹配染色体的命名规则。
    """
    return rname.startswith('chr')

def parse_args():
    """
    解析命令行参数。
    """
    parser = argparse.ArgumentParser(
        description='筛选出同时比对到染色体和 scaffold 的读，并列出所有相关比对。'
    )
    parser.add_argument('input_sam', help='输入的多重比对 SAM 文件')
    parser.add_argument('output_sam', help='输出的筛选后 SAM 文件')
    return parser.parse_args()

def main():
    args = parse_args()
    input_sam = args.input_sam
    output_sam = args.output_sam

    # 字典存储每个读的所有比对信息
    read_mappings = defaultdict(list)  # QNAME -> list of SAM lines
    read_chromosome = defaultdict(bool)  # QNAME -> 是否比对到染色体
    read_scaffold = defaultdict(bool)  # QNAME -> 是否比对到 scaffold
    header_lines = []  # 存储头部信息

    print("读取输入 SAM 文件...")
    with open(input_sam, 'r') as infile:
        for line in infile:
            if line.startswith('@'):
                header_lines.append(line)
                continue
            fields = line.strip().split('\t')
            if len(fields) < 11:
                # 非标准的 SAM 对齐行，跳过
                continue
            qname = fields[0]
            rname = fields[2]
            read_mappings[qname].append(line)
            if is_chromosome(rname):
                read_chromosome[qname] = True
            else:
                read_scaffold[qname] = True

    print("筛选出同时比对到染色体和 scaffold 的读...")
    # 确定哪些读同时比对到了染色体和 scaffold
    reads_to_keep = set()
    for qname in read_mappings:
        if read_chromosome.get(qname, False) and read_scaffold.get(qname, False):
            reads_to_keep.add(qname)

    print(f"总读数: {len(read_mappings)}")
    print(f"同时比对到染色体和 scaffold 的读数: {len(reads_to_keep)}")

    print("写入筛选后的 SAM 文件...")
    with open(output_sam, 'w') as outfile:
        # 写入头部信息
        for header in header_lines:
            outfile.write(header)
        # 写入符合条件的所有比对信息
        for qname in reads_to_keep:
            for mapping_line in read_mappings[qname]:
                outfile.write(mapping_line)

    print(f"筛选完成。结果已写入到 '{output_sam}'。")

if __name__ == '__main__':
    main()

