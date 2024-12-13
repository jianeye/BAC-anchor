#!/usr/bin/env python3
###python 11_right uniq_map.file multiple_map.file oup_prefix
###python 11_uniq_multiple.py  pob.filter.novel.1.sam.uniq_unpair.txt  pob.filter.novel.1.sam.multiple pob.filter.novel.1.uniq_unpair_multiple

import sys
import argparse
from collections import defaultdict

def parse_args():
    parser = argparse.ArgumentParser(description='Categorize read pairs based on haplotype and chromosome mappings.')
    parser.add_argument('uniq_map', help='Unique mapping SAM file (no header)')
    parser.add_argument('multiple_map', help='Multiple mapping SAM file (with header)')
    parser.add_argument('output_prefix', help='Output prefix for the result files')
    parser.add_argument('--debug', action='store_true', help='Enable debug mode with detailed output')
    return parser.parse_args()

def parse_sam_line(line):
    fields = line.strip().split('\t')
    if len(fields) < 4:
        return None, None, None
    qname = fields[0]
    rname = fields[2]
    pos = fields[3]
    return qname, rname, pos

def get_pair_read_name(qname):
    # 将 '_F' 替换为 '_R'，反之亦然
    if qname.endswith('_F'):
        return qname[:-2] + '_R'
    elif qname.endswith('_R'):
        return qname[:-2] + '_F'
    else:
        return None  # 如果没有 '_F' 或 '_R' 后缀

def extract_haplotype(rname):
    # 假设单倍型信息总是最后一个下划线后的部分
    if '_' in rname:
        parts = rname.rsplit('_', 1)
        return parts[0], parts[1]  # chromosome, haplotype
    return rname, '0'  # default haplotype if not present

def read_unique_map(uniq_map_file, debug=False):
    read_pairs = {}
    with open(uniq_map_file, 'r') as f:
        for line_num, line in enumerate(f, 1):
            if not line.strip():
                continue
            qname, rname, pos = parse_sam_line(line)
            if qname is None:
                if debug:
                    print(f"Skipping malformed line {line_num} in unique_map.sam")
                continue
            pair_read = get_pair_read_name(qname)
            if pair_read is None:
                if debug:
                    print(f"Read name '{qname}' does not end with '_F' or '_R'. Skipping.")
                continue
            read_pairs[qname] = {
                'pair_read': pair_read,
                'rname': rname,
                'pos': pos
            }
            if debug and line_num <= 10:
                print(f"[Unique Map] Line {line_num}: QNAME={qname}, Pair_Read={pair_read}, RNAME={rname}, POS={pos}")
    return read_pairs

def read_multiple_map(multiple_map_file, debug=False):
    multiple_mappings = defaultdict(list)
    mapping_entries = []
    headers = []
    with open(multiple_map_file, 'r') as f:
        for line_num, line in enumerate(f, 1):
            if line.startswith('@'):
                headers.append(line.rstrip('\n'))
                continue  # skip header
            qname, rname, pos = parse_sam_line(line)
            if qname is None:
                if debug:
                    print(f"Skipping malformed line {line_num} in multiple_map.sam")
                continue
            multiple_mappings[qname].append((rname, pos, line.rstrip('\n')))
            mapping_entries.append(line.rstrip('\n'))
            if debug and line_num <= 10:
                print(f"[Multiple Map] Line {line_num}: QNAME={qname}, RNAME={rname}, POS={pos}")
    return multiple_mappings, mapping_entries, headers

def categorize_reads(read_pairs, multiple_mappings, debug=False):
    # Prepare output data structures
    same_hap = []
    dif_hap = []
    diff_chr = []
    else_best = []

    # Counters for debugging
    both_mapped = 0
    only_f_mapped = 0
    only_r_mapped = 0

    # Set to track used mappings (full_line strings)
    used_mappings = set()

    for base_num, (qname, info) in enumerate(read_pairs.items(), 1):
        pair_read = info['pair_read']
        read_rname = info['rname']
        read_pos = info['pos']

        if debug and base_num <= 10:
            print(f"\n[Categorize] Processing read pair {base_num}: {qname} / {pair_read}")

        # Check if pair_read exists in multiple_mappings
        if pair_read in multiple_mappings:
            pair_mappings = multiple_mappings[pair_read]
            both_mapped += 1
            if debug and base_num <= 10:
                print(f"  Found {len(pair_mappings)} mappings for pair read '{pair_read}'")

            # Extract haplotype information from unique map read
            chr_f_unique, hap_f_unique = extract_haplotype(read_rname)

            # Categorize based on pair_read's multiple mappings
            categorized = False
            for rname, pos, full_line in pair_mappings:
                chr_map, hap_map = extract_haplotype(rname)
                if chr_map == chr_f_unique and hap_map == hap_f_unique:
                    # Same haplotype
                    distance = abs(int(pos) - int(read_pos))
                    same_hap.append(f"{pair_read}\t{rname}\t{pos} : {qname}\t{read_rname}\t{read_pos} : {distance}")
                    used_mappings.add(full_line)
                    categorized = True
                    if debug and base_num <= 10:
                        print(f"  Assigned to same_hap: {pair_read} {rname} {pos} : {qname} {read_rname} {read_pos} : {distance}")
                    break  # Prioritize first match

            if not categorized:
                # Check for different haplotype but same chromosome
                for rname, pos, full_line in pair_mappings:
                    chr_map, hap_map = extract_haplotype(rname)
                    if chr_map == chr_f_unique and hap_map != hap_f_unique:
                        #distance = abs(int(pos) - int(read_pos))
                        dif_hap.append(f"{pair_read}\t{rname}\t{pos} : {qname}\t{read_rname}\t{read_pos}")
                        used_mappings.add(full_line)
                        categorized = True
                        if debug and base_num <= 10:
                            print(f"  Assigned to dif_hap: {pair_read} {rname} {pos} : {qname} {read_rname} {read_pos} : {distance}")
                        break

            if not categorized:
                # Different chromosome
                for rname, pos, full_line in pair_mappings:
                    chr_map, hap_map = extract_haplotype(rname)
                    if chr_map != chr_f_unique:
                        #distance = abs(int(pos) - int(read_pos))
                        diff_chr.append(f"{pair_read}\t{rname}\t{pos} : {qname}\t{read_rname}\t{read_pos}")
                        used_mappings.add(full_line)
                        categorized = True
                        if debug and base_num <= 10:
                            print(f"  Assigned to diff_chr: {pair_read} {rname} {pos} : {qname} {read_rname} {read_pos} : {distance}")
                        break

            if not categorized:
                # Else, assign to else_best with the first mapping
                rname, pos, full_line = pair_mappings[0]
                #distance = abs(int(pos) - int(read_pos))
                else_best.append(f"{pair_read}\t{rname}\t{pos} : {qname}\t{read_rname}\t{read_pos} ")
                used_mappings.add(full_line)
                if debug and base_num <= 10:
                    print(f"  Assigned to else_best (best mapping): {pair_read} {rname} {pos} : {qname} {read_rname} {read_pos} : {distance}")

        else:
            # Pair read not found in multiple_mappings
            else_best.append(f"{pair_read}: {qname}\t{read_rname}\t{read_pos}")
            if debug and base_num <= 10:
                print(f"  Assigned to else_best (pair read not found): {pair_read}  : {qname}\t{read_rname}\t{read_pos}")

    return same_hap, dif_hap, diff_chr, else_best, both_mapped, used_mappings

def write_output(same_hap, dif_hap, diff_chr, else_best, prefix, debug=False):
    with open(f"{prefix}_same_hap.sam", 'w') as f1, \
         open(f"{prefix}_dif_hap.sam", 'w') as f2, \
         open(f"{prefix}_diff_chr.sam", 'w') as f3, \
         open(f"{prefix}_else.sam", 'w') as f4:
        
        for line in same_hap:
            f1.write(line + '\n')
        for line in dif_hap:
            f2.write(line + '\n')
        for line in diff_chr:
            f3.write(line + '\n')
        for line in else_best:
            f4.write(line + '\n')
    if debug:
        print(f"[Write Output] Written {len(same_hap)} to {prefix}_same_hap.sam")
        print(f"[Write Output] Written {len(dif_hap)} to {prefix}_dif_hap.sam")
        print(f"[Write Output] Written {len(diff_chr)} to {prefix}_diff_chr.sam")
        print(f"[Write Output] Written {len(else_best)} to {prefix}_else.sam")

def write_unused_mappings(unused_mappings, headers, prefix, debug=False):
    with open(f"{prefix}_unused.sam", 'w') as f:
        for header in headers:
            f.write(header + '\n')
        for line in unused_mappings:
            f.write(line + '\n')
    if debug:
        print(f"[Write Output] Written {len(unused_mappings)} unused mappings to {prefix}_unused.sam")

def main():
    args = parse_args()

    if args.debug:
        print("[Main] Debug mode enabled.")

    print("Reading unique map...")
    read_pairs = read_unique_map(args.uniq_map, debug=args.debug)
    print(f"Total unique read pairs: {len(read_pairs)}")

    print("Reading multiple map...")
    multiple_mappings, mapping_entries, headers = read_multiple_map(args.multiple_map, debug=args.debug)
    print(f"Total reads with multiple mappings: {len(multiple_mappings)}")

    print("Categorizing reads...")
    same_hap, dif_hap, diff_chr, else_best, both_mapped, used_mappings = categorize_reads(read_pairs, multiple_mappings, debug=args.debug)
    print(f"Same haplotype: {len(same_hap)}")
    print(f"Different haplotype: {len(dif_hap)}")
    print(f"Different chromosome: {len(diff_chr)}")
    print(f"Else: {len(else_best)}")
    print(f"Read pairs with both 'F' and 'R' mapped: {both_mapped}")

    print("Writing output files...")
    write_output(same_hap, dif_hap, diff_chr, else_best, args.output_prefix, debug=args.debug)

    # Determine unused mappings
    print("Determining unused mappings...")
    used_mappings_set = set(used_mappings)
    unused_mappings = [line for line in mapping_entries if line not in used_mappings_set]
    write_unused_mappings(unused_mappings, headers, args.output_prefix, debug=args.debug)
    print("Done.")

if __name__ == "__main__":
    main()

