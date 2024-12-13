import pandas as pd
import sys

def main(inp1, inp2, oup_prefix):
    # 读取文件1，跳过头部行
    with open(inp1) as f:
        data1 = [line for line in f if not line.startswith('@')]
    
    # 将数据转换为 DataFrame
    data1 = pd.DataFrame([x.strip().split('\t') for x in data1])
    print("文件1的形状:", data1.shape)

    # 根据实际列数设置列名
    if data1.shape[1] >= 11:
        data1.columns = ['read', 'flag', 'chr', 'position', 'mapq', 'cigar', 'mrnm', 'mpos', 'isize', 'seq', 'qual'] + [f'opt{i}' for i in range(data1.shape[1] - 11)]
    else:
        print("文件1的列数不匹配:", data1.shape[1])
        return

    # 读取文件2
    data2 = pd.read_csv(inp2, sep='\t', header=None, on_bad_lines='skip')
    print("读取文件2的前几行：")
    print(data2.head())

    # 提取要去掉的reads名称
    reads_to_remove = data2[0].str.split('/').str[1]

    # 去掉已配对的reads
    filtered_data = data1[~data1['read'].isin(reads_to_remove)]
    print("过滤后的数据行数:", filtered_data.shape[0])

    # 结果列表
    same_hap_results = []
    same_chr_results = []
    diff_chr_results = []

    # 处理同一单倍型的结果
    for chr_name in filtered_data['chr'].unique():
        chr_data = filtered_data[filtered_data['chr'] == chr_name]
        paired_reads = chr_data[chr_data['read'].str.endswith('F')].merge(
            chr_data[chr_data['read'].str.endswith('R')],
            on='read', suffixes=('_F', '_R')
        )

        for _, row in paired_reads.iterrows():
            pos_distance = abs(int(row['position_R']) - int(row['position_F']))
            same_hap_results.append(f"{row['read_F']}\t{row['chr']}\t{row['position_F']} : {row['read_R']}\t{row['chr']}\t{row['position_R']} : {pos_distance}")

    # 处理同一染色体的结果
    for chr_name in filtered_data['chr'].unique():
        chr_data = filtered_data[filtered_data['chr'] == chr_name]
        for _, row in chr_data.iterrows():
            if row['read'].endswith('F'):
                read_F = row['read']
                pos_F = row['position']
                same_chr_results.append(f"{read_F}\t{row['chr']}\t{pos_F}")

    # 处理不同染色体的结果
    for _, row in filtered_data.iterrows():
        if not row['chr'].startswith('chr'):
            diff_chr_results.append(f"{row['read']}\t{row['chr']}\t{row['position']}")

    # 输出结果
    try:
        with open(f'{oup_prefix}_same_hap.txt', 'w') as f:
            for line in same_hap_results:
                f.write(line + '\n')

        with open(f'{oup_prefix}_same_chr.txt', 'w') as f:
            for line in same_chr_results:
                f.write(line + '\n')

        with open(f'{oup_prefix}_diff_chr.txt', 'w') as f:
            for line in diff_chr_results:
                f.write(line + '\n')

        print("分类完成，结果已保存。")
    except Exception as e:
        print("保存文件时出错:", e)

if __name__ == "__main__":
    if len(sys.argv) != 4:
        print("用法: python xx.py <inp1> <inp2> <oup.prefix>")
    else:
        main(sys.argv[1], sys.argv[2], sys.argv[3])
