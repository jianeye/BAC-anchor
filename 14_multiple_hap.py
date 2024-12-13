import re
import sys
from collections import defaultdict

# 输入文件路径和输出文件路径
input_file = sys.argv[1]  # pob.filter.novel.1.uniq_unpair_multiple_unused.sam  
output_file1 = input_file + ".mul.hap"
output_file2 = input_file + ".mul.hap.60"

# 用于存储每个read名称对应的ccs_F和ccs_R映射位置信息
read_dict = defaultdict(lambda: {"F": None, "R": None})

# 读取SAM文件
with open(input_file, "r") as infile:
    for line in infile:
        # 跳过SAM文件头信息
        if line.startswith("@"):
            continue
        
        columns = line.strip().split()
        read_name = columns[0]
        chrom = columns[2]
        pos = columns[3]
        
        # 提取read名称中的标记（_F 或 _R）
        match = re.search(r"(.+)(_F|_R)", read_name)
        if match:
            read_core = match.group(1)  # 取得read名称的主体部分
            direction = match.group(2)[1]  # 取得方向标记 'F' 或 'R'
            
            # 将chrom, pos, 和direction信息添加到字典中
            read_dict[read_core][direction] = (chrom, pos)

# 将符合条件的成对reads写入输出文件1
with open(output_file1, "w") as outfile1:
    for read_core, mappings in read_dict.items():
        # 确保同时存在F和R端
        if mappings["F"] and mappings["R"]:
            # 确保F和R端映射到相同的染色体（单倍型）
            if mappings["F"][0] == mappings["R"][0]:  # 判断染色体是否相同
                # 计算距离
                distance = abs(int(mappings['F'][1]) - int(mappings['R'][1]))
                
                # 格式化输出数据
                output_line = (
                    f"{read_core}_F\t{mappings['F'][0]}\t{mappings['F'][1]}:"
                    f"{read_core}_R\t{mappings['R'][0]}\t{mappings['R'][1]}:{distance}\n"
                )
                outfile1.write(output_line)

# 将符合条件的成对reads写入输出文件2
with open(output_file2, "w") as outfile2:
    for read_core, mappings in read_dict.items():
        # 确保同时存在F和R端
        if mappings["F"] and mappings["R"]:
            # 确保F和R端映射到相同的染色体（单倍型）
            if mappings["F"][0] == mappings["R"][0]:  # 判断染色体是否相同
                # 计算距离
                distance = abs(int(mappings['F'][1]) - int(mappings['R'][1]))
                
                # 判断距离条件
                if 60000 <= distance <= 100000:
                    output_line = (
                        f"{read_core}_F\t{mappings['F'][0]}\t{mappings['F'][1]}:"
                        f"{read_core}_R\t{mappings['R'][0]}\t{mappings['R'][1]}:{distance}\n"
                    )
                    outfile2.write(output_line)

print(f"已完成，符合条件的成对reads已写入 {output_file1} 和 {output_file2}")
