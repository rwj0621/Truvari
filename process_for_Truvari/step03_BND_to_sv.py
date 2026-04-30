import sys
import re
import os
import argparse

def parse_bnd_rigorous(chrom1, pos1, alt):
    """
    严谨判定 BND 方向逻辑：
    基于 VCF 4.2 规范，结合括号方向与碱基位置判定 SV 类型
    """
    # 1. 解析正则：提取括号、染色体、位置
    match = re.search(r'([\[\]])(.+?):(\d+)([\[\]])', alt)
    if not match:
        return None, None, None
    
    bracket1, chrom2, pos2, bracket2 = match.groups()
    pos2 = int(pos2)
    
    # 2. 跨染色体判定
    if chrom1 != chrom2:
        return "TRA", chrom2, pos2

    # 3. 判定 N (碱基) 的位置
    # 如果第一个字符是碱基且不是括号，说明 N 在左侧 (e.g., N[...[ )
    is_n_first = alt[0].upper() in "ATGCN"

    # 4. 严谨方向判定
    if bracket1 == '[' and bracket2 == '[':
        if is_n_first:
            return "DEL", chrom2, pos2  # N[...[ -> 3'接5' (缺失)
        else:
            return "INV", chrom2, pos2  # [...[N -> 5'接5' (倒位)
            
    elif bracket1 == ']' and bracket2 == ']':
        if is_n_first:
            return "INV", chrom2, pos2  # N]...] -> 3'接3' (倒位)
        else:
            return "DUP", chrom2, pos2  # ]...]N -> 5'接3' (重复)
            
    return "BND", chrom2, pos2

def process_vcf(input_path, output_path):
    print(f"[*] Reading: {input_path}")
    count = 0
    
    with open(input_path, 'r') as f, open(output_path, 'w') as out:
        for line in f:
            # 处理 Header
            if line.startswith('#'):
                # 注入必要的元数据描述，防止 SURVIVOR 报错
                if "##INFO=<ID=SVTYPE" in line:
                    out.write(line)
                    out.write('##INFO=<ID=SVLEN,Number=1,Type=Integer,Description="Difference in length between REF and ALT variants">\n')
                    out.write('##INFO=<ID=END,Number=1,Type=Integer,Description="End position of the variant described by this record">\n')
                else:
                    out.write(line)
                continue
            
            cols = line.split('\t')
            chrom1 = cols[0]
            pos1 = int(cols[1])
            alt = cols[4]
            info = cols[7]

            # 仅处理包含 SVTYPE=BND 的记录
            if "SVTYPE=BND" in info:
                svtype, chrom2, pos2 = parse_bnd_rigorous(chrom1, pos1, alt)
                
                if svtype:
                    # 更新 SVTYPE
                    new_info = info.replace("SVTYPE=BND", f"SVTYPE={svtype}")
                    
                    # 如果是同染色体变异，计算长度和 END
                    if svtype != "TRA" and svtype != "BND":
                        svlen = abs(pos2 - pos1)
                        end_pos = max(pos1, pos2)
                        # 确保 INFO 中包含 SVLEN 和 END
                        if "SVLEN=" not in new_info:
                            new_info += f";SVLEN={svlen}"
                        if "END=" not in new_info:
                            new_info += f";END={end_pos}"
                    
                    cols[7] = new_info
                    out.write("\t".join(cols))
                    count += 1
                else:
                    out.write(line)
            else:
                out.write(line)
    
    print(f"[+] Done! Converted {count} BND records to Symbolic SVs.")
    print(f"[+] Output saved to: {output_path}")

if __name__ == "__main__":
    # 创建命令行参数解析器
    parser = argparse.ArgumentParser(description="根据括号方向和碱基位置，将 BND 记录解析为具体的 SV 类型 (DEL/INV/DUP/TRA)。")
    
    # 添加参数：-i 代表 input，-o 代表 output
    parser.add_argument("-i", "--input", required=True, help="输入的 VCF 文件路径")
    parser.add_argument("-o", "--output", required=True, help="输出的 VCF 文件路径")

    args = parser.parse_args()

    # 检查输入文件是否存在
    if not os.path.exists(args.input):
        print(f"错误: 找不到输入文件 {args.input}")
        sys.exit(1)

    # 调用处理函数
    process_vcf(args.input, args.output)