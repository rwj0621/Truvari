import sys
import argparse
import os

def process_vcf(input_vcf, output_vcf):
    seen_mates = set()
    count_removed = 0

    # 检查输入文件是否存在
    if not os.path.exists(input_vcf):
        print(f"错误: 找不到输入文件 {input_vcf}")
        return

    with open(input_vcf, 'r') as f, open(output_vcf, 'w') as out:
        for line in f:
            # 保留 VCF 头部信息
            if line.startswith('#'):
                out.write(line)
                continue
            
            cols = line.split('\t')
            # 基础校验，防止处理非 VCF 行
            if len(cols) < 8:
                continue

            current_id = cols[2]
            info = cols[7]
            
            # 提取 MATEID
            mate_id = ""
            for field in info.split(';'):
                if field.startswith('MATEID=') or field.startswith('MATE_ID='):
                    mate_id = field.split('=')[1]
                    break # 找到就停止当前循环
            
            # 核心去重逻辑
            if current_id in seen_mates:
                count_removed += 1
                continue
            else:
                if mate_id:
                    seen_mates.add(mate_id)
                out.write(line)

    print(f"--- 处理完成 ---")
    print(f"输入文件: {input_vcf}")
    print(f"输出文件: {output_vcf}")
    print(f"已移除冗余 BND 记录数: {count_removed}")

if __name__ == "__main__":
    # 创建命令行参数解析器
    parser = argparse.ArgumentParser(description="将 VCF 中的双条目 BND 转换为单条目表示。")
    
    # 添加参数：-i 代表 input，-o 代表 output
    parser.add_argument("-i", "--input", required=True, help="输入的 VCF 文件路径")
    parser.add_argument("-o", "--output", required=True, help="输出的 VCF 文件路径")

    args = parser.parse_args()

    # 调用处理函数
    process_vcf(args.input, args.output)
