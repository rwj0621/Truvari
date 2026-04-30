import sys
import os
import re
import argparse

# 定义标准染色体集合
STANDARD_CHROMS = set([f"chr{i}" for i in range(1, 23)] + ["chrX", "chrY"])

def process_vcf(input_vcf, output_dir):
    print(f"[*] 正在读取输入文件: {input_vcf}")
    
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
        print(f"[*] 已创建输出文件夹: {output_dir}")

    headers = []
    all_svs = []
    all_variants_no_tra_no_bnd = []
    svs_by_type = {}
    
    count_total = 0
    count_filtered = 0

    with open(input_vcf, 'r') as f:
        for line in f:
            # ================= 1. 过滤 Header =================
            if line.startswith('#'):
                if line.startswith('##contig='):
                    match = re.search(r'ID=([^,>]+)', line)
                    if match and match.group(1) not in STANDARD_CHROMS:
                        continue 
                headers.append(line)
                continue
                
            # ================= 2. 过滤变异数据行 =================
            count_total += 1
            cols = line.strip('\n').split('\t')
            
            if len(cols) < 8:
                continue
                
            chrom1 = cols[0]
            filter_status = cols[6]
            info = cols[7]
            alt = cols[4]

            # 提取 SVTYPE
            svtype_match = re.search(r'SVTYPE=([^;]+)', info)
            svtype = svtype_match.group(1).upper() if svtype_match else "UNKNOWN"

            # --- 条件 1: 检查染色体是否为标准染色体 ---
            if chrom1 not in STANDARD_CHROMS:
                continue
                
            is_std_chrom2 = True
            chr2_match = re.search(r'CHR2=([^;]+)', info)
            if chr2_match:
                if chr2_match.group(1) not in STANDARD_CHROMS:
                    is_std_chrom2 = False
            else:
                # 兼容原始的 BND 格式，例如 C[chr12:10372[
                alt_chr_match = re.search(r'[\[\]](.+?):\d+[\[\]]', alt)
                if alt_chr_match and alt_chr_match.group(1) not in STANDARD_CHROMS:
                    is_std_chrom2 = False
                    
            if not is_std_chrom2:
                continue

            # --- 条件 2: FILTER 必须是 PASS ---
            if filter_status.upper() != 'PASS':
                continue

            # --- 条件 3: VAF > 0.05 (从 INFO 列中提取) ---
            vaf_match = re.search(r'\bVAF=([\d\.]+)', info)
            if vaf_match:
                vaf_value = float(vaf_match.group(1))
                if vaf_value <= 0.05:
                    continue  # VAF 太低，丢弃
            else:
                continue # 找不到 VAF，为了严谨也丢弃

            # --- 条件 4: 变异大小 > 50bp (TRA 和 BND 豁免) ---
            if svtype != 'TRA' and svtype != 'BND':
                sv_len = -1
                
                # 优先匹配 SVLEN 或 SVINSLEN (abs 处理负数)
                len_match = re.search(r'\bSVLEN=(-?\d+)', info)
                if len_match:
                    sv_len = abs(int(len_match.group(1)))
                else:
                    inslen_match = re.search(r'\bSVINSLEN=(-?\d+)', info)
                    if inslen_match:
                        sv_len = abs(int(inslen_match.group(1)))
                
                # 如果没写 SVLEN，通过坐标 (END - POS) 自己算物理长度
                if sv_len == -1:
                    end_match = re.search(r'\bEND=(\d+)', info)
                    if end_match:
                        pos1 = int(cols[1])
                        end_pos = int(end_match.group(1))
                        calc_len = abs(end_pos - pos1)
                        if calc_len > 0:
                            sv_len = calc_len

                # 最终裁决：如果长度有效且 <= 50bp，丢弃
                if sv_len != -1 and sv_len <= 50:
                    continue

            # ================= 3. 数据归类 =================
            count_filtered += 1
            all_svs.append(line)
            if svtype != 'TRA' and svtype != 'BND':
                all_variants_no_tra_no_bnd.append(line)
            
            if svtype not in svs_by_type:
                svs_by_type[svtype] = []
            svs_by_type[svtype].append(line)

    # ================= 4. 输出并执行压缩/索引 =================
    print("\n[*] 开始输出、压缩并建立索引...")
    
    # 辅助函数：输出文件并执行 bgzip 和 tabix
    def write_and_compress(filename, lines_to_write):
        out_path = os.path.join(output_dir, filename)
        
        # 写入纯文本 vcf
        with open(out_path, 'w') as f_out:
            f_out.writelines(headers)
            f_out.writelines(lines_to_write)
            
        print(f"✅ 生成文本文件: {filename}")
        
        # 执行 bgzip 和 tabix 命令
        cmd = f"bgzip -c {out_path} > {out_path}.gz && tabix -f -p vcf {out_path}.gz"
        os.system(cmd)
        
        # 压缩完成后，可以考虑删除原始的纯文本 .vcf (可选，这里保留了以供检查)
        # os.remove(out_path)
        
        print(f"   --> 已成功压缩并建索引: {filename}.gz")

    # 处理 ALL 文件
    all_filename = f"all_SV_{len(all_svs)}.vcf"
    write_and_compress(all_filename, all_svs)

    # 处理不含 TRA 和 BND 的文件
    if all_variants_no_tra_no_bnd:
        write_and_compress(f"non_TRA_BND_{len(all_variants_no_tra_no_bnd)}.vcf", all_variants_no_tra_no_bnd)

    # 处理各个亚型文件
    for t_svtype, lines in svs_by_type.items():
        type_filename = f"{t_svtype}_{len(lines)}.vcf"
        write_and_compress(type_filename, lines)

    print(f"\n[+] 过滤完成！")
    print(f"    原始记录数: {count_total}")
    print(f"    保留记录数: {count_filtered} (已剔除杂染色体、非PASS、低VAF及 <=50bp 的变异)")
    print(f"    输出目录: {output_dir}")
    print(f"    (所有文件均已生成对应的 .vcf.gz 和 .tbi 索引)")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="多条件严格过滤 VCF 文件，并按 SVTYPE 拆分输出并压缩建索引。")
    parser.add_argument("-i", "--input", required=True, help="输入的 VCF 文件路径")
    parser.add_argument("-o", "--outdir", required=False, help="输出文件夹路径（可选，如果不填，默认在输入文件同级目录下创建 step04_filtered）")
    
    args = parser.parse_args()

    if args.outdir:
        final_out_dir = os.path.abspath(args.outdir)
    else:
        input_dir = os.path.dirname(os.path.abspath(args.input))
        final_out_dir = os.path.join(input_dir, "step04_filtered")
    
    process_vcf(args.input, final_out_dir)