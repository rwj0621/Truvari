#!/usr/bin/env python3
import sys
import os
import re
import glob

# =========================================================================
#  核心配置区
# =========================================================================

# 1. 批量输入的基础目录 (脚本会自动在这个目录下寻找所有的 somaticSV_with_VAF.vcf)
INPUT_BASE_DIR = "/data/renweijie/Software/SV_tools/Manta/HCC1395_Batch_Run"

# 2. 批量输出的基础目录 (脚本会在这里以样本名为名创建子文件夹)
OUTPUT_BASE_DIR = "/data/renweijie/Software/SV_tools/Truvari/process_vcf/HCC1395/Manta"

# 3. VAF 过滤阈值
MIN_VAF = 0.05

# 4. 是否保留中间步骤的临时 VCF 文件 (True 或 False)
KEEP_TMP_FILES = False

STANDARD_CHROMS = set([f"chr{i}" for i in range(1, 23)] + ["chrX", "chrY"])

def step1_dedup_bnd(input_vcf, output_vcf):
    seen_mates = set()
    count_removed = 0
    with open(input_vcf, 'r') as f, open(output_vcf, 'w') as out:
        for line in f:
            if line.startswith('#'):
                out.write(line)
                continue
            cols = line.split('\t')
            if len(cols) < 8: continue

            current_id = cols[2]
            info = cols[7]
            mate_id = ""
            for field in info.split(';'):
                if field.startswith('MATEID=') or field.startswith('MATE_ID='):
                    mate_id = field.split('=')[1]
                    break 
            
            if current_id in seen_mates:
                count_removed += 1
                continue
            else:
                if mate_id: seen_mates.add(mate_id)
                out.write(line)
    return count_removed

def parse_bnd_rigorous(chrom1, pos1, alt):
    match = re.search(r'([\[\]])(.+?):(\d+)([\[\]])', alt)
    if not match: return None, None, None
    
    bracket1, chrom2, pos2, bracket2 = match.groups()
    pos2 = int(pos2)
    if chrom1 != chrom2: return "TRA", chrom2, pos2

    is_n_first = alt[0].upper() in "ATGCN"
    if bracket1 == '[' and bracket2 == '[':
        if is_n_first: return "DEL", chrom2, pos2  
        else: return "INV", chrom2, pos2  
    elif bracket1 == ']' and bracket2 == ']':
        if is_n_first: return "INV", chrom2, pos2 
        else: return "DUP", chrom2, pos2 
    return "BND", chrom2, pos2

def step2_convert_bnd(input_vcf, output_vcf):
    count = 0
    with open(input_vcf, 'r') as f, open(output_vcf, 'w') as out:
        for line in f:
            if line.startswith('#'):
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

            if "SVTYPE=BND" in info:
                svtype, chrom2, pos2 = parse_bnd_rigorous(chrom1, pos1, alt)
                if svtype:
                    new_info = info.replace("SVTYPE=BND", f"SVTYPE={svtype}")
                    if svtype != "TRA" and svtype != "BND":
                        svlen = abs(pos2 - pos1)
                        end_pos = max(pos1, pos2)
                        if "SVLEN=" not in new_info: new_info += f";SVLEN={svlen}"
                        if "END=" not in new_info: new_info += f";END={end_pos}"
                    
                    cols[7] = new_info
                    out.write("\t".join(cols) + "\n" if not cols[-1].endswith("\n") else "\t".join(cols))
                    count += 1
                else:
                    out.write(line)
            else:
                out.write(line)
    return count

def step3_filter_and_split(input_vcf, output_dir):
    headers, all_svs, all_variants_no_tra_no_bnd = [], [], []
    svs_by_type = {}
    
    with open(input_vcf, 'r') as f:
        for line in f:
            if line.startswith('#'):
                if line.startswith('##contig='):
                    match = re.search(r'ID=([^,>]+)', line)
                    if match and match.group(1) not in STANDARD_CHROMS: continue 
                headers.append(line)
                continue
                
            cols = line.strip('\n').split('\t')
            if len(cols) < 8: continue
                
            chrom1 = cols[0]
            filter_status = cols[6]
            info = cols[7]
            alt = cols[4]

            svtype_match = re.search(r'SVTYPE=([^;]+)', info)
            svtype = svtype_match.group(1).upper() if svtype_match else "UNKNOWN"

            if chrom1 not in STANDARD_CHROMS: continue
                
            is_std_chrom2 = True
            chr2_match = re.search(r'CHR2=([^;]+)', info)
            if chr2_match:
                if chr2_match.group(1) not in STANDARD_CHROMS: is_std_chrom2 = False
            else:
                alt_chr_match = re.search(r'[\[\]](.+?):\d+[\[\]]', alt)
                if alt_chr_match and alt_chr_match.group(1) not in STANDARD_CHROMS: is_std_chrom2 = False
            if not is_std_chrom2: continue

            if filter_status.upper() != 'PASS': continue

            vaf_match = re.search(r'\bVAF=([\d\.]+)', info)
            if vaf_match:
                if float(vaf_match.group(1)) <= MIN_VAF: continue

            if svtype != 'TRA' and svtype != 'BND':
                sv_len = -1
                len_match = re.search(r'\bSVLEN=(-?\d+)', info)
                if len_match: sv_len = abs(int(len_match.group(1)))
                else:
                    inslen_match = re.search(r'\bSVINSLEN=(-?\d+)', info)
                    if inslen_match: sv_len = abs(int(inslen_match.group(1)))
                
                if sv_len == -1:
                    end_match = re.search(r'\bEND=(\d+)', info)
                    if end_match:
                        calc_len = abs(int(end_match.group(1)) - int(cols[1]))
                        if calc_len > 0: sv_len = calc_len

                if sv_len != -1 and sv_len <= 50: continue

            all_svs.append(line)
            if svtype != 'TRA' and svtype != 'BND': all_variants_no_tra_no_bnd.append(line)
            if svtype not in svs_by_type: svs_by_type[svtype] = []
            svs_by_type[svtype].append(line)

    def write_and_compress(filename, lines_to_write):
        out_path = os.path.join(output_dir, filename)
        with open(out_path, 'w') as f_out:
            f_out.writelines(headers)
            f_out.writelines(lines_to_write)
        os.system(f"bgzip -f -c {out_path} > {out_path}.gz && tabix -f -p vcf {out_path}.gz")

    write_and_compress(f"all_SV_{len(all_svs)}.vcf", all_svs)
    if all_variants_no_tra_no_bnd:
        write_and_compress(f"non_TRA_BND_{len(all_variants_no_tra_no_bnd)}.vcf", all_variants_no_tra_no_bnd)
    for t_svtype, lines in svs_by_type.items():
        write_and_compress(f"{t_svtype}_{len(lines)}.vcf", lines)
        
    return svs_by_type

def step4_generate_stats(svs_by_type, output_dir, sample_name):
    stat_file = os.path.join(output_dir, f"{sample_name}.tsv")
    core_types = ['DEL', 'DUP', 'INS', 'INV', 'TRA']
    counts = {t: 0 for t in core_types}
    total_all, total_no_tra = 0, 0
    
    for svtype, lines in svs_by_type.items():
        count = len(lines)
        counts[svtype] = count 
        total_all += count
        if svtype != 'TRA': total_no_tra += count
            
    with open(stat_file, 'w') as f:
        # 去掉了 Tool 列
        f.write("Sample\tSVTYPE\tCount\n")
        for t in core_types: f.write(f"{sample_name}\t{t}\t{counts.get(t, 0)}\n")
        for t, c in counts.items():
            if t not in core_types: f.write(f"{sample_name}\t{t}\t{c}\n")
        f.write(f"{sample_name}\tALL\t{total_all}\n")
        f.write(f"{sample_name}\tALL_no_TRA\t{total_no_tra}\n")

# =========================================================================
# 批处理主流程
# =========================================================================
def main():
    # 查找所有目标 VCF 文件
    search_pattern = os.path.join(INPUT_BASE_DIR, "*", "results", "variants", "somaticSV_with_VAF.vcf")
    vcf_files = glob.glob(search_pattern)

    if not vcf_files:
        print(f" 错误: 未在 {INPUT_BASE_DIR} 找到任何 somaticSV_with_VAF.vcf 文件！")
        sys.exit(1)

    print(f"============================================================")
    print(f"[*] 扫描到 {len(vcf_files)} 个 Manta 样本，开始批量过滤并输出至:")
    print(f"[*] {OUTPUT_BASE_DIR}")
    print(f"============================================================\n")

    success_count = 0

    for input_vcf in vcf_files:
        try:
            # 解析样本名 (从路径中向上回溯 3 层拿到样本名，例如 WGS_EA_T_1)
            sample_dir = os.path.dirname(os.path.dirname(os.path.dirname(input_vcf)))
            sample_name = os.path.basename(sample_dir)

            # 拼接并创建该样本对应的专属输出目录
            sample_out_dir = os.path.join(OUTPUT_BASE_DIR, sample_name)
            os.makedirs(sample_out_dir, exist_ok=True)

            print(f"🚀 处理样本: {sample_name} -> {sample_out_dir}")

            tmp_step1 = os.path.join(sample_out_dir, "tmp_step1_dedup.vcf")
            tmp_step2 = os.path.join(sample_out_dir, "tmp_step2_converted.vcf")

            # 1. 去重 BND
            rm_count = step1_dedup_bnd(input_vcf, tmp_step1)
            # 2. 转换 BND
            cv_count = step2_convert_bnd(tmp_step1, tmp_step2)
            # 3. 终极过滤并拆分
            filtered_svs_dict = step3_filter_and_split(tmp_step2, sample_out_dir)
            # 4. 生成统计表 (无Tool版)
            step4_generate_stats(filtered_svs_dict, sample_out_dir, sample_name)
            
            # 清理临时文件
            if not KEEP_TMP_FILES:
                for f in [tmp_step1, tmp_step2]:
                    if os.path.exists(f): os.remove(f)

            print(f"   完成！(去重: {rm_count}, 转换: {cv_count})")
            success_count += 1

        except Exception as e:
            print(f"   样本处理失败 ({input_vcf}): {e}")

    print(f"\n============================================================")
    print(f"🎉 批量处理完毕！成功: {success_count} / {len(vcf_files)}")
    print(f"============================================================")

if __name__ == "__main__":
    main()