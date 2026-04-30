import sys
import os
import re
import subprocess
import gzip
from collections import defaultdict

# ==============================================================================
# 全局参数配置
# ==============================================================================
# 支持 .vcf 或 .vcf.gz 输入
INPUT_VCF = "/data/renweijie/Software/SV_tools/GRIDSS/HCC1395/WGS_IL_3/WGS_IL_T_3.gripss.filtered.vcf.gz"
OUTPUT_DIR = "/data/renweijie/Software/SV_tools/Truvari/process_vcf/HCC1395/GRIDSS/WGS_IL_T_3"
SAMPLE_PREFIX = "WGS_IL_T_3"

# 过滤阈值
TAF_THRESHOLD = 0.05
STANDARD_CHROMS = set([f"chr{i}" for i in range(1, 23)] + ["chrX", "chrY"])
# ==============================================================================

def chrom_key(chrom):
    """
    将染色体名称映射为可排序的整数，用于排序。
    """
    c = chrom.replace("chr", "").upper()
    if c == "X": return 23
    if c == "Y": return 24
    if c in ["M", "MT"]: return 25
    try: return int(c)
    except ValueError: return 999

def extract_taf(info_field):
    """从 INFO 字段解析 TAF 标签的值"""
    match = re.search(r'TAF=([-+]?\d*\.\d+|\d+)', info_field)
    if match:
        return float(match.group(1))
    return 1.0  # 若无 TAF 标签，则不因频率过滤

def get_locus_key(chrom1, pos1, chrom2, pos2):
    """生成唯一配对 Key，确保同一对 BND 拥有相同的 Key"""
    k1, k2 = chrom_key(chrom1), chrom_key(chrom2)
    if k1 < k2 or (k1 == k2 and pos1 <= pos2):
        return f"{chrom1}_{pos1}_{chrom2}_{pos2}"
    else:
        return f"{chrom2}_{pos2}_{chrom1}_{pos1}"

def run_bcftools_pipeline(vcf_path, out_gz_name):
    """调用 bcftools 进行排序、压缩并建立索引"""
    out_gz_path = os.path.join(OUTPUT_DIR, out_gz_name)
    try:
        cmd_sort = f"bcftools sort -O z -o {out_gz_path} {vcf_path}"
        subprocess.run(cmd_sort, shell=True, check=True, stderr=subprocess.PIPE)
        
        cmd_tabix = f"tabix -f -p vcf {out_gz_path}"
        subprocess.run(cmd_tabix, shell=True, check=True, stderr=subprocess.PIPE)
        
        print(f"   已生成: {os.path.basename(vcf_path)} 以及 {out_gz_name} (+.tbi)")
    except subprocess.CalledProcessError as e:
        print(f"  [!] bcftools 报错: {e.stderr.decode()}")

def process_vcf():
    # 1. 自动检查并建立输出目录
    if not os.path.exists(OUTPUT_DIR):
        print(f"[*] 目录 {OUTPUT_DIR} 不存在，正在自动创建...")
        os.makedirs(OUTPUT_DIR, exist_ok=True)

    if not os.path.exists(INPUT_VCF):
        print(f"[!] 错误: 找不到输入文件 {INPUT_VCF}")
        sys.exit(1)

    headers = []
    bnd_groups = defaultdict(list)
    general_records = []
    
    count_total = 0
    count_pass = 0

    print(f"[*] 正在读取输入并应用过滤 (PASS, TAF >= {TAF_THRESHOLD}, 人类染色体): {INPUT_VCF}")
    
    open_func = gzip.open if INPUT_VCF.endswith(".gz") else open
    mode = 'rt' if INPUT_VCF.endswith(".gz") else 'r'

    with open_func(INPUT_VCF, mode) as f:
        for line in f:
            if line.startswith('#'):
                if "##INFO=<ID=SVTYPE" in line:
                    headers.append(line)
                    headers.append('##INFO=<ID=SVLEN,Number=1,Type=Integer,Description="Difference in length between REF and ALT variants">\n')
                    headers.append('##INFO=<ID=END,Number=1,Type=Integer,Description="End position of the variant described by this record">\n')
                else:
                    headers.append(line)
                continue

            count_total += 1
            cols = line.strip('\n').split('\t')
            if len(cols) < 8: continue
            
            chrom1, pos1, vid, alt, filter_col, info = cols[0], int(cols[1]), cols[2], cols[4], cols[6], cols[7]
            
            # --- 过滤逻辑 ---
            # 1. 过滤非 PASS
            if filter_col.upper() != "PASS": continue
            # 2. 过滤非标准人类染色体
            if chrom1 not in STANDARD_CHROMS: continue
            # 3. 过滤 TAF < 阈值
            if extract_taf(info) < TAF_THRESHOLD: continue

            count_pass += 1

            svtype_match = re.search(r'SVTYPE=([^;]+)', info)
            original_svtype = svtype_match.group(1).upper() if svtype_match else "UNKNOWN"

            if original_svtype == "BND" or "SVTYPE=BND" in info:
                match = re.search(r'([\[\]])(.+?):(\d+)([\[\]])', alt)
                is_sgl = False
                if match:
                    _, chrom2, pos2_str, _ = match.groups()
                    pos2 = int(pos2_str)
                    if chrom2 not in STANDARD_CHROMS: continue
                elif "." in alt:
                    chrom2, pos2, is_sgl = chrom1, pos1, True
                else: continue

                # 提取 EVENT 或 MATEID 标识用于配对
                event_m = re.search(r'\bEVENT=([^;]+)', info)
                mate_m = re.search(r'\bMATEID=([^;]+)|\bMATE_ID=([^;]+)', info)
                mid = mate_m.group(1) or mate_m.group(2) if mate_m else ""

                if event_m: pair_key = f"EVENT_{event_m.group(1)}"
                elif vid and mid and not is_sgl: pair_key = f"MATE_{'_'.join(sorted([vid, mid]))}"
                else: pair_key = f"LOCUS_{get_locus_key(chrom1, pos1, chrom2, pos2)}"

                bnd_groups[pair_key].append({
                    'line': line, 'cols': cols, 'chrom1': chrom1, 'pos1': pos1, 
                    'chrom2': chrom2, 'pos2': pos2, 'alt': alt, 'info': info, 'is_sgl': is_sgl
                })
            else:
                general_records.append({'line': line + "\n", 'svtype': original_svtype})

    print(f"[*] 过滤统计: 总数据行 {count_total}, 通过过滤保留 {count_pass}")
    print(f"[*] 正在合并记录并解析类型 (执行染色体核型优先级保留逻辑)...")
    
    final_records_by_type = defaultdict(list)
    all_records = []
    no_tra_records = []

    # 处理常规变异 (非 BND)
    for rec in general_records:
        t = rec['svtype']
        final_records_by_type[t].append(rec['line'])
        all_records.append(rec['line'])
        if t not in ['TRA', 'BND']: no_tra_records.append(rec['line'])

    # 处理 BND 记录并去重
    for pair_key, records in bnd_groups.items():
        records.sort(key=lambda x: (chrom_key(x['chrom1']), x['pos1']))
        primary = records[0]
        c1, p1, c2, p2, alt, info, is_sgl = primary['chrom1'], primary['pos1'], primary['chrom2'], primary['pos2'], primary['alt'], primary['info'], primary['is_sgl']
        
        if not is_sgl:
            match = re.search(r'([\[\]])(.+?):(\d+)([\[\]])', alt)
            b1, _, _, b2 = match.groups()
            is_n = alt[0].upper() in "ATGCN"
            ins_len = len(re.sub(r'[\[\]].+?:\d+[\[\]]', '', alt)) - 1
            del_len = abs(p2 - p1)
            
            if c1 != c2: 
                parsed_t = "TRA"
            else:
                if b1 == '[' and b2 == '[': parsed_t = "DEL" if is_n else "INV"
                elif b1 == ']' and b2 == ']': parsed_t = "INV" if is_n else "DUP"
                else: parsed_t = "BND"
                if parsed_t in ["DEL", "DUP"] and ins_len > del_len: parsed_t = "INS"
        else:
            parsed_t, ins_len, del_len = "INS", len(alt.replace('.', '')) - 1, 0

        # --- 修正标签 ---
        new_info = re.sub(r'SVTYPE=BND', f'SVTYPE={parsed_t}', info)
        if "SVTYPE=" not in new_info: new_info += f";SVTYPE={parsed_t}"
        
        # 核心逻辑：解决 END < POS 索引问题
        # 如果是跨染色体变异 (TRA) 或单端变异，END 必须等于 POS，否则 tabix 索引会报错
        if c1 != c2 or parsed_t == "TRA":
            end_val = p1
        else:
            end_val = max(p1, p2)

        new_info = re.sub(r'END=\d+', f'END={end_val}', new_info)
        if "END=" not in new_info: new_info += f";END={end_val}"
        
        if parsed_t != "TRA":
            svlen = ins_len if parsed_t == "INS" else del_len
            if "SVLEN=" not in new_info: new_info += f";SVLEN={svlen}"

        primary['cols'][7] = new_info
        final_line = "\t".join(primary['cols']).strip() + "\n"
        
        final_records_by_type[parsed_t].append(final_line)
        all_records.append(final_line)
        if parsed_t not in ['TRA', 'BND']: no_tra_records.append(final_line)

    # --- 输出逻辑 ---
    print(f"[*] 正在分类型输出结果文件...")
    stats_file = os.path.join(OUTPUT_DIR, f"{SAMPLE_PREFIX}.tsv")
    
    with open(stats_file, 'w') as f_stat:
        f_stat.write("SVTYPE\tCount\n")
        # 1. 各个子类型文件
        for t in sorted(final_records_by_type.keys()):
            lines = final_records_by_type[t]
            f_stat.write(f"{t}\t{len(lines)}\n")
            vname = f"{SAMPLE_PREFIX}_{t}_{len(lines)}.vcf"
            vpath = os.path.join(OUTPUT_DIR, vname)
            with open(vpath, 'w') as out:
                out.writelines(headers); out.writelines(lines)
            run_bcftools_pipeline(vpath, vname + ".gz")

        # 2. 汇总文件 (ALL 和 ALL_no_TRA)
        for name, recs in [("ALL", all_records), ("ALL_no_TRA", no_tra_records)]:
            f_stat.write(f"{name}\t{len(recs)}\n")
            vname = f"{SAMPLE_PREFIX}_{name}_{len(recs)}.vcf"
            vpath = os.path.join(OUTPUT_DIR, vname)
            with open(vpath, 'w') as out:
                out.writelines(headers); out.writelines(recs)
            run_bcftools_pipeline(vpath, vname + ".gz")

    print(f"[*] 所有任务完成！统计表已保存至: {stats_file}")

if __name__ == "__main__":
    process_vcf()