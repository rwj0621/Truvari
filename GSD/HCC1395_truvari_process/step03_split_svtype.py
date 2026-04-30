import os
import re

# ================= 1. 路径配置 =================
base_dir = "/data/renweijie/data/HCC1395/HCC1395_truvari_process/HCC1395_1788_to_truvari/output"
input_vcf = os.path.join(base_dir, "step02_1788_TRA.vcf")

if not os.path.exists(input_vcf):
    raise FileNotFoundError(f"找不到输入文件: {input_vcf}，请确认上一步是否生成了纯文本 VCF。")

# ================= 2. 读取并分类数据 =================
headers = []
variants_by_type = {}

all_variants_with_bnd = []
all_variants_no_bnd = []
all_variants_no_bnd_no_tra = []

print(f"[*] 正在读取并解析: {input_vcf}")
with open(input_vcf, 'r') as f:
    for line in f:
        # 如果是表头行（以 # 开头），单独存起来，每个拆分后的文件都需要
        if line.startswith('#'):
            headers.append(line)
        else:
            cols = line.strip().split('\t')
            if len(cols) < 8:
                continue
                
            info_col = cols[7]
            
            # 使用正则表达式从 INFO 列中精准提取 SVTYPE 的值 (例如 SVTYPE=DEL)
            match = re.search(r'SVTYPE=([^;]+)', info_col)
            if match:
                svtype = match.group(1).upper()
            else:
                svtype = "UNKNOWN"
            
            # --- 分类归档 ---
            
            # 1. 按具体 SV 类型归档
            if svtype not in variants_by_type:
                variants_by_type[svtype] = []
            variants_by_type[svtype].append(line)
            
            # 2. 存入 "包含 BND" 的汇总列表
            all_variants_with_bnd.append(line)
            
            # 3. 存入 "不含 BND" 的汇总列表
            if svtype != 'BND':
                all_variants_no_bnd.append(line)
            # 4. 在不是BND 的前提下，进一步排除TRA

            if svtype != 'BND':
                if svtype != 'TRA':
                    all_variants_no_bnd_no_tra.append(line)


# ================= 3. 定义输出与压缩的辅助函数 =================
def write_and_compress(filename, lines_to_write):
    """把变异行写入 VCF 文件，并自动调用 bcftools 排序、压缩和建索引"""
    output_path = os.path.join(base_dir, filename)
    count = len(lines_to_write)
    
    with open(output_path, 'w') as out_f:
        out_f.writelines(headers)
        out_f.writelines(lines_to_write)
        
    print(f" 成功生成: {filename} (包含 {count} 个变异)")
    
    # 使用 bcftools sort 自动修复乱序并压缩建索引
    cmd = f"bcftools sort {output_path} -Oz -o {output_path}.gz && tabix -f -p vcf {output_path}.gz"
    os.system(cmd)

# ================= 4. 输出所有拆分和汇总文件 =================
print("\n[*] 开始生成输出文件并执行排序压缩...")

# 1. 循环输出各个亚型的独立文件 (step03_DEL.xxx.vcf 等)
for svtype, lines in variants_by_type.items():
    sv_filename = f"step03_{svtype}.{len(lines)}.vcf"
    write_and_compress(sv_filename, lines)

# 2. 输出包含 BND 的完整版
all_with_bnd_filename = f"step03_ALL_withBND.{len(all_variants_with_bnd)}.vcf"
write_and_compress(all_with_bnd_filename, all_variants_with_bnd)

# 3. 输出不包含 BND 的干净版
all_no_bnd_filename = f"step03_ALL_noBND.{len(all_variants_no_bnd)}.vcf"
write_and_compress(all_no_bnd_filename, all_variants_no_bnd)

# 4.不包含 BND 且不包含 TRA
write_and_compress(f"step03_ALL_noBND_noTRA.{len(all_variants_no_bnd_no_tra)}.vcf", all_variants_no_bnd_no_tra)

print("\n[+] 拆分及汇总任务全部完成！")