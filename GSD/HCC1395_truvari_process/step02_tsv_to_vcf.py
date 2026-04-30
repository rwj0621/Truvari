import pandas as pd
import os

# ================= 1. 路径配置 =================
base_dir = "/data/renweijie/data/HCC1395/HCC1395_truvari_process/HCC1395_1788_to_truvari/output"
input_file = os.path.join(base_dir, "step01_1788.tsv")
output_vcf = os.path.join(base_dir, "step02_1788.vcf")
fai_file = "/data/renweijie/data/GRCh38/GRCh38.d1.vd1.fa.fai"

# ================= 2. 定义标准染色体及排序 =================
# 强制规定生物学顺序：chr1 到 chr22, 然后 chrX, chrY
CHROM_ORDER = [f"chr{i}" for i in range(1, 23)] + ["chrX", "chrY"]

# ================= 3. 读取并完美排序数据 =================
df = pd.read_csv(input_file, sep='\t')

if df.empty:
    raise ValueError(f"输入文件 {input_file} 为空，请检查文件内容！")

df = df.dropna(subset=['Chrom1', 'Pos1', 'Chrom2', 'Pos2', 'SV_type'])

# 去除可能存在的空格
df['Chrom1'] = df['Chrom1'].astype(str).str.strip()
df['Chrom2'] = df['Chrom2'].astype(str).str.strip()

# 1. 将 Chrom1 设置为具有严格大小顺序的分类类型（确保 chr2 排在 chr10 前面）
df['Chrom1'] = pd.Categorical(df['Chrom1'], categories=CHROM_ORDER, ordered=True)
# 2. 确保 Pos1 是数字，方便后续按坐标大小排序
df['Pos1'] = pd.to_numeric(df['Pos1'], errors='coerce')
# 3. 按染色体先后、再按坐标从小到大进行排序
df = df.sort_values(['Chrom1', 'Pos1'])

size_col_name = [col for col in df.columns if 'SV_Size' in col][0]

# ================= 4. 写入 VCF =================
with open(output_vcf, 'w') as f:
    f.write('##fileformat=VCFv4.2\n')
    f.write('##FILTER=<ID=PASS,Description="All filters passed">\n')
    f.write('##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type of structural variant">\n')
    f.write('##INFO=<ID=END,Number=1,Type=Integer,Description="End position of the variant">\n')
    f.write('##INFO=<ID=SVLEN,Number=1,Type=Integer,Description="Length of the variant">\n')
    f.write('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n')
    
    # --- 写入 Contig (字典匹配 + 排序输出) ---
    fai_dict = {}
    if os.path.exists(fai_file):
        with open(fai_file, 'r') as fai:
            for line in fai:
                cols = line.strip().split('\t')
                if len(cols) >= 2:
                    fai_dict[cols[0]] = cols[1]  # 先把 fai 存成字典
                    
    # 严格按照 CHROM_ORDER 的顺序写入 Header，彻底过滤且排序完美
    for chrom in CHROM_ORDER:
        if chrom in fai_dict:
            f.write(f"##contig=<ID={chrom},length={fai_dict[chrom]}>\n")
    
    f.write('#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE\n')
    
    # --- 写入 Body ---
    for _, row in df.iterrows():
        chrom1 = str(row['Chrom1'])
        chrom2 = str(row['Chrom2'])
        
        p1 = int(row['Pos1'])
        p2 = int(float(row['Pos2']))
        svtype = str(row['SV_type']).strip().upper()
        
        if chrom1 == chrom2:
            vcf_pos = min(p1, p2)
            vcf_end = max(p1, p2)
        else:
            vcf_pos = p1
            vcf_end = p1 
            
        sv_size_raw = row[size_col_name]
        sv_size = abs(int(float(sv_size_raw))) if pd.notna(sv_size_raw) else 0
        svlen = -sv_size if svtype == 'DEL' else sv_size
        
        info = f"SVTYPE={svtype};END={vcf_end};SVLEN={svlen}"
        line = f"{chrom1}\t{vcf_pos}\t.\tN\t<{svtype}>\t.\tPASS\t{info}\tGT\t1/1\n"
        f.write(line)

print(f"VCF 纯文本文件已生成 (Header 和数据均已完美排序): {output_vcf}")

# ================= 5. 执行排序、压缩和索引 =================
cmd = f"bcftools sort {output_vcf} -Oz -o {output_vcf}.gz && tabix -f -p vcf {output_vcf}.gz"
print(f"正在执行命令: {cmd}")
os.system(cmd)