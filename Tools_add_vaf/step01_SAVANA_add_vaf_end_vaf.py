import pysam
import numpy as np
import os

# ==========================================
# 1. 输入文件路径
# ==========================================
INPUT_VCF = "/data/renweijie/Softwares/SV_tools/savana/HCC1395_HIFI/HCC1395.classified.somatic.vcf"

def add_vaf_and_end_savana(vcf_in_path):
    # 定义输出文件名
    out_vcf = vcf_in_path.replace('.vcf', '_vaf_end.vcf')
    
    # 1. 打开 VCF 文件
    vcf_in = pysam.VariantFile(vcf_in_path, "r")
    
    # 2. 修改 Header (声明字段)
    if "VAF" not in vcf_in.header.info:
        vcf_in.header.info.add("VAF", 1, "Float", "Variant Allele Frequency")
    if "END" not in vcf_in.header.info:
        vcf_in.header.info.add("END", 1, "Integer", "Stop position of the interval")
    
    # 3. 创建输出文件，继承修改后的 Header
    vcf_out = pysam.VariantFile(out_vcf, 'w', header=vcf_in.header)

    # 4. 遍历处理每一行记录
    for var in vcf_in:
        # --- A. 写入 VAF 数据 ---
        # 提取 TUMOUR_AF (Savana 输出是元组，如 (0.5, 0.6)) 并取均值
        af_values = var.info.get('TUMOUR_AF')
        if af_values:
            var.info['VAF'] = float(np.mean(af_values))
        else:
            var.info['VAF'] = 0.0
        
        # 将修改后的记录写入新文件
        vcf_out.write(var)

    # 5. 关闭文件
    vcf_in.close()
    vcf_out.close()
    print(f"Done! 结果已保存至: {out_vcf}")

if __name__ == "__main__":
    if os.path.exists(INPUT_VCF):
        add_vaf_and_end_savana(INPUT_VCF)
    else:
        print(f"错误: 找不到文件 {INPUT_VCF}")
