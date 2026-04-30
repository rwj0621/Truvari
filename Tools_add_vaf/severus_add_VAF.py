import pysam
import os

def add_vaf_severus(vcf_in_path):
    out_vcf = vcf_in_path.replace('.vcf', '_vaf.vcf')
    
    # 打开输入文件
    vcf_reader = pysam.VariantFile(vcf_in_path, "r")
    
    # 在 Header 中添加 VAF 定义
    vcf_reader.header.info.add("VAF", 1, "Float", "variant_allele_frequency")
    
    # 打开输出文件
    vcf_out = pysam.VariantFile(out_vcf, 'w', header=vcf_reader.header)
    
    # 逐行处理
    for var in vcf_reader:
        # 1. 获取 sample_id
        sample_id = list(var.samples.keys())[0]
        # 2. 将 FORMAT 中的 VAF 赋值给 INFO 中的 VAF
        var.info['VAF'] = var.samples[sample_id]['VAF']
        # 3. 写入新文件
        vcf_out.write(var)
    
    # 关闭文件
    vcf_out.close()
    vcf_reader.close()
    print(f"成功！处理后的文件保存在: {out_vcf}")

# ==========================================
# 输入文件路径
# ==========================================
if __name__ == "__main__":
    # 修改下面这一行的路径
    input_file = "/data/renweijie/Softwares/SV_tools/severus/HCC1395_Somatic_SV_output/somatic_SVs/severus_somatic.vcf"
    
    if os.path.exists(input_file):
        add_vaf_severus(input_file)
    else:
        print(f"错误：找不到文件 {input_file}")