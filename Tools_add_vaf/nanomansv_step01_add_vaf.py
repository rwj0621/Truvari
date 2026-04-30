import pysam

def add_vaf_nanomonsv(vcf_in_path):
    # 定义输出路径
    out_vcf = vcf_in_path.replace('.vcf', '_vaf.vcf')
    
    # 读入 VCF 并添加 Header 信息
    vcf_in = pysam.VariantFile(vcf_in_path, "r")
    vcf_in.header.info.add("VAF", 1, "Float", "variant_allele_frequency")
    
    # 写入新 VCF
    vcf_out = pysam.VariantFile(out_vcf, 'w', header=vcf_in.header)
    
    for var in vcf_in:
        # 核心逻辑：直接从 TUMOR 样本中读取 VR 和 TR 计算 VAF
        sample = var.samples['TUMOR']
        var.info['VAF'] = sample['VR'] / (sample['TR'] + sample['VR'])
        vcf_out.write(var)
        
    vcf_in.close()
    vcf_out.close()
    print("Done! Result: " + out_vcf)

if __name__ == "__main__":
    # 在这里填入你的文件路径
    add_vaf_nanomonsv("/data/renweijie/Softwares/SV_tools/nanomonsv/HCC1395_PacBio_output/HCC1395_tumor/HCC1395_tumor_PacBio.nanomonsv.result.vcf")