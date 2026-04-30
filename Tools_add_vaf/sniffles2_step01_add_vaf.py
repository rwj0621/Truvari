import pysam
import argparse
import os

def add_vaf_sniffles(input_vcf):
    if not os.path.exists(input_vcf):
        return

    out_vcf = input_vcf.replace('.vcf', '_vaf.vcf')

    vcf_in = pysam.VariantFile(input_vcf, "r")

    if "VAF" not in vcf_in.header.info:
        vcf_in.header.info.add("VAF", 1, "Float", "Variant Allele Frequency calculated from DV and DR")

    vcf_out = pysam.VariantFile(out_vcf, 'w', header=vcf_in.header)

    sample_names = list(vcf_in.header.samples)
    if not sample_names:
        return
        
    # ================= 修复的核心逻辑 =================
    # 优先找名字里带 tumor/output 的，或者直接拿第二个样本
    if len(sample_names) >= 2:
        # 如果有两个样本，Tumor 通常是第二个 (Index 1)
        tumor_sample = sample_names[1] 
    else:
        tumor_sample = sample_names[0]
    
    print(f"[*] 实际用于计算 VAF 的样本列名是: {tumor_sample}")
    # =================================================

    count = 0
    for var in vcf_in:
        count += 1
        try:
            dv_raw = var.samples[tumor_sample].get('DV')
            dr_raw = var.samples[tumor_sample].get('DR')

            dv = dv_raw[0] if isinstance(dv_raw, tuple) else dv_raw
            dr = dr_raw[0] if isinstance(dr_raw, tuple) else dr_raw

            if dv is None: dv = 0
            if dr is None: dr = 0

            total_reads = dv + dr
            if total_reads > 0:
                vaf = dv / total_reads
            else:
                vaf = 0.0

            var.info['VAF'] = round(vaf, 4)

        except Exception:
            pass
            
        vcf_out.write(var)

    vcf_in.close()
    vcf_out.close()

    print(f"[+] 处理完成！共计算了 {count} 条变异记录的 VAF。")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="为 Sniffles2 的 VCF 计算 VAF 并写入 INFO 列。")
    parser.add_argument("-i", "--input", required=True, help="输入的 Sniffles2 VCF 文件路径")
    
    args = parser.parse_args()
    add_vaf_sniffles(args.input)