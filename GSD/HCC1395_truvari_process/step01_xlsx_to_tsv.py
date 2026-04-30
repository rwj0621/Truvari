import pandas as pd

# 定义输入和输出文件路径 (注意输出后缀改为了 .tsv)
input_excel = '/data/renweijie/data/HCC1395/HCC1395_truvari_process/13059_2022_2816_MOESM4_ESM.xlsx'
output_tsv = '/data/renweijie/data/HCC1395/HCC1395_truvari_process/HCC1395_1788_to_truvari/output/step01_1788.tsv'

# 读取 Excel 文件
df = pd.read_excel(input_excel)

# 保存为 TSV 格式：添加 sep='\t' 参数指定制表符分隔，index=False 表示不保存行索引
df.to_csv(output_tsv, sep='\t', index=False)

print(f"文件已成功转换为: {output_tsv}")
