# [Truvari](https://github.com/ACEnglish/truvari)
* 文件路径 A2
  * HCC1395 1788SV
 
          /data/renweijie/data/HCC1395/HCC1395_truvari_process/13059_2022_2816_MOESM4_ESM.xlsx
## 一、把金标准处理成Truvari可用的格式
### 1. xlsx转换成tsv
* 进入Truvari环境

        conda activate Truvari
* 安装 openpyxl 库

        pip install openpyxl
* 运行转换脚本

        python /data/renweijie/data/HCC1395/HCC1395_truvari_process/HCC1395_1788_to_truvari/step01_xlsx_to_tsv.py
### 2. tsv转换成vcf
头部需要去除标准染色体之外的信息
* 运行转换脚本

        python /data/renweijie/data/HCC1395/HCC1395_truvari_process/HCC1395_1788_to_truvari/step02_tsv_to_vcf.py
### 3. 分SV类型提取

        python /data/renweijie/data/HCC1395/HCC1395_truvari_process/HCC1395_1788_to_truvari/step03_split_svtype.py

