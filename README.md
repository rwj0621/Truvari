# [Truvari](https://github.com/ACEnglish/truvari)
## 一、转换脚本文件目录
### 1. 金标准结果转换

        /data/renweijie/data/HCC1395/HCC1395_truvari_process
### 2.nanomansv结果转换

        /data/renweijie/Softwares/SV_tools/nanomonsv/HCC1395_PacBio_output/HCC1395_tumor/filtered_vcf
### 3.severus 结果转换
### 4.sniffles2 结果转换
## 二、把金标准处理成Truvari可用的格式
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
## 三、过滤工具输出结果
* 进入脚本所在目录

        cd /data/renweijie/Softwares/SV_tools/nanomonsv/HCC1395_PacBio_output/HCC1395_tumor/filtered_vcf
### 1.添加vaf

        python step01_add_vaf.py
### 2. 将BND改为单条目

        python step02_merge_BND.py \
        -i /data/renweijie/Softwares/SV_tools/nanomonsv/HCC1395_PacBio_output/HCC1395_tumor/HCC1395_tumor_PacBio.nanomonsv.result_vaf.vcf \
        -o /data/renweijie/Softwares/SV_tools/nanomonsv/HCC1395_PacBio_output/HCC1395_tumor/filtered_vcf/step02_single_BND.vcf
### 3.将BND转换为具体的SV类型

       python step03_BND_to_sv.py \
       -i /data/renweijie/Softwares/SV_tools/nanomonsv/HCC1395_PacBio_output/HCC1395_tumor/filtered_vcf/step02_single_BND.vcf \
       -o /data/renweijie/Softwares/SV_tools/nanomonsv/HCC1395_PacBio_output/HCC1395_tumor/filtered_vcf/step03_SVType.vcf
### 4.过滤SV 并生成总的和分SV类型的结果
* 只保留标准染色体
* 只保留变异大小 >50bp
* 只保留 filter pass
* 过滤掉 vaf<0.05
* 运行脚本

        python step04_filter.py \
       -i /data/renweijie/Softwares/SV_tools/nanomonsv/HCC1395_PacBio_output/HCC1395_tumor/filtered_vcf/step03_SVType.vcf 
        

   

