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
* nanomansv
  
        python step01_add_vaf.py
* sniffles2

        python /data/renweijie/Softwares/SV_tools/sniffles2/filtered_vcf/step01_add_vaf.py \
        -i /data/renweijie/Softwares/SV_tools/sniffles2/HCC1395_somatic.vcf
### 2. 将BND改为单条目
* nanomansv

       python step02_merge_BND.py \
        -i /data/renweijie/Softwares/SV_tools/nanomonsv/HCC1395_PacBio_output/HCC1395_tumor/HCC1395_tumor_PacBio.nanomonsv.result_vaf.vcf \
        -o /data/renweijie/Softwares/SV_tools/nanomonsv/HCC1395_PacBio_output/HCC1395_tumor/filtered_vcf/step02_single_BND.vcf
* severus

        python step02_merge_BND.py \
        -i /data/renweijie/Softwares/SV_tools/severus/HCC1395_Somatic_SV_output/somatic_SVs/severus_somatic_vaf.vcf \
        -o /data/renweijie/Softwares/SV_tools/severus/HCC1395_Somatic_SV_output/somatic_SVs/filtered_vcf/step02_single_BND.vcf
* sniffles2

        python step02_merge_BND.py \
        -i /data/renweijie/Softwares/SV_tools/sniffles2/HCC1395_somatic_vaf.vcf \
        -o /data/renweijie/Softwares/SV_tools/sniffles2/filtered_vcf/step02_single_BND.vcf
### 3.将BND转换为具体的SV类型
* nanomansv

       python step03_BND_to_sv.py \
       -i /data/renweijie/Softwares/SV_tools/nanomonsv/HCC1395_PacBio_output/HCC1395_tumor/filtered_vcf/step02_single_BND.vcf \
       -o /data/renweijie/Softwares/SV_tools/nanomonsv/HCC1395_PacBio_output/HCC1395_tumor/filtered_vcf/step03_SVType.vcf
* severus

       python step03_BND_to_sv.py \
       -i /data/renweijie/Softwares/SV_tools/severus/HCC1395_Somatic_SV_output/somatic_SVs/filtered_vcf/step02_single_BND.vcf \
       -o /data/renweijie/Softwares/SV_tools/severus/HCC1395_Somatic_SV_output/somatic_SVs/filtered_vcf/step03_SVType.vcf
* sniffles2

       python step03_BND_to_sv.py \
       -i /data/renweijie/Softwares/SV_tools/sniffles2/filtered_vcf/step02_single_BND.vcf \
       -o /data/renweijie/Softwares/SV_tools/sniffles2/filtered_vcf/step03_SVType.vcf
### 4.过滤SV 并生成总的和分SV类型的结果
* 只保留标准染色体
* 只保留变异大小 >50bp
* 只保留 filter pass
* 过滤掉 vaf<0.05
* nanomansv

        python step04_filter.py \
       -i /data/renweijie/Softwares/SV_tools/severus/HCC1395_Somatic_SV_output/somatic_SVs/filtered_vcf/step03_SVType.vcf
* severus

        python step04_filter.py \
       -i /data/renweijie/Softwares/SV_tools/severus/HCC1395_Somatic_SV_output/somatic_SVs/filtered_vcf/step03_SVType.vcf
* sniffles2

        python step04_filter.py \
       -i /data/renweijie/Softwares/SV_tools/sniffles2/filtered_vcf/step03_SVType.vcf
## 四、Truvari两两比较
只需要修改一下几个参数
* -b 金标准vcf
* -c 工具vcf
* -o 输出文件路径

        truvari bench \
  -b /data/renweijie/Softwares/SV_tools/severus/HCC1395_Somatic_SV_output/somatic_SVs/filtered_vcf/step04_filtered/TRA_134.vcf.gz \
  -c /data/renweijie/Softwares/SV_tools/sniffles2/filtered_vcf/step04_filtered/TRA_131.vcf.gz \
  -o /data/renweijie/Softwares/Truvari/severus_sniffles2_output/HCC1395_TRA \
  --refdist 500 \
  --typeignore \
  --pctsize 0.7 \
  -p 0 \
  --dup-to-ins \
  -s 0 \
  -S 0 \
  --sizemax 100000000 \
  --passonly
## 五、绘制断点偏移图
根据Truvari结果绘制断点偏移图 分SV大小
### 1.配置绘图环境

        pip install matplotlib seaborn pandas
### 1. DEL
只需要修改一下参数
* -b tp-base.vcf.gz
* -c tp-comp.vcf.gz
* -p 输出图的前缀
* 运行脚本

        python /data/renweijie/Softwares/Truvari/breakpoint_shift/plot_breakpoint_shift_DEL.py \
        -b /data/renweijie/Softwares/Truvari/nanomansv_severus_output/HCC1395_ALL/tp-base.vcf.gz \
        -c /data/renweijie/Softwares/Truvari/nanomansv_severus_output/HCC1395_ALL/tp-comp.vcf.gz \
        -d /data/renweijie/Softwares/Truvari/breakpoint_shift/nanomansv_severus \
        -p Nano_vs_severus



  

  
        

   

