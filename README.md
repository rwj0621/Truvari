# [Truvari](https://github.com/ACEnglish/truvari)
## 一、转换脚本文件目录
### 1. 金标准结果转换

        /data/renweijie/data/HCC1395/HCC1395_truvari_process
### 2.nanomansv结果转换

        /data/renweijie/Softwares/SV_tools/nanomonsv/HCC1395_PacBio_output/HCC1395_tumor/filtered_vcf
## 二、把金标准处理成Truvari可用的格式
### 1. [xlsx转换成tsv](https://github.com/rwj0621/Truvari/blob/main/GSD/HCC1395_truvari_process/step01_xlsx_to_tsv.py)
* 进入Truvari环境

        conda activate Truvari
* 安装 openpyxl 库

        pip install openpyxl
* 运行转换脚本

        python /data/renweijie/data/HCC1395/HCC1395_truvari_process/HCC1395_1788_to_truvari/step01_xlsx_to_tsv.py
### 2. [tsv转换成vcf](https://github.com/rwj0621/Truvari/blob/main/GSD/HCC1395_truvari_process/step02_tsv_to_vcf.py)
头部需要去除标准染色体之外的信息
* 运行转换脚本

        python /data/renweijie/data/HCC1395/HCC1395_truvari_process/HCC1395_1788_to_truvari/step02_tsv_to_vcf.py
### 3. [分SV类型提取](https://github.com/rwj0621/Truvari/blob/main/GSD/HCC1395_truvari_process/step03_split_svtype.py)
        python /data/renweijie/data/HCC1395/HCC1395_truvari_process/HCC1395_1788_to_truvari/step03_split_svtype.py
## 三、过滤工具输出结果
* 进入脚本所在目录

        cd /data/renweijie/Softwares/SV_tools/nanomonsv/HCC1395_PacBio_output/HCC1395_tumor/filtered_vcf
### 1.添加vaf
* [nanomansv](https://github.com/rwj0621/Truvari/blob/main/Tools_add_vaf/nanomansv_step01_add_vaf.py)
  
        python step01_add_vaf.py
* sniffles2

        python /data/renweijie/Softwares/SV_tools/sniffles2/filtered_vcf/step01_add_vaf.py \
        -i /data/renweijie/Softwares/SV_tools/sniffles2/HCC1395_somatic.vcf
* SAVANA
SAVANA 里没有专门记录end，因此需要添加end

        python /data/renweijie/Softwares/SV_tools/savana/step01_add_vaf_end_vaf.py
* [severus](https://github.com/rwj0621/Truvari/blob/main/Tools_add_vaf/severus_add_VAF.py)

        python /data/renweijie/Softwares/SV_tools/severus/HCC1395_Somatic_SV_output/preprocess_for_truvari/severus_add_VAF.py
### 2. 将BND改为单条目
* 进入脚本所在目录

        cd /data/renweijie/Softwares/SV_tools/nanomonsv/HCC1395_PacBio_output/HCC1395_tumor/filtered_vcf
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
* SAVANA
 * 自己的结果

   
        python step02_merge_BND.py \
        -i /data/renweijie/Softwares/SV_tools/savana/HCC1395_HIFI/HCC1395.classified.somatic_vaf_end.vcf \
        -o /data/renweijie/Softwares/SV_tools/savana/filtered_vcf/step02_single_BND.vcf
 * 文章结果

         python step02_merge_BND.py \
        -i /data/renweijie/data/HCC1395/HCC1395_severus_variants_calls/savana_H1395.haplotagged.classified.somatic_HIFI_vaf.vcf \
        -o /data/renweijie/Softwares/SV_tools/savana/filtered_vcf/step02_articel_single_BND.vcf
* GRIDSS

         python /data/renweijie/Software/SV_tools/GRIDSS/HCC1395/EA_T_1/step02_merge_BND.py \
        -i /data/renweijie/Software/SV_tools/GRIDSS/HCC1395/EA_T_1/WGS_EA_T_1.gripss.filtered.vcf \
        -o /data/renweijie/Software/SV_tools/GRIDSS/HCC1395/EA_T_1/step02_articel_single_BND.vcf
        
        
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
* SAVANA
  * 自己的结果

       python step03_BND_to_sv.py \
       -i /data/renweijie/Softwares/SV_tools/savana/filtered_vcf/step02_single_BND.vcf \
       -o /data/renweijie/Softwares/SV_tools/savana/filtered_vcf/step03_SVType.vcf
   * 文章结果

        python step03_BND_to_sv.py \
       -i /data/renweijie/Softwares/SV_tools/savana/filtered_vcf/step02_articel_single_BND.vcf \
       -o /data/renweijie/Softwares/SV_tools/savana/filtered_vcf/step03_article_SVType.vcf
* GRIDSS

        python /data/renweijie/Software/SV_tools/GRIDSS/HCC1395/EA_T_1/step03_BND_to_sv.py \
       -i /data/renweijie/Software/SV_tools/GRIDSS/HCC1395/EA_T_1/step02_articel_single_BND.vcf \
       -o /data/renweijie/Software/SV_tools/GRIDSS/HCC1395/EA_T_1/step03_SVType.vcf
      
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
* SAVANA

        python step04_filter.py \
       -i /data/renweijie/Softwares/SV_tools/savana/filtered_vcf/step03_SVType.vcf
* GRIDSS

        python /data/renweijie/Software/SV_tools/GRIDSS/HCC1395/EA_T_1/step04_filter.py \
       -i /data/renweijie/Software/SV_tools/GRIDSS/HCC1395/EA_T_1/step03_SVType.vcf
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
### 2.运行脚本 
只需要修改一下参数
* -b tp-base.vcf.gz
* -c tp-comp.vcf.gz
* -p 输出图的前缀
* -d 根输出目录
* 运行脚本

        python /data/renweijie/Softwares/Truvari/breakpoint_shift/plot_breakpoint_shift.py \
        -b /data/renweijie/Softwares/Truvari/nanomansv_severus_output/HCC1395_ALL/tp-base.vcf.gz \
        -c /data/renweijie/Softwares/Truvari/nanomansv_severus_output/HCC1395_ALL/tp-comp.vcf.gz \
        -p Nano_vs_severus \
        -d /data/renweijie/Softwares/Truvari/breakpoint_shift/nanomansv_severus
## 六、合并结果
### 1. 解压缩需要的文件

        gunzip -c /data/renweijie/Softwares/Truvari/nanomansv_severus_output/HCC1395_ALL/tp-base.vcf.gz > /data/renweijie/Softwares/Truvari/nanomansv_severus_output/HCC1395_ALL/tp-base.vcf
### 2.Survari合并文件路径
* nanomansv_severus 保留severus 

        /data/renweijie/Softwares/Truvari/nanomansv_severus_output/HCC1395_ALL/tp-comp.vcf
* severus_sniffles2 保留severus

        /data/renweijie/Softwares/Truvari/severus_sniffles2_output/HCC1395_ALL/tp-base.vcf
* severus_nanomansv_sniffles2 保留severus

        /data/renweijie/Softwares/Truvari/severus_sniffles2_nanomansv_output/HCC1395_ALL/tp-base.vcf
* severus_nanomansv_sniffles2 保留nanomansv

        /data/renweijie/Softwares/Truvari/nanomansv_sniffles2_severus_output/HCC1395_ALL/tp-base.vcf
* nanomansv_sniffles2 保留nanomansv

        /data/renweijie/Softwares/Truvari/nanomansv_sniffles2_output/HCC1395_ALL/tp-base.vcf
* 经金标准验证过的 severus TP

        /data/renweijie/Softwares/SV_tools/severus/HCC1395_Somatic_SV_output/somatic_SVs/Truvari_output/severus_HCC1395_all/tp-comp.vcf
* 经金标准验证过的 nanomansv TP

        /data/renweijie/Softwares/SV_tools/nanomonsv/HCC1395_PacBio_output/HCC1395_tumor/Truvari_output/nanomansv_HCC1395_all/tp-comp.vcf
* 经金标准验证过的 sniffles2 TP

        /data/renweijie/Softwares/SV_tools/sniffles2/Truvari_output/sniffles2_HCC1395_all/tp-comp.vcf
        
### 3.合并severus相关结果
* 激活survivor环境

         conda activate survivor
* 进入SURVIVOR所在目录

        cd /data/renweijie/Softwares/Survivor/SURVIVOR-master/Debug
* 创建合并vcf列表
* 删除所有 Windows 换行符
  
        sed -i 's/\r//g' /data/renweijie/Softwares/Truvari/HCC1395_merged_SV/step01_severus_merged.txt

* 确保文件末尾有一个换行符
  
        sed -i '$a\' /data/renweijie/Softwares/Truvari/HCC1395_merged_SV/step01_severus_merged.txt
* 运行survivor

        ./SURVIVOR merge /data/renweijie/Softwares/Truvari/HCC1395_merged_SV/step01_severus_merged.txt 10 1 1 0 0 50 /data/renweijie/Softwares/Truvari/HCC1395_merged_SV/step01_severus_consensus.vcf
### 4.找nanomansv和sniffles2共识SV
思路：合并 三者共识SV 与 nanomansv和sniffles2共识SV，仅筛选合并后SUPP=1的结果（应该是38）保留 nanomansv 断点
* 创建合并vcf列表
* 删除所有 Windows 换行符
  
        sed -i 's/\r//g' /data/renweijie/Softwares/Truvari/HCC1395_merged_SV/step02_nanomansv_merged.txt

* 确保文件末尾有一个换行符
  
        sed -i '$a\' /data/renweijie/Softwares/Truvari/HCC1395_merged_SV/step02_nanomansv_merged.txt
* 运行survivor

        ./SURVIVOR merge /data/renweijie/Softwares/Truvari/HCC1395_merged_SV/step02_nanomansv_merged.txt 10 1 1 0 0 50 /data/renweijie/Softwares/Truvari/HCC1395_merged_SV/step02_nanomansv_sniffles2_only_consensus.vcf
* 仅提取 SUPP=1的结果

        bcftools view -i 'INFO/SUPP=="1"' /data/renweijie/Softwares/Truvari/HCC1395_merged_SV/step02_nanomansv_sniffles2_only_consensus.vcf -Ov -o /data/renweijie/Softwares/Truvari/HCC1395_merged_SV/step03_nanomansv_sniffles2_no_severus.vcf
### 5.找仅nanomansv召回的 TP
* 创建合并vcf列表
* 删除所有 Windows 换行符
  
        sed -i 's/\r//g' /data/renweijie/Softwares/Truvari/HCC1395_merged_SV/step04_nanomansv_TP_consensus.txt

* 确保文件末尾有一个换行符
  
        sed -i '$a\' /data/renweijie/Softwares/Truvari/HCC1395_merged_SV/step04_nanomansv_TP_consensus.txt
* 运行survivor

        ./SURVIVOR merge /data/renweijie/Softwares/Truvari/HCC1395_merged_SV/step04_nanomansv_TP_consensus.txt 10 1 1 0 0 50 /data/renweijie/Softwares/Truvari/HCC1395_merged_SV/step04_nanomansv_TP_consencus.vcf
* 仅提取 SUPP=1的结果

        bcftools view -i 'INFO/SUPP=="1"' /data/renweijie/Softwares/Truvari/HCC1395_merged_SV/step04_nanomansv_TP_consencus.vcf -Ov -o /data/renweijie/Softwares/Truvari/HCC1395_merged_SV/step04_nanomansv_TP_only.vcf
### 5.找仅sniffles2召回的 TP
* 创建合并vcf列表
* 删除所有 Windows 换行符
  
        sed -i 's/\r//g' /data/renweijie/Softwares/Truvari/HCC1395_merged_SV/step05_sniffles2_TP_consensus.txt

* 确保文件末尾有一个换行符
  
        sed -i '$a\' /data/renweijie/Softwares/Truvari/HCC1395_merged_SV/step05_sniffles2_TP_consensus.txt
* 运行survivor

        ./SURVIVOR merge /data/renweijie/Softwares/Truvari/HCC1395_merged_SV/step05_sniffles2_TP_consensus.txt 10 1 1 0 0 50 /data/renweijie/Softwares/Truvari/HCC1395_merged_SV/step05_sniffles2_TP_consencus.vcf
* 仅提取 SUPP=1的结果

        bcftools view -i 'INFO/SUPP=="1"' /data/renweijie/Softwares/Truvari/HCC1395_merged_SV/step05_sniffles2_TP_consencus.vcf -Ov -o /data/renweijie/Softwares/Truvari/HCC1395_merged_SV/step05_sniffles2_TP_only.vcf


       



  

  
        

   

