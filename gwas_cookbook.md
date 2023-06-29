# **GWAS**

##  #配置流程软件：

```shell
#tree：
下载：http://mama.indstate.edu/users/ice/tree/
解压：tar -zxvf tree-1.8.0.tgz
进入安装目录： cd tree-1.8.0
安装：make install
测试：tree

#sratoolkit：
下载：https://trace.ncbi.nlm.nih.gov/Traces/sra/sra.cgi?view=software
解压：tar -zxvf sratoolkit.2.10.5-centos_linux64.tar.gz
修改环境变量：echo 'PATH=$PATH:~/bigdata/gaojy/biosoft/sratoolkit/bin' >> ~/.bashrc
修改生效：source ~/.bashrc
下载：nohup prefetch SRR824846 &  
格式转换：fastq-dump --gzip --split-files SRR824846.sra
公共数据可通过sratoolkit下载

#seqtk
下载：https://github.com/lh3/seqtk/releases/tag/v1.3
切目录：cd seqtk
安装：make

#fastQC:：py2.7
conda install fastqc

#multiQC：py3.6
conda install multiqc

#fastp：
下载；wget http://opengene.org/fastp/fastp.0.23.1
改名：mv fastp.0.23.1 fastp
添加执行权限：chmod a+x ./fastp
修改环境变量：echo 'PATH=$PATH:/bigdata/gaojy/biosoft/fastp/fastp/' >> ~/.bashrc
修改生效：source ~/.bashrc

#cutadapt：
conda create -n cutadaptenv cutadapt
conda activate cutadaptenv
cutadapt --version

#bwa：
下载：https://sourceforge.net/projects/bio-bwa/
解压: tar -jxvf bwa-0.7.15.tar.bz2
安装：cd bwa
安装：make
修改环境变量：echo 'PATH=$PATH:/bigdata/gaojy/biosoft/bwa/' >> ~/.bashrc
修改生效：source ~/.bashrc
建索引：bwa index -a bwtsw hg38.fna

#samtools：
下载：https://github.com/samtools/samtools/releases
解压：tar jxvf samtools.tar.bz2
安装：./configure --prefix=/bigdata/gaojy/biosoft/samtools
安装：make
安装：make install
测试：./samtools --help
修改环境变量：echo 'PATH=$PATH:/bigdata/gaojy/biosoft/samtools/' >> ~/.bashrc
修改生效：source ~/.bashrc

#jdk：
下载：http://www.oracle.com/technetwork/java/javase/downloads/jdk8-downloads-2133151.html
解压：tar zxvf jdk-8u131-linux-x64.tar.gz
修改环境变量：echo 'PATH=$PATH:/bigdata/gaojy/biosoft/jdk/bin/' >> ~/.bashrc
修改生效：source ~/.bashrc
测试：java -version

#gatk：
下载：https://github.com/broadinstitute/gatk/releases/tag/4.1.7.0
解压：unzip gatk.zip
修改环境变量：echo 'PATH=$PATH:/bigdata/gaojy/biosoft/gatk' >> ~/.bashrc
修改生效：source ~/.bashrc
测试：./gatk --help

#plink：
下载：http://www.cog-genomics.org/plink2/
解压：unzip plink.zip
修改环境变量：echo 'PATH=$PATH:/bigdata/gaojy/biosoft/plink' >> ~/.bashrc
修改生效：source ~/.bashrc

#tabix：
下载：wget https://sourceforge.net/projects/samtools/files/tabix/tabix-0.2.6.tar.bz2
解压：tar xjvf tabix-0.2.6.tar.bz2
切换目录：cd tabix
编译安装：make
修改环境变量：echo 'PATH=$PATH:/bigdata/gaojy/biosoft/tabix' >> ~/.bashrc
修改生效：source ~/.bashrc

#bcftools：
下载：wget https://github.com/samtools/bcftools/releases/download/1.8/bcftools-1.8.tar.bz2
解压：tar xjvf bcftools-1.8.tar.bz2
切换目录：cd bcftools
编译：./configure
编译：make
修改环境变量：echo 'PATH=$PATH:/bigdata/gaojy/biosoft/bcftools/' >> ~/.bashrc
修改生效：source ~/.bashrc

#bedtools：
下载：https://github.com/arq5x/bedtools2/archive/v2.25.0.tar.gz
解压：tar zxvf v2.25.0
进入目录：cd bedtools2-2.25.0
安装：make
修改环境变量：echo 'PATH=$PATH:/bigdata/gaojy/biosoft/bedtools/' >> ~/.bashrc
修改生效：source ~/.bashrc

#crossmap：py3.6
conda install -c bioconda pyBigWig
pip install CrossMap
which CrossMap.py#锁定位置

#beagle:
下载：https://faculty.washington.edu/browning/beagle/b5_2.html#download
千人组:https://bochet.gcc.biostat.washington.edu/beagle/1000_Genomes_phase3_v5a/
#直接可用

#shapeit:
下载：http://mathgen.stats.ox.ac.uk/genetics_software/shapeit/shapeit.html
解压：tar -zxvf shapeit.v2.r900.glibcv2.17.linux.tar.gz
千人组：https://mathgen.stats.ox.ac.uk/impute/data_download_1000G_phase1_integrated_SHAPEIT2_16-06-14.html
#直接使用

#impute2:
下载：wget https://mathgen.stats.ox.ac.uk/impute/impute_v2.3.2_x86_64_static.tgz
解压：tar xzvf impute_v2.3.2_x86_64_static.tgz
千人组:https://mathgen.stats.ox.ac.uk/impute/impute_v2.html#download
#直接可用

#glimpse_v2:
下载：https://github.com/odelaneau/glimpse/releases
解压：tar xzvf glimpse_v2.0.0-27-g0919952_20221207.tar.gz
#直接可用

#quilt:
conda activate




#gemma：
下载：https://github.com/genetics-statistics/GEMMA/releases
#提供binary file,解压直接用

#smr
下载：https://yanglab.westlake.edu.cn/software/smr/#Download
GTEx：https://yanglab.westlake.edu.cn/data/SMR/GTEx_V8_cis_eqtl_summary.html
#直接可用


```

##   #下载文件：

```shell
人类基因组：
人类基因组文件版本较多（NCBI、Ensembl、UCSC），这里选用NCBI的GRCh38： 版本较新，并且在google clould的genomic public data中有与之相对应的vcf参考数据集，方便后续的模型训练及BQSR与VQSR操作
#Reference：
https://zhuanlan.zhihu.com/p/112592962

GATK bundle：
https://console.cloud.google.com/storage/browser/genomics-public-data/resources/broad/hg38/v0;tab=objects?prefix=&forceOnObjectsSortingFiltering=false

#dbsnp_138.hg38.vcf.gz：研究者们发表了相关文章提交上来的变异，很多并没有做过严格的验证，甚至很多是假的，所以并不十分可信，这个数据集不可以作为模型的训练集，更不可以作为真集看待，唯一的用处就是标注我们的数据测得的变异中哪些是已知的变异

#hapmap_3.3.hg38.vcf.gz：数据来自国际人类单倍体型图计划，包含大量家系数据，有严格的质控和严密的实验验证，准确性目前公认最高

#1000G_omni2.5.hg38.vcf.gz：数据来源于Illumina的Omni基因型芯片，它的验证结果常常作为基因型的金标准，是一个高可信的变异结果

#Mills_and_1000G_gold_standard.indels.hg38.vcf.gz：被专门做过验证的indel数据集，较为可信，但indel并不容易验证，很多时候也都是挑一些比较容易验证的结果

snpdb下载地址：
https://ftp.ncbi.nlm.nih.gov/snp/.redesign/pre_build152/organisms/human_9606_b151_GRCh38p7/VCF/GATK

```

## #整理原始数据：

```shell
#数据地址：
oss://shulab/wangzhang/B20211027-4328_MGGZ20210512A1GUNK-4N/rawdata/
oss://shulab/wangzhang/data-release-A20211118/
COPD_multiomics_metagenome_MGGZ20190625A18A6E-3P\shenzhen

#拷贝、移动：
cp ./* /bigdata/gaojy/cohort1600
mv ./LD21*/*.gz ./
find -perm 777 -exec mv {} /bigdata/gaojy/cohort1600/ \;

#改名：
rename 's/...//' *.gz
#rename用法和sed类似，此方法较笨拙。单纯通过字符个数从前到后删除文件名的内容，直到删够为止（要求文件名长度相同）
#尝试使用正则表达式

#随机抽取指定数量的reads(制造测序深度梯度)
```

##   #数据质控：

```shell
#fastqcq评估序列质量 
nohup time fastqc ./*.fq.gz -t 20
#multiqc汇总评估报告
multiqc -f -d ./ -o ./

#Reference：
https://www.jianshu.com/p/bc7044ade693
https://www.jianshu.com/p/dc6820eb342e
https://zhuanlan.zhihu.com/p/20731723

#获得接头序列
fastp --detect_adapter_for_pe \
	-i ./999.R1.fq.gz \
	-I ./999.R2.fq.gz \
	-o ./999.clean.R1.fq.gz \
	-O ./999.clean.R2.fq.gz \
	-z 4 -q 20 -u 40 -n 6

#Reference：
https://www.cnblogs.com/dataanaly/p/13197991.html
https://www.cxyzjd.com/article/qq_43034844/87788346

#去接头(Cutadapt)
R1 3'：AGATCGGAAGAGCACACGTCTGAACTCCAGTCA
R2 3'：AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT

#conda activate cutadaptenv

#单样本:
cutadapt \
	--discard-trimmed \
	-O 18 \
	-j 20 \
	-b AGATCGGAAGAGCACACGTCTGAACTCCAGTCA \
	-B AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT \
	-o /bigdata/gaojy/cohort1600/cutadapt/999.R1.cutadp.fq.gz \
	-p /bigdata/gaojy/cohort1600/cutadapt/999.R2.cutadp.fq.gz \
	/bigdata/gaojy/cohort1600/999.R1.fq.gz \
	/bigdata/gaojy/cohort1600/999.R2.fq.gz

#多样本：
#!/bin/bash
for i in ./*R1.fq.gz
do
        SAMPLE=$(echo ${i} | sed "s/R1\.fq\.gz//")
        echo ${SAMPLE}R1.fq.gz ${SAMPLE}R2.fq.gz
        cutadapt \
        --discard-trimmed \
        -O 18 \
        -j 20 \
        -b AGATCGGAAGAGCACACGTCTGAACTCCAGTCA \
        -B AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT \
        -o ./cutadapt/${SAMPLE}R1.fq.gz \
        -p ./cutadapt/${SAMPLE}R2.fq.gz \
        ${SAMPLE}R1.fq.gz ${SAMPLE}R2.fq.gz
done

#Reference：
https://cutadapt.readthedocs.io/en/stable/index.html
https://zhuanlan.zhihu.com/p/57037645
https://www.bioinfo-scrounger.com/archives/250/
https://cloud.tencent.com/developer/article/1625189
https://www.cnblogs.com/xudongliang/p/6404958.html
http://www.biotrainee.com/thread-2164-1-1.html

#质控(Trimommatic)
java 
	-jar /bigdata/gaojy/biosoft/trimmomatic/trimmomatic-0.39.jar PE \
	-threads 20 \
	./LR210301-0064_S20210312-0079_C014LWG.R1.fq.gz \
	./LR210301-0064_S20210312-0079_C014LWG.R2.fq.gz \
	-baseout ./LR210301-0064_S20210312-0079_C014LWG.fq.gz \
	ILLUMINACLIP:/bigdata/gaojy/biosoft/trimmomatic/adapters/TruSeq3-PE.fa:2:30:10 \
	LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:50

#Reference：
```

## #序列比对：

```shell
bwa index 
	-a bwtsw 
	/bigdata/gaojy/ref/bwa/GRCh38_chr_LN/GRCh38_chr_LN.fasta 
#为参考基因组(GRCh38)建索引

Reference：
http://starsyi.github.io/2016/05/24/BWA-%E5%91%BD%E4%BB%A4%E8%AF%A6%E8%A7%A3/

#建立dict文件
gatk 
	CreateSequenceDictionary \
	-R /bigdata/gaojy/ref/bwa_samtools/GRCh38_chr_LN/GRCh38_chr_LN.fasta \
	-O /bigdata/gaojy/ref/bwa_samtools/GRCh38_chr_LN/GRCh38_chr_LN.dict 
#为参考基因组(GRCh38)建立dict文件

#序列比对
nohup bwa 
	mem 
	-t 20 
	-M \
	-R "@RG\tID:1004V4\tSM:1004V4\tLB:WGS\tPL:Illumina" \ 
#-R参数内容包含测序信息，必须添加，否则后续会影响多样本合并
	/bigdata/gaojy/ref/bwa_samtools/GRCh38_chr_LN/GRCh38_chr_LN.fasta \
	/bigdata/gaojy/pro_snp/1004V4_1.fastq.gz \
	/bigdata/gaojy/pro_snp/1004V4_2.fastq.gz \
	1> /bigdata/gaojy/pro_snp/1004V4.sam \
	2> /bigdata/gaojy/pro_snp/1004V4.log 
#序列比对到参考基因组（GRCh38）

#!/bin/bash
for i in ./*R1.fq.gz
do
	SAMPLE=$(echo ${i} | sed "s/R1\.fq\.gz//")
	echo ${SAMPLE}R1.fq.gz ${SAMPLE}R2.fq.gz
	bwa mem -t 20 -M \
	-R "@RG\tID:${SAMPLE}\tSM:${SAMPLE}\tLB:WGS\tPL:Illumina" \
	/bigdata/gaojy/ref/bwa_samtools/GRCh38_chr_LN/GRCh38_chr_LN.fasta \
	./${SAMPLE}R1.fq.gz \
	./${SAMPLE}R2.fq.gz \
	| samtools sort -@ 20 -n -O bam \
	-o ${SAMPLE}.sorted.bam
done
#Reference:
```

## **#bam文件解读**

### Header section

```shell
该部分全部以“@”开头，提供基本的软件版本，参考序列信息，排序信息等
@HD:这一行中有各种不同的标识
标识“VN”用以说明格式版本
标识“SO”用以说明比对排序的情况，有unknown (default)、unsorted、queryname和coordinate
对于coordinate：
排序的主键是Alignments section的第三列“RNAME”，其顺序由@SQ行的“SN”标识的顺序定义
次键是Alignments section的第四列“POS”字段。对于RNAME和POS相等的比对，排列顺序则是任意的
@SQ:“SN”标签是参考序列说明，它的值主要是用于Alignments section的第三列“RNAME”和第七列“MRNM”比对的记录
@PG:使用的程序说明；该行“ID”为程序记录标识符，“PN”为程序名字，“CL”为命令行
@CO:任意的说明信息
@RG:比对上的序列（read）说明
```

### Alignment section

```shell
QNAME: 序列的名字，也就是reads的名称
FLAG: 是一个标记的数字，是有需要转换成二进制才能知道代表的意思，各个数字分别代表:
1  read是pair中的一条（read表示本条read，mate表示pair中的另一条read）
2  pair一正一负完美的比对上
4  这条read没有比对上
8  mate没有比对上
16 这条read反向比对
32 mate反向比对
64 这条read是read1
128 这条read是read2
256 第二次比对
512 比对质量不合格
1024 read是PCR或光学副本产生
2048 辅助比对结果
假如说标记不是以上列举出的数字，比如说83=（64+16+2+1），就是这几种情况值和，可以使用二进制数来表示
RNAME: 参考序列的名字
POS: 在参考序列上的位置
MAPQ: mapping qulity 越高则位点越独特，比对的质量值
CIGAR:代表比对结果的CIGAR字符串，如37M1D2M1I，这段字符的意思是37个匹配，1个参考序列上的删除，2个匹配，1个参考序列上的插入。M代表的是alignment match(可以是错配)，可以理解为表示比对的具体情况
M：match/mismatch
I：插入 insertion(和参考基因组相比）
D: 删除 deletion(和参考基因组相比）
N：跳跃 skipped(和参考基因组相比）
S：软剪切 soft clipping ，（表示unaligned,）
H：硬剪切 hard clipping  （被剪切的序列不存在于序列中）
P：填充  padding(表示参考基因组没有，而reads里面含有位点
RNEXT: mate 序列所在参考序列的名称，mate一般指大的片段序列
PNEXT: mate 序列在参考序列上的位置
TLEN: 估计出的片段的长度，当mate 序列位于本序列上游时该值为负值。
SEQ: read的序列
QUAL: read序列对应的ASCII码格式的碱基质量值
```

## #基因分型前处理：

```shell
#转换成二进制
samtools 
	view 
	-b 
	-S 
	-@ 20 
	/bigdata/gaojy/pro_snp/1004V4.sam > /bigdata/gaojy/pro_snp/1004V4.bam 

#查看比对结果
samtools 
	flagstat 
	-@ 20 1004V4.bam > 1004V4.txt

#增加头文件
nohup gatk 
	AddOrReplaceReadGroups 
	-I=1040V4.bam 
	-O=1040V4C.bam 
	--RGID=1040V4 
	--RGSM=1040V4 
	--RGPU=1040V4 
	--RGLB=WGS 
	--RGPL=ILLUMINA

#排序
samtools 
	sort 
	-@ 20 
	-n 
	-O bam 
	-o /bigdata/gaojy/pro_snp/1004V4.sorted.bam 
	/bigdata/gaojy/pro_snp/1004V4.bam

#标记重复（GATK）
gatk 
	MarkDuplicates \
	-I /bigdata/gaojy/pro_snp/1004V4.sorted.bam \
	-O /bigdata/gaojy/pro_snp/1004V4.sorted.marked.bam \
	-M /bigdata/gaojy/pro_snp/1004V4.sorted.marked_metrics.txt
#GATK中的标记重复这一步过于消耗资源，因此考虑使用samtools中的相应功能进行标记

#SamTools标记pcr重复，共3步
samtools 
	fixmate 
	-m 
	-@ 20 
	/bigdata/gaojy/pro_snp/1004V4.sorted.bam 
	/bigdata/gaojy/pro_snp/1004V4.sorted.fixmate.bam
#添加mate score用于后续的标记重复

samtools 
	sort 
	-@ 20 
	-o /bigdata/gaojy/pro_snp/1004V4.positionsorted.bam 
	/bigdata/gaojy/pro_snp/1004V4.sorted.fixmate.bam
#对添加mate score的.bam文件做排序

samtools 
	markdup 
	-@ 20 
	/bigdata/gaojy/pro_snp/1004V4.positionsorted.bam 
	/bigdata/gaojy/pro_snp/1004V4.markdup.bam
#标记重复

#建索引.bai
samtools 
	index 
	1004V4.markdup.bam
#为标记重复的bam文件建立索引.bai文件用于后续的分析

#BaseRecal
nohup gatk 
	BaseRecalibrator \
	-I /bigdata/gaojy/pro_snp/1004V4.markdup.bam \
	-R /bigdata/gaojy/ref/bwa_samtools/GRCh38_chr_LN/GRCh38_chr_LN.fasta \
	-O /bigdata/gaojy/pro_snp/1004V4.BQSR.recal.table \
	--known-sites /bigdata/gaojy/ref/gatk/dbsnp_138.hg38.vcf.gz \
	--known-sites /bigdata/gaojy/ref/gatk/hapmap_3.3.hg38.vcf.gz \
	--known-sites /bigdata/gaojy/ref/gatk/1000G_omni2.5.hg38.vcf.gz \
	--known-sites /bigdata/gaojy/ref/gatk/hg38_v0_1000G_phase3_v4_20130502.sites.hg38.vcf 
#计算出了所有需要进行重校正的read和特征值，然后把这些信息输出为一份校准表文件，用于接下来的BQSR。主要是通过机器学习的方法构建测序碱基的错误率模型，然后对这些碱基的质量值进行相应的调整

#BQSR
nohup gatk 
	ApplyBQSR \
	-bqsr /bigdata/gaojy/pro_snp/1004V4.BQSR.recal.table \
	-R /bigdata/gaojy/ref/bwa_samtools/GRCh38_chr_LN/GRCh38_chr_LN.fasta \
	-I /bigdata/gaojy/pro_snp/1004V4.markdup.bam \
	-O /bigdata/gaojy/pro_snp/1004V4.markdup.BQSR.bam
#利用第一步得到的校准表文件重新调整原来BAM文件中的碱基质量值，并使用这个新的质量值重新输出一份新的BAM文件同时生成索引文件.bai
```

## #基因分型：

```shell
#Haplotype
nohup gatk 
	HaplotypeCaller \
	-R /bigdata/gaojy/ref/bwa_samtools/GRCh38_chr_LN/GRCh38_chr_LN.fasta \
	-ERC GVCF \
	--native-pair-hmm-threads 20 \
	-I /bigdata/gaojy/pro_snp/1004V4.markdup.BQSR.bam \
	-O /bigdata/gaojy/pro_snp/1004V4.markdup.BQSR.g.vcf.gz
#生成单倍体件及其.index用于接下来的变异检测及之后的质控（耗时）

#制作文件列表
find ./ -name "*.g.vcf.gz" > input.list
#生成多样本单倍体文件列表

#合并
nohup gatk 
	CombineGVCFs \
	-R /bigdata/gaojy/ref/bwa_samtools/GRCh38_chr_LN/GRCh38_chr_LN.fasta \
	-V ./input.list
	-O ./COPD.g.vcf.gz

#Genotype（vcf）
nohup gatk 
	GenotypeGVCFs \
	-R /bigdata/gaojy/ref/bwa_samtools/GRCh38_chr_LN/GRCh38_chr_LN.fasta \
	-dbsnp /bigdata/gaojy/ref/gatk/dbsnp_138.hg38.vcf.gz\
	-V /bigdata/gaojy/pro_snp/COPD.g.vcf.gz \
	-O /bigdata/gaojy/pro_snp/COPD.genotype.vcf.gz
#检测并标记出变异位点信息及rs号

#vcf文件解读：
https://blog.csdn.net/genome_denovo/article/details/78697679
https://blog.csdn.net/u012150360/article/details/70666213
https://www.jianshu.com/p/754f6f2bd061
```

## **#vcf文件解读**

![image-20211220170332543](C:\Users\Vincent_Bioinfo\AppData\Roaming\Typora\typora-user-images\image-20211220170332543.png)

```shell
1.CHROM：记录染色体编号
2.POS：记录变异位点在参考基因组中的位置。如果是SNP的话，POS即SNP的位置；如果是INDEL的话，位置是INDEL的第一个碱基位置。
3.ID：SNP/INDEL的ID, 如在dbSNP中有该SNP的id，则会在此行给出；若没有，则用'.'表示其为一个novel variant, dbSNP编号通常以rs开头，一般只有人类基因组才有dbSNP编号；INDEL指的是在基因组的某个位置上所发生的small deletion,small inverion小片段序列的插入或者删除，其长度通常在50bp以下。
4.REF：参考基因组该位置碱基类型，必须是A,C,G,T,N    N表示不确定碱基，SNP应该一个位点就是一个碱基。
5.ALT：与参考序列比较，发生突变的变异碱基类型，必须是A,C,G,T,N，.    多个用逗号分割。"." 表示这个地方没有reads覆盖为缺失。
6.QUAL：变异位点检测质量值，越高越可靠。表示在该位点存在variant的可能性，该值越高，则variant的可能性越大。
7.FILTER：如果该位点通过过滤标准那么我们可以在该列标记为"PASS",说明该列质量值高
```

## #基因型VQSR:

```shell
#SNP
nohup gatk 
	SelectVariants \
	-R /bigdata/gaojy/ref/bwa_samtools/GRCh38_chr_LN/GRCh38_chr_LN.fasta \
	--select-type-to-include SNP \
	-V /bigdata/gaojy/pro_snp/1004V4.markdup.BQSR.genotype.vcf \
	-O /bigdata/gaojy/pro_snp/1004V4.markdup.BQSR.genotype.snp.vcf

#Indel
nohup gatk 
	SelectVariants \
	-R /bigdata/gaojy/ref/bwa_samtools/GRCh38_chr_LN/GRCh38_chr_LN.fasta \
	--select-type-to-include INDEL \
	-V /bigdata/gaojy/pro_snp/1004V4.markdup.BQSR.genotype.vcf \
	-O /bigdata/gaojy/pro_snp/1004V4.markdup.BQSR.genotype.indel.vcf

#VariantRecal_SNP
nohup gatk 
	VariantRecalibrator \
	-R /bigdata/gaojy/ref/bwa_samtools/GRCh38_chr_LN/GRCh38_chr_LN.fasta \
	-V /bigdata/gaojy/pro_snp/1004V4.markdup.BQSR.genotype.snp.vcf \
	--resource:hapmap,known=false,training=true,truth=true,prior=15.0 \
	/bigdata/gaojy/ref/gatk/hapmap_3.3.hg38.vcf.gz \
	--resource:omini,known=false,training=true,truth=false,prior=12.0 \
	/bigdata/gaojy/ref/gatk/1000G_omni2.5.hg38.vcf.gz \
	--resource:dbsnp,known=true,training=false,truth=false,prior=2.0 \
	/bigdata/gaojy/ref/gatk/dbsnp_138.hg38.vcf.gz \
	--resource:1000g,known=false,training=true,truth=false,prior=10.0
	bigdata/gaojy/ref/gatk/hg38_v0_1000G_phase3_v4_20130502.sites.hg38.vcf \
	-an DP 
	-an QD 
	-an FS 
	-an MQ 
	-an SOR 
	-an BaseQRankSum 
	-an ReadPosRankSum 
	-an MQRankSum \
	-mode SNP \
	-tranche 100.0 -tranche 99.9 -tranche 99.0 -tranche 95.0 -tranche 90.0 \
	--rscript-file /bigdata/gaojy/pro_snp/1004V4.markdup.BQSR.genotype.snp.plots.R \
	--tranches-file /bigdata/gaojy/pro_snp/1004V4.markdup.BQSR.genotype.snp.tranches \
	--tmp-dir /bigdata/gaojy/tmp \
	-O /bigdata/gaojy/pro_snp/1004V4.markdup.BQSR.genotype.snp.recal.table \
#校准genotype.snp.vcf，生成校准表文件.recal

#VQSR_SNP
nohup gatk 
	ApplyVQSR \
	-R /bigdata/gaojy/ref/bwa_samtools/GRCh38_chr_LN/GRCh38_chr_LN.fasta \
	-V /bigdata/gaojy/pro_snp/1004V4.markdup.BQSR.genotype.snp.vcf \
	--tranches-file /bigdata/gaojy/pro_snp/1004V4.markdup.BQSR.genotype.snp.tranches \
	--recal-file /bigdata/gaojy/pro_snp/1004V4.markdup.BQSR.genotype.snp.recal.table \
	-mode SNP \
	-O /bigdata/gaojy/pro_snp/1004V4.markdup.BQSR.genotype.snp.VQSR.vcf \

#VariantRecal_Indel
nohup gatk 
	VariantRecalibrator \
	-R /bigdata/gaojy/ref/bwa_samtools/GRCh38_chr_LN/GRCh38_chr_LN.fasta \
	-V /bigdata/gaojy/pro_snp/1004V4_GRCh38_chr_LN.markdup.recal.BQSR.genotype.indel.vcf \
	-resource:mills,known=true,trainig=true,truth=true,prior=12.0
	/bigdata/gaojy/ref/gatk/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz \
	-an DP 
	-an QD 
	-an FS 
	-an SOR 
	-an ReadPosRankSum 
	-an MQRankSum \
	-mode INDEL \
	--rscript-file /bigdata/gaojy/pro_snp/1004V4.markdup.recal.BQSR.genotype.indel.plots.R \
	--tranches-file /bigdata/gaojy/pro_snp/1004V4.markdup.recal.BQSR.genotype.indel.tranches \
	-O /bigdata/gaojy/pro_snp/1004V4.markdup.recal.BQSR.genotype.indel.recal
#原理同snp，只是训练集换成了针对indel的Mills_and_1000G_gold_standard.indels.hg38.vcf.gz

#VQSR_Indel
nohup gatk 
	ApplyVQSR \
	-R /bigdata/gaojy/ref/bwa_samtools/GRCh38_chr_LN/GRCh38_chr_LN.fasta \
	-V /bigdata/gaojy/pro_snp/1004V4.markdup.recal.BQSR.genotype.indel.vcf \
	--tranches-file /bigdata/gaojy/pro_snp/1004V4.markdup.recal.BQSR.genotype.indel.tranches \
	--recal-file /bigdata/gaojy/pro_snp/1004V4.markdup.recal.BQSR.genotype.indel.recal \
	-mode SNP \
	-O /bigdata/gaojy/pro_snp/1004V4_GRCh38_chr_LN.markdup.recal.BQSR.genotype.indel.VQSR.vcf
```

## #数据清洗：

```shell
sed 		'/_random/d;/chrUn/d;/_alt/d;/chrX/d;/chrY/d;/chrM/d;/_decoy/d;/chrEBV/d;/##bcftools_merge/d;/##GATKCommandLine/d;/HLA/d' 
cohort135.vcf > cohort.washed.vcf
#摘取目标染色体chr1~chr22

grep '^#' merge.vcf > header.vcf
grep 'PASS' merge.vcf >> header.vcf
#过滤出PASS的位点

nohup gatk 
	SortVcf 
	--TMP_DIR /bigdata/gaojy/tmp/ \
	-I ./cohort135.snp.VQSR.washed.vcf \
	-O ./cohort135.snp.VQSR.washed.sorted.vcf
#排序
```

## #Plink质控：

```shell
plink 
	--vcf cohort135.arranged.vcf 
	--recode 
	--allow-no-sex 
	--out cohort135.arranged
#转换格式

plink 
	--file chr1 
	--missing 
	--out chr1
#检测缺失率

plink 
	--file cohort135 
	-geno 0.05  
	--recode 
	--out cohort135.geno
#保留检出率>95%的位点（该snp在95%以上的样本中存在）

plink 
	--file cohort135 
	-mind 0.05  
	--recode 
	--out cohort135.geno.mind
#保留检出率>95%的样本(该样本含有的snp占总snp的95%以上)

plink 
	--file cohort135 
	-maf 0.05  
	--recode 
	--out cohort135.geno.mind.maf
#保留次要等位基因频率>0.05的变异位点

plink 
	--file cohort135 
	-hwe 0.0001  
	--recode
	--out cohort135.geno.mind.maf.hwe
#删除不在HWE中的位点:hwe<0.0001

#Reference:
https://zhuanlan.zhihu.com/p/347958906
https://dengfei2013.gitee.io/plink-cookbook/plink-%E8%BD%AF%E4%BB%B6%E4%BB%8B%E7%BB%8D.html
```

## #转换坐标系:

```shell
bgzip -c C132V8H.vcf > C132V8H.vcf.gz
tabix -p vcf K147V0H.vcf.gz
#压缩建索引
#这里一定要用bgzip，不能用gzip，否则拆分时会报错

for i in {1..22}
do
	/data/gaojy/biosoft/bcftools-1.8/bcftools view \
	-r chr${i} cohort135.arranged.vcf.gz \
	-o chr${i}.vcf.gz \
	-O z
done
#拆分染色体

#转换坐标：
#首先制作.bed文件，crossmap只接受.bed文件，且只定义前三列：chr、start、end，文件无表头：
cat chr1.vcf | sed -e '/^#/d' > chr1.noheader.vcf
sed -i 's/^/chr&/' chr1.noheader.vcf #添加chr
awk '{print $1, $2, $2+1, $3}' chr1.noheader.vcf > chr1.bed

#为.bed文件添加行号，方便将转换后的位点与.vcf信息相对应
awk '$0=NR" "$0' chr1.bed > chr1.NR.bed
awk '{print $2, $3, $4, $5, $1}' chr1.NR.bed > chr1.bed

#使用crossmap进行坐标转换：
conda activate py36
which CrossMap.py
python /bigdata/gaojy/miniconda2/envs/py36/bin/CrossMap.py \ 
	bed \ 
	../GRCh38_to_GRCh37.chain.gz \ 
	chr1.bed > chr1.37.bed

#挑出成功转换的位点：
awk '{if($7=="chr1"){print $0}}' chr1.37.bed > chr1.hg37.bed 
awk '{print $5}' chr1.hg37.bed > chr1.hg37.position.txt 
awk '$0=NR" "$0' chr1.noheader.vcf > chr1.noheader.NR.vcf 

awk 'NR==FNR{a[$1];next}$1 in a' chr1.hg37.position.txt chr1.noheader.NR.vcf > chr1.noheader.ne
w.NR.vcf 
#NR==FNR，表示将要处理第一个文件
#a[$1]，以第一个文件第1列为下标，建立数组a
#next，表示要跳过第2个文件，不执行a[$1]
#$1 in a，判断如果第2个文件中的第1列存在于数组a就会打印输出
#根据位置文件筛选出在.vcf文件中缺失率小于0.1的snp位点

#接着需要给新的.vcf文件附上新的坐标并删除行号
awk '{print $8}' chr1.hg37.bed > chr1.hg37.newposition.txt 
awk 'FNR==NR{a[NR]=$1;next}{$3=a[FNR]}1' chr1.hg37.newposition.txt chr1.noheader.new.NR.vcf > ch
r1.hg37.noheader.vcf 
awk '{$1="";print $0}' chr1.hg37.noheader.vcf > chr1.hg37.noheader.noNR.vcf 

#加入header
grep "^#" chr1.vcf > chr1.hg37.vcf 
cat chr1.hg37.noheader.noNR.vcf >> chr1.hg37.vcf 

#将各列用制表符隔开
awk 'BEGIN{ FS=" ";OFS="\t" }{ print $1,$2,$3,$4,....}' chr1.hg37.vcf > chr1.hg37.tab.vcf 

#删掉存在“*”的行
sed -e '/*/d' chr1.hg37.tab.vcf > chr1.hg37.ready.vcf 
```

## #构建EAS参考基因集

```shell
#获得最新版本的千人基因组.vcf格式文件(n=3202)
#转格式
for i in {1..22}
do
	/data/gaojy/biosoft/plink_linux_x86_64/plink 
	--vcf chr${i}.filtered.shapeit2.phased.vcf 
	--recode 
	--allow-no-sex 
	--out chr${i}.filtered.shapeit2.phased
done

#转二进制
for i in {1..22}
do
	/data/gaojy/biosoft/plink_linux_x86_64/plink 
	--file chr${i}.filtered.shapeit2.phased 
	--make-bed 
	--out chr${i}.filtered.shapeit2.phased
done

#提取EAS(n=463)
for i in {1..22}
do
	/data/gaojy/biosoft/plink_linux_x86_64/plink 
	--bfile chr${i}.filtered.shapeit2.phased 
	--noweb 
	--keep eas.txt 
	--recode 
	--make-bed 
	--out chr${i}.filtered.shapeit2.phased.eas
done

#转换格式为.vcf
for i in {1..22}
do
/data/gaojy/biosoft/plink_linux_x86_64/plink 
	--bfile chr${i}.filtered.shapeit2.phased.eas 
	--recode vcf-iid 
	--out chr${i}.filtered.shapeit2.phased.eas
done

#统计各染色体snp情况
for i in {1..22}
do
	/data/gaojy/biosoft/bcftools-1.8/bcftools stats 
	chr${i}.filtered.shapeit2.phased.vcf > chr${i}.filtered.shapeit2.phased.stats
done

#取交集
for i in {1..22}
do
	/data/gaojy/biosoft/bcftools-1.8/bcftools isec 
	chr${i}.filtered.shapeit2.phased.eas.vcf.gz 
	chr${i}.filtered.shapeit2.phased.vcf.gz 
	-p ./chr${i}
done

```



## #定相填补：

```shell
1):Beagle：
java 
	-jar /bigdata/gaojy/biosoft/beagle/beagle.28Jun21.220.jar \ 
	ref=/bigdata/gaojy/biosoft/beagle/chr1.1kg.phase3.v5a.vcf.gz \ 
	map=/bigdata/gaojy/biosoft/beagle/plink.chr1.GRCh37.map \ 
	gt=chr1.hg37.qc.vcf \ 
	out=chr1.hg37.imputed \ 
	nthreads=10 \ 
	ne=20000 \ 
	window=100 \ 
	seed=-99999

#得到22条填补后的染色体文件，将他们合并
cat chr1.hg37.imputed | sed -e '/^#/d' > chr1.hg37.imputed.noheader.vcf
sed -i 's/^/chr&/' chr1.hg37.imputed.noheader.vcf 
cat {........} > cohort135.hg37.imputed.noheader.vcf

#接着进行post-imputation-qc,保留DR2>0.7的snp
#首先提取出INFO列
awk '{print $8}' cohort135.hg37.imputed.noheader.vcf > INFO.csv 
awk 'BEGIN{ FS=";";OFS="\t" }{ print $1,$2,$3}' INFO.csv > INFO.tab.csv 
awk '{print $1}' INFO.tab.csv > DR2.csv 
sed "s/DR2=//" DR2.csv > DR2.txt 
paste DR2.txt cohort135.hg37.imputed.noheader.vcf > cohort135.hg37.imputed.noheader.DR.vcf 

#删除含有多个DR2值的行
sed -e '/,/d' cohort135.hg37.imputed.noheader.DR.vcf > cohort135.hg37.imputed.noheader.DR.singl
e.vcf 

#筛选出DR2>0.7的snp
awk '$1>0.7{print $0}' cohort135.hg37.imputed.noheader.DR.single.vcf > cohort135.hg37.imputed.no
header.DR.single.0.7.vcf

#接着我们使用PLINK对DR质控后的snp进行进一步的post-imputation qc
cut -f 2-160 cohort135.hg37.imputed.noheader.DR.single.0.7.vcf > cohort135.hg37.imputed.noheade
r.DR.single.0.7.clean.vcf 
cat header.vcf cohort135.hg37.imputed.noheader.DR.single.0.7.clean.vcf > cohort135.hg37.imputed.
DR.0.7.vcf 

2):Shapeit+Impute:

3):GLIMPSE_v2:

1.Extract sites from the reference panel
/data/gaojy/biosoft/bcftools-1.8/bcftools view 
	-G 
	-Oz 
	-o ./1000GP.chr22.noNA12878.sites.vcf.gz 
	./1000GP.chr22.noNA12878.bcf

/data/gaojy/biosoft/bcftools-1.8/bcftools index 
	-f ./1000GP.chr22.noNA12878.sites.vcf.gz

2. Split the genome into chunks
/data/gaojy/biosoft/GLIMPSE-2.0.0/tutorial/bin/GLIMPSE2_chunk_static 
	--input /data/gaojy/ref/glimpse2/chr22.filtered.shapeit2.phased.vcf.gz 
	--region chr22 
	--sequential 
	--threads 20 
	--output chunks.chr22.txt

3. Create binary reference panel
VCF=/data/gaojy/gwas/cohort1644/hg38/chr22.vcf
REF=/data/gaojy/ref/glimpse2/chr22.filtered.shapeit2.phased.vcf.gz
MAP=/data/gaojy/ref/glimpse2/chr22.b38.gmap.gz

while IFS="" read -r LINE || [ -n "$LINE" ];
do
	printf -v ID "%02d" $(echo $LINE | cut -d" " -f1)
	IRG=$(echo $LINE | cut -d" " -f3)
	ORG=$(echo $LINE | cut -d" " -f4)
	/data/gaojy/biosoft/GLIMPSE-2.0.0/tutorial/bin/GLIMPSE2_split_reference_static \
	--reference ${REF} \
	--map ${MAP} \
	--input-region ${IRG} \
	--output-region ${ORG} \
	--output ./chr22.imputed
done < chunks.chr22.txt

bgzip -c chr22.vcf > chr22.vcf.gz
tabix -p vcf chr22.vcf.gz
#养成建索引的好习惯

4. Running GLIMPSE2
REF=/data/gaojy/gwas/cohort1644/hg38/chr22.imputed
VCF=/data/gaojy/gwas/cohort1644/hg38/chr22.vcf.gz

while IFS="" read -r LINE || [ -n "$LINE" ]; 
do   
	printf -v ID "%02d" $(echo $LINE | cut -d" " -f1)
	IRG=$(echo $LINE | cut -d" " -f3)
	ORG=$(echo $LINE | cut -d" " -f4)
	CHR=$(echo ${LINE} | cut -d" " -f2)
	REGS=$(echo ${IRG} | cut -d":" -f 2 | cut -d"-" -f1)
	REGE=$(echo ${IRG} | cut -d":" -f 2 | cut -d"-" -f2)
	OUT=./NA12878_imputed
	/data/gaojy/biosoft/GLIMPSE-2.0.0/tutorial/bin/GLIMPSE2_phase_static 
	--bam-file ${BAM} 
	--reference ${REF}_${CHR}_${REGS}_${REGE}.bin 
	--output ${OUT}_${CHR}_${REGS}_${REGE}.bcf
done < chunks.chr22.txt

5. Ligate chunks of the the same chromosome

LST=list.chr22.txt
ls -1v chr22_imputed_*.bcf > ${LST}

OUT=chr22.ligated.bcf
/data/gaojy/biosoft/GLIMPSE-2.0.0/tutorial/bin/GLIMPSE2_ligate_static --input ${LST} --output $OUT

#转回.vcf
/data/gaojy/biosoft/bcftools-1.8/bcftools view chr22.ligated.bcf > chr22.ligated.vcf

#QC
......

4):QUILT:



```

## **#填补后质控**：

```shell
plink 
	--vcf cohort135.hg37.imputed.DR.0.7.vcf 
	--recode 
	--allow-no-sex 
	--out cohort135.hg37.imputed.DR.0.7 
plink 
	--file cohort135.hg37.imputed.DR.0.7 
	--maf 0.05 
	--recode 
	--out cohort135.hg37.imputed.DR.0.7.maf0.05 

#删掉离群的两个样本：
plink 
	--bfile cohort135.hg37.imputed.DR.0.7.maf0.05 
	--remove sample.txt 
	--recode 
	--make-bed 
	--out ./sampleQC/cohort135.hg37.imputed.DR.0.7.maf0.05.sampleQC 
```

##   #基因型—表型关联:

```shell
#首先计算每个样本的前十个PCA值
plink 
	--bfile cohort135.hg37.imputed.DR.0.7.maf0.05 
	--indep-pairwise 50 5 0.2 
	--out cohort135.hg37.imputed.DR.0.7.maf0.05.indepSNP 
plink 
	--bfile cohort135.hg37.imputed.DR.0.7.maf0.05 
	--extract cohort135.hg37.imputed.DR.0.7.maf0.05.indepSNP.prune.in 
	--pca 10 
	--out cohort135.hg37.imputed.DR.0.7.maf0.05.indepSNP.pca 
	
#GWAS中系统偏差的一个重要来源是种群分层。已经表明， 即使是单一种群个体中，也可能存在细微程度的分层。因此，测试和控制种群分层的存在是不可或缺的质量控制步骤。有多种方法可以校正种群分层。PLINK中包含的一种方法：多维缩放（multimensional scaling, MDS）法。该方法计算样本中任何一对个体之间共享的等位基因的全基因组平均比例，以生成每个个体遗传变异的定量指数（成分，components）。可以绘制单个成分的分数以探索在遗传意义上是否存在比预期更相似的个体组。	
#主成分可以被视为反映 由于祖先不同而造成的样本遗传变异的连续方差轴。对于某一主成分，有相似值的个体有着相似的对应祖先

1):PLINK:
#纳入协变量后，直接进行关联
plink 
	--bfile cohort135.hg37.imputed.DR.0.7.maf0.05 
	--linear 
	--pheno 135.txt 
	--pheno-name FEV1/FVC 
	--covar 135co.txt 
	--allow-no-sex 
	--out cohort135.hg37.imputed.DR.0.7.maf0.05.linear.FVC 

2):GEMMA:
#首先将fam文件中最后一行的表型值换成真正的表型值（因为在gemma0.94和0.97中缺失表型值无法进行亲缘关系矩阵的计算）
#计算亲缘关系矩阵：
/bigdata/gaojy/biosoft/gemma/gemma 
	-bfile cohort135.hg37.imputed.DR.0.7.maf0.05 
	-gk 2 
	-o cohort135.hg37.imputed.DR.0.7.maf0.05

#接着进行关联分析：
/bigdata/gaojy/biosoft/gemma/gemma 
	-bfile ./cohort135.hg37.imputed.DR.0.7.maf0.05 
	-k ./cohort135.hg37.imputed.DR.0.7.maf0.05.sXX.txt 
	-lmm 4 
	-c ./cohort135.co.txt 
	-o cohort135.hg37.imputed.DR.0.7.maf0.05 
#注：GEMMA的协变量文件需要无表头

3):MicrobiomeGWAS:
#首先使用R中vegan包的vegdist函数计算样本间的pray-curtis距离矩阵
#输入文件：.bed, .bim, .fam 以及距离矩阵文件.txt, 协变量文件.co.txt
Rscript /bigdata/gaojy/biosoft/microbiomeGWAS/R/microbiomeGWAS.R 
	-r /bigdata/gaojy/biosoft/microbiomeGWAS/ 
	-p ./cohort135.hg37.imputed.DR.0.7.maf0.05 
	-d ./135.bray.txt 
	-c ./135.co.txt  
```

## #基因注释：

```shell

```

# **Post-GWAS**

## #孟德尔随机化：

```R
library(TwoSampleMR) #加载R包
copd <- extract_instruments(outcomes='bbj-a-103', clump=TRUE, r2=0.02, kb=5000, access_token=NULL) #获取暴露数据
head(copd) #查看暴露数据
dim(copd) #查看IV个数

#计算PVE
submicro$pve = (2*(submicro$beta^2*submicro$af*(1-submicro$af)))/(2*submicro$beta*submicro$af*(1-submicro$af) + submicro$se^2*2*97*submicro$af*(1-submicro$af))#计算PVE
write.table(submicro,"D:/R_data/mr/taxa/pve.txt",sep="\t",quote= F,col.names = NA)
pve = read.table('D:/R_data/mr/taxa/pve.txt', header=1)

#计算F统计值
pve$F = pve$pve*(97-2)/(1-pve$pve)
pve
write.table(pve,"D:/R_data/mr/taxa/F.txt",sep="\t",quote= F,col.names = NA)

#读入结局文件
t2d = read.table('D:/R_data/mr/Streptococcus_intermedius.txt.gz', header=1)

# 添加phenotype列
t2d$phenotype <- 'Streptococcus_intermedius' 

#提取IV在结局中的信息
t2d_out <- format_data(dat=t2d,
                       type = "outcome",
                       snps = copd$SNP,
                       header = TRUE,
                       phenotype_col = "phenotype",
                       snp_col = "rs",
                       beta_col = "beta",
                       se_col = "se",
                       eaf_col = "af",
                       effect_allele_col = "allele1",
                       other_allele_col = "allele0",
                       pval_col = "p_wald")
dim(t2d_out)

mydata <- harmonise_data(exposure_dat=copd, outcome_dat=t2d_out, action= 2)
#将IV的效应等位基因(effect allele)对齐
dim(mydata)
mydata
write.table(mydata,"D:/R_data/mr/IV.txt",sep="\t",quote= F,col.names = NA)

res <- mr(mydata)#计算MR
res
dim(res)
write.table(res,"D:/R_data/mr/taxa/Strep.COPD.EAS.mr.5.txt",sep="\t",quote= F,col.names = NA)#保存结果

het <- mr_heterogeneity(mydata)#异质性检验
het
mr_method_list()
mr(mydata,method_list=c('mr_ivw_fe')) #使用固定效应模型

pleio <- mr_pleiotropy_test(mydata)#多效性检验
pleio
write.table(pleio,"D:/R_data/pleiotropy.txt",sep="\t",quote= F,col.names = NA)

single <- mr_leaveoneout(mydata)#逐个剔除法敏感性检验
mr_leaveoneout_plot(single)

p1 <- mr_scatter_plot(res, mydata)#绘制散点图
p1

methods = all_method=c('mr_ivw', 'mr_egger_regression', 'mr_weighted_mode', 'mr_weighted_median')
res_single <- mr_singlesnp(mydata, all_method=methods)
p2 <- mr_forest_plot(res_single)#绘制森林图
p2

res_single <- mr_singlesnp(mydata)
p4 <- mr_funnel_plot(res_single)#绘制漏斗图
p4

#reverse
microbiome <-read.table('D:/R_data/mr/Streptococcus_intermedius.txt.gz', header=T)#读入暴露文件
submicro <- subset(microbiome,microbiome$p_wald<1e-5)
dim(submicro)
write.table(submicro,"D:/R_data/Stre.sub.txt",sep="\t",quote= F,col.names = NA)

exp_dat <- read_exposure_data(filename = 'D:/R_data/Stre.sub.txt',
                              clump = FALSE,
                              sep= "\t",
                              snp_col = "rs",
                              beta_col = "beta",
                              se_col = "se",
                              eaf_col = "af",
                              effect_allele_col = "allele1",
                              other_allele_col = "allele0",
                              pval_col = "p_wald")

submicro_clump <-clump_data(exp_dat,clump_r2=0.02, clump_kb=5000, pop='EAS')#获取工具变量

t2d_out <- extract_outcome_data(snps=submicro_clump$SNP, outcomes='bbj-a-103', 
                                proxies = FALSE, 
                                maf_threshold = 0.001, 
                                access_token = NULL)#提取IV在结局中的信息
mydata <- harmonise_data(exposure_dat=submicro_clump, outcome_dat=t2d_out, action= 2)
#将IV的效应等位基因（effect allele）对齐
mydata
dim(mydata)
res <- mr(mydata)#计算MR
res
```

## #eQTL\sQTL\mQTL共定位分析：

```shell
eQTL共定位分析属于Post-GWAS的一项重要工作，旨在GWAS结果的基础上鉴定与表型相关的eQTL位点。**传统的GWAS是将全基因组范围内的常见变异进行关联分析，鉴定与表型相关的基因座，但鉴定出来的位点大多数位于基因间隔区，其如何通过基因或者通路影响表型很难被阐述。**基于此，开发了eQTL共定位分析方法。其原理是利用已有数据库公布的eQTL位点，结合GWAS summary数据，鉴定与表型相关的eQTL位点。

表达数量性状基因座是对上述概念的进一步深化，它指的是染色体上一些能特定调控mRNA和蛋白质表达水平的区域，其mRNA/蛋白质的表达水平量与数量性状成比例关系。

eQTL可分为顺式作用eQTL和反式作用eQTL：

- 顺式作用eQTL就是某个基因的eQTL定位到该基因所在的基因组区域，表明可能是该基因本身的差别引起的mRNA水平变化
- 反式作用eQTL是指某个基因的eQTL定位到其他基因组区域，表明其他基因的差别控制该基因mRNA水平的差异

**eQTL就是把基因表达作为一种性状,研究遗传突变与基因表达的相关性**

简单地说, 遗传学研究经常发现一些致病或易感突变, 这些突变怎样导致表型有时候不太直观; 所以用某个基因的差异表达作为过渡: 

**突变A>B ---> 基因表达变化 ---> 表型改变**

然而，从基因的改变到疾病等现象的出现，中间缺失了重要的一环，那就是基因的表达。也许在测序中，我们可以看到某一个基因上某一个位置的变化（比如说SNP单核苷酸变化），**但是这种变化并不一定会影响mRNA的产生或者蛋白的改变。也就有可能不会影响到疾病或其他生物学过程。**于是科学家想到了另一个指标——mRNA的序列数据。因为只有被表转译到mRNA上的基因，才可能进一步表达为蛋白

但是要怎么搞清DNA改变是怎么影响mRNA的出现呢？这一过程被称为Expression quantitative trait loci（eQTL） 分析，目的在于得到单个DNA突变与单个基因表达量之间的相关性。**与单个基因mRNA表达量相关的DNA突变，就被称为eQTL**

我们首先通过全基因组测序获得每个个体的DNA全序，然后以同种族的其他个体作为参照，标记出该个体所有的DNA变异位点， 称为SNP位点。同时，我们通过全基因组mRNA表达量测序得到该个体的特定组织样本中的基因表达量。**以全部DNA变异位点为自变量，轮流以每种mRNA表达量为因变量，用大量的个体数据做样本进行线性回归，就可以得到每一个SNP位点和每一个mRNA表达量之间的关系**
```

**制作GWAS summary文件：**.ma

```shell
zcat Burkholderia_contaminans.txt.gz | awk '{print $2,$5,$6,$7,$8,$9,$13,"NA"}' > Burkholderia_contaminans.eQTL.ma
sed -i '1i\SNP A1 A2 freq b se p n' Campylobacter_concisus.eQTL.ma
sed -i '2d' Campylobacter_concisus.eQTL.ma
```

![image-20221004154900835](C:\Users\Vincent_Bioinfo\AppData\Roaming\Typora\typora-user-images\image-20221004154900835.png)

```shell
/data/gaojy/biosoft/smr-1.3.1-linux-x86_64/smr-1.3.1 
--bfile ../COPD.clean 
--gwas-summary ./Neisseria_meningitidis.eQTL.ma 
--beqtl-summary /data/gaojy/biosoft/smr-1.3.1-linux-x86_64/Lung.lite.eQTL/Lung.lite 
--thread-num 20 
--diff-freq-prop 0.05 
--out ./Neisseria_meningitidis.smr

awk -F'\t' '{if($2==123) print $0}'  data.txt

for chr in {1..22}; do ./command chr$chr.txt; done 
awk -F'\t' '{if($19<0.005&&$20>0.05&&$21>10) print $0}'  149.sQTL.smr > 149.sQTL.smr.PASS.txt

for chr in {1..22}; 
do 
	/data/gaojy/biosoft/smr-1.3.1-linux-x86_64/smr-1.3.1 
	--bfile ../COPD 
	--gwas-summary ./235.eQTL.ma 
	--beqtl-summary /data/gaojy/biosoft/smr-1.3.1-linux-x86_64/Lung.lite.sQTL/sQTL_Lung.lite.chr$chr 
	--out ./plot/ENSG00000116288.12
	--plot
	--probe chr1:7961735:7962763:clu_55378:ENSG00000116288.12 
	--probe-wind 500 
	--gene-list /data/gaojy/biosoft/smr-1.3.1-linux-x86_64/glist-hg19; 
done
```

## **#GCTA遗传力分析：**

"GWAS研究中发现的显著SNP只能解释人类群体中**复杂性状很小一部分**的遗传变异。那么剩余的遗传力在哪里？ **很多时候这部分遗传力并没有丢失，而是由于部分snp效应太小以至于无法达到显著水平而没有被检测到"**

"实际操作中我们首先需要**排除掉有亲缘关系的个体**，主要原因是该模型的目的是估计所有SNP解释的方差，如果纳入有亲缘关系的个体，会由于表型关联（ phenotypic correlations ）而造成偏差，例如由于共同的环境因素。即使没有上述的偏差，那么对它的解读也会有别于 “无关联的”个体：一个基于家系的估计值捕捉的是**所有因果变异的贡献**，而该方法捕捉的则是**与被基因分型的SNP处于LD的因果变异的贡献**"

heritability,翻译为遗传力， 用来描述表型变异中遗传变异的比例。众所周知，表型(P)由基因型(G)和环境因素(E)共同控制

即：**P = G  + E**

**遗传力就是基因G所占的比例**，具体的，通过方差来描述遗传变异和表型变异，则遗传力的公式如下

![img](https://cdn.nlark.com/yuque/0/2022/png/21865910/1650446466789-b1ec88dd-2018-4e18-8521-ebdb7f3a4cfb.png)

“分子为一组样本基因型的方差，分母为表型的方差。方差表征的是一组样本的离散程度，所以遗传力是一个针对群体的概念，通过该公式计算出来的遗传力也称之为广义遗传力”

```shell
GCTA估计GRM：
/bigdata/gaojy/biosoft/gcta/gcta 
	--bfile cohort135.hg37.imputed.DR.0.7.maf0.05.sampleQC 
	--autosome 
	--maf 0.05 
	--make-grm 
	--out cohort135.hg37.imputed.DR.0.7.maf0.05.sampleQC 
	--thread-num 20

GCTA-GREML估计Heritability：
/bigdata/gaojy/biosoft/gcta/gcta 
	--grm cohort135.hg37.imputed.DR.0.7.maf0.05.sampleQC 
	--pheno 133gcta.txt 
	--reml 
	--out cohort135.hg37.imputed.DR.0.7.maf0.05.sampleQC 
	--thread-num 20
```

## **#PRS多基因风险评估：**

传统的GWAS研究只计算单个SNP位点与表型之间的关联性，再用Bonferroni校正，通过给定的阈值，筛选出显著的SNP位点。 这样会存在两个问题：

**第一**：Bonferroni校正非常严格，很多对表型有贡献的位点会因为达不到阈值而被过滤掉

**第二**：单个位点对表型的解释度是很低的，尤其是对于高血压这种多基因控制的表型，用一个个单独的位点解释高血压患病风险，就显得很单薄。

**因此，开发一个能让我们直观的感受患某种疾病的风险多高的工具，显然是非常有必要**











