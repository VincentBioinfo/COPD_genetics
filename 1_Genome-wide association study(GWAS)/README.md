1. Software and database required

Software:

cutadapt:
bwa：https://sourceforge.net/projects/bio-bwa/
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
Database:
