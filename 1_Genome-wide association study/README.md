

# 1. Software and database required

------

## Software:

cutadapt:https://github.com/marcelm/cutadapt/

bwa:https://sourceforge.net/projects/bio-bwa/

samtools:https://github.com/samtools/samtools/releases

jdk:http://www.oracle.com/technetwork/java/javase/downloads/jdk8-downloads-2133151.html

gatk:https://github.com/broadinstitute/gatk/releases/tag/4.1.7.0

plink:http://www.cog-genomics.org/plink2/

tabix:https://sourceforge.net/projects/samtools/files/tabix/tabix-0.2.6.tar.bz2

bcftools:https://github.com/samtools/bcftools/releases/download/1.8/bcftools-1.8.tar.bz2

bedtools:https://github.com/arq5x/bedtools2/archive/v2.25.0.tar.gz

crossmap:https://crossmap.sourceforge.net/

beagle:https://faculty.washington.edu/browning/beagle/b5_2.html#download

gemma:https://github.com/genetics-statistics/GEMMA/releases



## Database:

Human genome(hg38):https://www.ncbi.nlm.nih.gov/genome/51

GATK bundle(hg38):https://console.cloud.google.com/storage/browser/genomics-public-data/resources/broad/hg38/v0;tab=objects?prefix=&forceOnObjectsSortingFiltering=false

Human reference panel for imputation:https://faculty.washington.edu/browning/beagle/beagle.html

# 2.Data preparation

------

```shell
mkdir 00_rawdata 01_sequence_QC 02_alignment 03_genotyping 04_genotype_QC 05_Imputation 06_Association
mv *.fastq 00_rawdata
cd 00_rawdata/
```

# 3.Sequence quality control

------



```shell
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

nohup time fastqc ./*.fq.gz -t 20
multiqc -f -d ./ -o ./
```

# 4.Sequence alignment

------



```shell
bwa index 
	-a bwtsw 
	/bigdata/gaojy/ref/bwa/GRCh38_chr_LN/GRCh38_chr_LN.fasta 

gatk 
	CreateSequenceDictionary \
	-R /bigdata/gaojy/ref/bwa_samtools/GRCh38_chr_LN/GRCh38_chr_LN.fasta \
	-O /bigdata/gaojy/ref/bwa_samtools/GRCh38_chr_LN/GRCh38_chr_LN.dict 

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
```

# 5.Genotyping

------



```shell
samtools view -b -S -@ 20 
	/02_alignment/SAMPLE.sam > /02_alignment/SAMPLE.bam 

samtools sort -@ 20 -n -O bam 
	-o /02_alignment/SAMPLE.sorted.bam 
	/02_alignment/SAMPLE.bam

samtools fixmate -m -@ 20 
	/02_alignment/SAMPLE.sorted.bam 
	/02_alignment/SAMPLE.fixmate.bam

samtools sort -@ 20 
	-o /02_alignment/SAMPLE.positionsorted.bam 
	/02_alignment/SAMPLE.fixmate.bam

samtools markdup -@ 20 
	/02_alignment/SAMPLE.positionsorted.bam 
	/02_alignment/SAMPLE.markdup.bam

samtools index SAMPLE.markdup.bam

nohup gatk 
	BaseRecalibrator \
	-I 02_alignment/SAMPLE.markdup.bam \
	-R /bigdata/gaojy/ref/bwa_samtools/GRCh38_chr_LN/GRCh38_chr_LN.fasta \
	-O /03_genotyping/SAMPLE.BQSR.recal.table \
	--known-sites /bigdata/gaojy/ref/gatk/dbsnp_138.hg38.vcf.gz \
	--known-sites /bigdata/gaojy/ref/gatk/hapmap_3.3.hg38.vcf.gz \
	--known-sites /bigdata/gaojy/ref/gatk/1000G_omni2.5.hg38.vcf.gz \
	--known-sites /bigdata/gaojy/ref/gatk/hg38_v0_1000G_phase3_v4_20130502.sites.hg38.vcf 

nohup gatk 
	ApplyBQSR \
	-bqsr /03_genotyping/SAMPLE.BQSR.recal.table \
	-R /bigdata/gaojy/ref/bwa_samtools/GRCh38_chr_LN/GRCh38_chr_LN.fasta \
	-I /02_alignment/SAMPLE.markdup.bam \
	-O /03_genotyping/SAMPLE.markdup.BQSR.bam

nohup gatk 
	HaplotypeCaller \
	-R /bigdata/gaojy/ref/bwa_samtools/GRCh38_chr_LN/GRCh38_chr_LN.fasta \
	-ERC GVCF \
	--native-pair-hmm-threads 20 \
	-I /03_genotyping/SAMPLE.markdup.BQSR.bam \
	-O /03_genotyping/SAMPLE.g.vcf.gz

find /03_genotyping/ -name "*.g.vcf.gz" > input.list

nohup gatk 
	CombineGVCFs \
	-R /bigdata/gaojy/ref/bwa_samtools/GRCh38_chr_LN/GRCh38_chr_LN.fasta \
	-V /03_genotyping/input.list
	-O /03_genotyping/ALL_SAMPLE.g.vcf.gz

nohup gatk 
	GenotypeGVCFs \
	-R /bigdata/gaojy/ref/bwa_samtools/GRCh38_chr_LN/GRCh38_chr_LN.fasta \
	-dbsnp /bigdata/gaojy/ref/gatk/dbsnp_138.hg38.vcf.gz\
	-V /03_genotyping/ALL_SAMPLE.g.vcf.gz \
	-O /03_genotyping/ALL_SAMPLE.genotype.vcf.gz
```

# 6.Genotype QC

------

```shell
nohup gatk 
	SelectVariants \
	-R /bigdata/gaojy/ref/bwa_samtools/GRCh38_chr_LN/GRCh38_chr_LN.fasta \
	--select-type-to-include SNP \
	-V /03_genotyping/ALL_SAMPLE.genotype.vcf.gz \
	-O /03_genotyping/ALL_SAMPLE.genotype.snp.vcf.gz

nohup gatk 
	VariantRecalibrator \
	-R /bigdata/gaojy/ref/bwa_samtools/GRCh38_chr_LN/GRCh38_chr_LN.fasta \
	-V /03_genotyping/ALL_SAMPLE.genotype.snp.vcf.gz \
	--resource:hapmap,known=false,training=true,truth=true,prior=15.0 \
	/bigdata/gaojy/ref/gatk/hapmap_3.3.hg38.vcf.gz \
	--resource:omini,known=false,training=true,truth=false,prior=12.0 \
	/bigdata/gaojy/ref/gatk/1000G_omni2.5.hg38.vcf.gz \
	--resource:dbsnp,known=true,training=false,truth=false,prior=2.0 \
	/bigdata/gaojy/ref/gatk/dbsnp_138.hg38.vcf.gz \
	--resource:1000g,known=false,training=true,truth=false,prior=10.0
	bigdata/gaojy/ref/gatk/hg38_v0_1000G_phase3_v4_20130502.sites.hg38.vcf \
	-an DP -an QD -an FS -an MQ -an SOR -an BaseQRankSum -an ReadPosRankSum -an MQRankSum \
	-mode SNP \
	-tranche 100.0 -tranche 99.9 -tranche 99.0 -tranche 95.0 -tranche 90.0 \
	--rscript-file /03_genotyping/ALL_SAMPLE.genotype.snp.plots.R \
	--tranches-file /03_genotyping/ALL_SAMPLE.genotype.snp.tranches \
	--tmp-dir /bigdata/gaojy/tmp \
	-O /03_genotyping/ALL_SAMPLE.genotype.snp.recal.table

nohup gatk 
	ApplyVQSR \
	-R /bigdata/gaojy/ref/bwa_samtools/GRCh38_chr_LN/GRCh38_chr_LN.fasta \
	-V /03_genotyping/ALL_SAMPLE.genotype.snp.vcf.gz \
	--tranches-file /03_genotyping/ALL_SAMPLE.genotype.snp.tranches \
	--recal-file /03_genotyping/ALL_SAMPLE.genotype.snp.recal.table \
	-mode SNP \
	-O /04_genotype_QC/ALL_SAMPLE.genotype.snp.VQSR.vcf.gz
	
sed 		'/_random/d;/chrUn/d;/_alt/d;/chrX/d;/chrY/d;/chrM/d;/_decoy/d;/chrEBV/d;/##bcftools_merge/d;/##GATKCommandLine/d;/HLA/d' 
/04_genotype_QC/ALL_SAMPLE.genotype.snp.VQSR.vcf > /04_genotype_QC/ALL_SAMPLE.genotype.snp.VQSR.washed.vcf

grep '^#' /04_genotype_QC/ALL_SAMPLE.genotype.snp.VQSR.washed.vcf > /04_genotype_QC/ALL_SAMPLE.genotype.PASS.vcf
grep 'PASS' /04_genotype_QC/ALL_SAMPLE.genotype.snp.VQSR.washed.vcf >> /04_genotype_QC/ALL_SAMPLE.genotype.PASS.vcf

plink 
	--vcf /04_genotype_QC/ALL_SAMPLE.genotype.PASS.vcf
	--recode 
	--allow-no-sex 
	--out /04_genotype_QC/ALL_SAMPLE.genotype.PASS

plink 
	--file /04_genotype_QC/ALL_SAMPLE.genotype.PASS
	-geno 0.05  
	--recode 
	--out /04_genotype_QC/ALL_SAMPLE.genotype.PASS.geno

plink 
	--file /04_genotype_QC/ALL_SAMPLE.genotype.PASS.geno 
	-mind 0.05  
	--recode 
	--out /04_genotype_QC/ALL_SAMPLE.genotype.PASS.geno.mind

plink 
	--file /04_genotype_QC/ALL_SAMPLE.genotype.PASS.geno.mind
	-maf 0.05  
	--recode 
	--out /04_genotype_QC/ALL_SAMPLE.genotype.PASS.geno.mind.maf

plink 
	--file /04_genotype_QC/ALL_SAMPLE.genotype.PASS.geno.mind.maf 
	-hwe 0.0001  
	--recode
	--out /04_genotype_QC/ALL_SAMPLE.genotype.QC
```

# 7.Imputation

------

```shell
bgzip -c /04_genotype_QC/ALL_SAMPLE.genotype.QC.vcf > /04_genotype_QC/ALL_SAMPLE.genotype.QC.vcf.gz
tabix -p /04_genotype_QC/vcf ALL_SAMPLE.genotype.QC.vcf.gz

for i in {1..22}
do
	/data/gaojy/biosoft/bcftools-1.8/bcftools view \
	-r chr${i} /04_genotype_QC/ALL_SAMPLE.genotype.QC.vcf.gz \
	-o /05_Imputation/chr${i}.vcf.gz \
	-O z
done

cat /05_Imputation/chr1.vcf | sed -e '/^#/d' > /05_Imputation/chr1.noheader.vcf
sed -i 's/^/chr&/' /05_Imputation/chr1.noheader.vcf
awk '{print $1, $2, $2+1, $3}' /05_Imputation/chr1.noheader.vcf > /05_Imputation/chr1.bed
awk '$0=NR" "$0' /05_Imputation/chr1.bed > /05_Imputation/chr1.NR.bed
awk '{print $2, $3, $4, $5, $1}' /05_Imputation/chr1.NR.bed > /05_Imputation/chr1.bed

conda activate py36
which CrossMap.py
python /bigdata/gaojy/miniconda2/envs/py36/bin/CrossMap.py \ 
	bed \ 
	/05_Imputation/GRCh38_to_GRCh37.chain.gz \ 
	/05_Imputation/chr1.bed > /05_Imputation/chr1.37.bed

awk '{if($7=="chr1"){print $0}}' /05_Imputation/chr1.37.bed > /05_Imputation/chr1.hg37.bed 
awk '{print $5}' /05_Imputation/chr1.hg37.bed > /05_Imputation/chr1.hg37.position.txt 
awk '$0=NR" "$0' /05_Imputation/chr1.noheader.vcf > /05_Imputation/chr1.noheader.NR.vcf 

awk 'NR==FNR{a[$1];next}$1 in a' /05_Imputation/chr1.hg37.position.txt /05_Imputation/chr1.noheader.NR.vcf > /05_Imputation/chr1.noheader.new.NR.vcf 

awk '{print $8}' /05_Imputation/chr1.hg37.bed > /05_Imputation/chr1.hg37.newposition.txt 
awk 'FNR==NR{a[NR]=$1;next}{$3=a[FNR]}1' /05_Imputation/chr1.hg37.newposition.txt /05_Imputation/chr1.noheader.new.NR.vcf > /05_Imputation/chr1.hg37.noheader.vcf 
awk '{$1="";print $0}' /05_Imputation/chr1.hg37.noheader.vcf > /05_Imputation/chr1.hg37.noheader.noNR.vcf 

grep "^#" /05_Imputation/chr1.vcf > /05_Imputation/chr1.hg37.vcf 
cat /05_Imputation/chr1.hg37.noheader.noNR.vcf >> /05_Imputation/chr1.hg37.vcf 

awk 'BEGIN{ FS=" ";OFS="\t" }{ print $1,$2,$3,$4}' /05_Imputation/chr1.hg37.vcf > /05_Imputation/chr1.hg37.tab.vcf 

sed -e '/*/d' /05_Imputation/chr1.hg37.tab.vcf > /05_Imputation/chr1.hg37.ready.vcf 

java 
	-jar /bigdata/gaojy/biosoft/beagle/beagle.28Jun21.220.jar \ 
	ref=/bigdata/gaojy/biosoft/beagle/chr1.1kg.phase3.v5a.vcf.gz \ 
	map=/bigdata/gaojy/biosoft/beagle/plink.chr1.GRCh37.map \ 
	gt=/05_Imputation/chr1.hg37.ready.vcf \ 
	out=/05_Imputation/chr1.hg37.imputed \ 
	nthreads=10 \ 
	ne=20000 \ 
	window=100 \ 
	seed=-99999
	
cat /05_Imputation/chr1.hg37.imputed | sed -e '/^#/d' > /05_Imputation/chr1.hg37.imputed.noheader.vcf
sed -i 's/^/chr&/' /05_Imputation/chr1.hg37.imputed.noheader.vcf 

awk '{print $8}' /05_Imputation/All.genotype.snp.hg37.imputed.noheader.vcf > /05_Imputation/INFO.csv 
awk 'BEGIN{ FS=";";OFS="\t" }{ print $1,$2,$3}' INFO.csv > /05_Imputation/INFO.tab.csv 
awk '{print $1}' /05_Imputation/INFO.tab.csv > /05_Imputation/DR2.csv 
sed "s/DR2=//" /05_Imputation/DR2.csv > /05_Imputation/DR2.txt 
paste /05_Imputation/DR2.txt /05_Imputation/All.genotype.snp.hg37.imputed.noheader.vcf > /05_Imputation/All.genotype.snp.hg37.noheader.DR.vcf 

sed -e '/,/d' /05_Imputation/All.genotype.snp.hg37.noheader.DR.vcf > /05_Imputation/All.genotype.snp.hg37.noheader.DR.singl.vcf

awk '$1>0.7{print $0}' /05_Imputation/All.genotype.snp.hg37.noheader.DR.singl.vcf > /05_Imputation/All.genotype.snp.imputed.vcf

plink 
	--vcf /05_Imputation/All.genotype.snp.imputed.vcf
	--recode 
	--allow-no-sex 
	--out /05_Imputation/All.genotype.snp.imputed
plink 
	--file /05_Imputation/All.genotype.snp.imputed 
	--maf 0.05 
	--recode 
	--out /06_Association/All.genotype.snp.imputed.0.05
```

# 8.Association

------

```shell
/bigdata/gaojy/biosoft/gemma/gemma 
	-bfile /06_Association/All.genotype.snp.imputed.0.05
	-gk 2 
	-o /06_Association/All.genotype.snp.imputed.0.05

/bigdata/gaojy/biosoft/gemma/gemma 
	-bfile /06_Association/All.genotype.snp.imputed.0.05
	-k /06_Association/All.genotype.snp.imputed.0.05.sXX.txt 
	-lmm 4 
	-c /06_Association/co.txt
	-o /06_Association/All.genotype.snp.imputed.0.05 
```

