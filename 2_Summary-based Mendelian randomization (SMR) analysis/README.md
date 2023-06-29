# Summary-based Mendelian randomization (SMR) analysis

------

```shell
#eQTL
zcat /06_Association/All.genotype.snp.imputed.0.05.assoc.txt.gz | awk '{print $2,$5,$6,$7,$8,$9,$13,"NA"}' > /07_SMR/All.genotype.snp.imputed.0.05.ma
sed -i '1i\SNP A1 A2 freq b se p n' /07_SMR/All.genotype.snp.imputed.0.05.ma
sed -i '2d' /07_SMR/All.genotype.snp.imputed.0.05.ma

/data/gaojy/biosoft/smr-1.3.1-linux-x86_64/smr-1.3.1 
--bfile /06_Association/All.genotype.snp.imputed.0.05
--gwas-summary /07_SMR/All.genotype.snp.imputed.0.05.ma 
--beqtl-summary /data/gaojy/biosoft/smr-1.3.1-linux-x86_64/Lung.lite.eQTL/Lung.lite 
--thread-num 20 
--diff-freq-prop 0.05 
--out /07_SMR/All.genotype.snp.imputed.0.05.smr

awk -F'\t' '{if($19<0.005&&$20>0.05&&$21>10) print $0}' /07_SMR/All.genotype.snp.imputed.0.05.smr > /07_SMR/All.genotype.snp.imputed.0.05.PASS.smr

#sQTL
for chr in {1..22}; 
do 
	/data/gaojy/biosoft/smr-1.3.1-linux-x86_64/smr-1.3.1 
	--bfile /06_Association/All.genotype.snp.imputed.0.05
	--gwas-summary /07_SMR/All.genotype.snp.imputed.0.05.ma 
	--beqtl-summary /data/gaojy/biosoft/smr-1.3.1-linux-x86_64/Lung.lite.sQTL/sQTL_Lung.lite.chr$chr 
	--out ./plot/ENSG00000116288.12
	--plot
	--probe chr1:7961735:7962763:clu_55378:ENSG00000116288.12 
	--probe-wind 500 
	--gene-list /data/gaojy/biosoft/smr-1.3.1-linux-x86_64/glist-hg19; 
done
```

