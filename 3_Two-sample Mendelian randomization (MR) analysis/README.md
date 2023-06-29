# Two-sample Mendelian randomization (MR) analysis

------

```R
library(TwoSampleMR)
copd <- extract_instruments(outcomes='bbj-a-103', clump=TRUE, r2=0.02, kb=5000, access_token=NULL) #获取暴露数据
head(copd) 
dim(copd)

#Calculate PVE
submicro$pve = (2*(submicro$beta^2*submicro$af*(1-submicro$af)))/(2*submicro$beta*submicro$af*(1-submicro$af) + submicro$se^2*2*N*submicro$af*(1-submicro$af))
write.table(submicro,"D:/R_data/mr/taxa/pve.txt",sep="\t",quote= F,col.names = NA)

#Calculate F value
pve = read.table('D:/R_data/mr/taxa/pve.txt', header=1)
pve$F = pve$pve*(97-2)/(1-pve$pve)
write.table(pve,"D:/R_data/mr/taxa/F.txt",sep="\t",quote= F,col.names = NA)


t2d = read.table('D:/R_data/mr/Outcome_GWAS_summary.txt.gz', header=1)
t2d$phenotype <- 'outcome' 
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

mydata <- harmonise_data(exposure_dat=copd, outcome_dat=t2d_out, action= 2)
res <- mr(mydata)
write.table(res,"D:/R_data/mr.txt",sep="\t",quote= F,col.names = NA)

#reverse MR
microbiome <-read.table('D:/R_data/mr/Exposure_GWAS_summary.txt.gz', header=T)#读入暴露文件
submicro <- subset(microbiome,microbiome$p_wald<1e-5)
write.table(submicro,"D:/R_data/Exposure_GWAS_summary.sub.txt",sep="\t",quote= F,col.names = NA)

exp_dat <- read_exposure_data(filename = 'D:/R_data/Exposure_GWAS_summary.sub.txt',
                              clump = FALSE,
                              sep= "\t",
                              snp_col = "rs",
                              beta_col = "beta",
                              se_col = "se",
                              eaf_col = "af",
                              effect_allele_col = "allele1",
                              other_allele_col = "allele0",
                              pval_col = "p_wald")

submicro_clump <-clump_data(exp_dat,clump_r2=0.02, clump_kb=5000, pop='EAS')
t2d_out <- extract_outcome_data(snps=submicro_clump$SNP, outcomes='bbj-a-103', 
                                proxies = FALSE, 
                                maf_threshold = 0.001, 
                                access_token = NULL)
mydata <- harmonise_data(exposure_dat=submicro_clump, outcome_dat=t2d_out, action= 2)
res <- mr(mydata)

#heterogeneity test
het <- mr_heterogeneity(mydata)
mr(mydata,method_list=c('mr_ivw_fe'))

#pleiotropy test
pleio <- mr_pleiotropy_test(mydata)

single <- mr_leaveoneout(mydata)
mr_leaveoneout_plot(single)

#plot
p1 <- mr_scatter_plot(res, mydata)

methods = all_method=c('mr_ivw', 'mr_egger_regression', 'mr_weighted_mode', 'mr_weighted_median')
res_single <- mr_singlesnp(mydata, all_method=methods)
p2 <- mr_forest_plot(res_single)

res_single <- mr_singlesnp(mydata)
p3 <- mr_funnel_plot(res_single)

```

