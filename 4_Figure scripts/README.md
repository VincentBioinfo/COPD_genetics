# Figure scripts

------

## PERMANOVA

```R

library("vegan")
dis <- read.delim('D:/R_data/bray.curtis.txt',header=T,sep="\t",comment.char="",check.names=F)
dis <- as.dist(dis)
group <- read.delim('D:/R_data/Adonis.group.txt', sep = '\t', stringsAsFactors = FALSE)
adonis_result_dis <- adonis2(dis~GeneticPC1, group, permutations = 999) 
write.table(adonis_result_dis,'D:/R_data/Adonis.GeneticPC1.txt',sep='\t')

library(ggplot2)
library(dplyr)
library(data.table)
setwd("D:/R_data/adonis/")
list.files()
adonis_df <-  fread("adonis.csv") %>% arrange(adonis.R)
adonis_df <- adonis_df %>% mutate(color = sapply(adonis.p,function(x) ifelse(x < 0.05,"Sig","None")))
Colors_p <- c("#AD1E25","#BEBDBD")
names(Colors_p) <- c("Sig","None")
adonis_df$color <- factor(adonis_df$color,levels = c("Sig","None"))
adonis_df$clinical_Vars <- factor(adonis_df$clinical_Vars,levels =adonis_df$clinical_Vars)
gg <- ggplot(adonis_df,aes(x=clinical_Vars,y=adonis.R,fill =color))+geom_bar(stat="identity",width = 0.8) +
  scale_y_continuous(breaks=seq(0, 0.3, 0.05)) +
  coord_flip() +
  geom_text(aes(label=clinical_Vars), hjust = -0.1) +
  scale_fill_manual(values = Colors_p)+theme_light() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank()) 

ggsave(gg,filename = "adonis.COPD.pdf",width = 10, height = 7, units = "in", dpi = 300)
```

## Spearman correlation

```R
library('Hmisc')
data=read.table("D:/spearman.corr/sample.txt",
                header = 1,
                sep="\t",
                row.names = 1,
                comment.char = "",
                check.names = F,
                na.strings = "NA") 
res <- rcorr(as.matrix(data),type="spearman")
r <- res$r[1,]
write.table(r,'D:/spearman.corr/out/output.txt',
            sep='\t',
            append=T,
            quote=F,
            col.names=T)
p <- res$P[1,]
write.table(p,'D:/spearman.corr/out/output.txt',
            sep='\t',
            append=T,
            quote=F,
            col.names=T)

library('dplyr')
library('ggplot2')
library('reshape2')
library('ggpubr')
library('ggpmisc')
library('ggstatsplot')
library('ggcharts')

setwd("'D:/spearman.corr/")
data<-read.delim('sample.txt',header = T,sep = "\t",row.names = 1,comment.char = "",check.names = F)
ggscatterstats(data,
               x = column1,
               y = column2,
               xlab = "column1",
               ylab = "column2",
               type = "spearman",
               ggtheme = theme_bw(),
               ggstatsplot.layer = FALSE,
               line.color = "red",
               centrality.para = "mean",
               xfill = "#009E73",
               yfill = "#D55E00",
               marginal.type = "density",
               messages = FALSE,
               title = "Relationship between Host Genetic and Airway Microbiome Feature")
ggsave(filename = "correlation.pdf",width = 7, height = 7, units = "in", dpi = 300)
```













