#加载包
library(limma)
#导入数据
mRNA_FPKM_count<- read.table("GSE140275_mRNA_FPKM.txt",header = T, stringsAsFactors = F,row.names = 1) 
dat1<- data.frame(mRNA_FPKM_count[4:9])

#查看差异化的箱线图
boxplot(dat1)
#dat4 <- log(dat1 + 1) 如果差异很大取log归一化

#2.将FPKM转换为TPM
expMatrix <- dat1
fpkmToTpm <- function(fpkm)
{
  exp(log(fpkm) - log(sum(fpkm)) + log(1e6))
}
tpms <- apply(expMatrix,2,fpkmToTpm)
tpms[1:3,]
colSums(tpms)
tpms[1:3,]

#3.差异分析

group_list=c(rep('ctr',3),rep('stroke',3))
## 强制限定顺序
group_list <- factor(group_list,levels = c("ctr","stroke"),ordered = F)

#表达矩阵数据校正
exprSet <- tpms
boxplot(exprSet,outline=FALSE, notch=T,col=group_list, las=2)
library(limma) 
exprSet=normalizeBetweenArrays(exprSet)
boxplot(exprSet,outline=FALSE, notch=T,col=group_list, las=2)
#判断数据是否需要转换
exprSet <- log2(exprSet+1)
#差异分析：
dat <- exprSet
design=model.matrix(~factor( group_list ))
fit=lmFit(dat,design)
fit=eBayes(fit)
options(digits = 4)
topTable(fit,coef=2,adjust='BH')
bp=function(g){
  library(ggpubr)
  df=data.frame(gene=g,stage=group_list)
  p <- ggboxplot(df, x = "stage", y = "gene",
                 color = "stage", palette = "jco",
                 add = "jitter")
  #  Add p-value
  p + stat_compare_means()
}
deg=topTable(fit,coef=2,adjust='BH',number = Inf)
head(deg) 
#save(deg,file = 'deg.Rdata')

#标记基因上调或者下调
deg[which(deg$logFC >= 1 & deg$adj.P.Val < 0.05),'direction'] <- 'up'
deg[which(deg$logFC <= -1 & deg$adj.P.Val  < 0.05),'direction'] <- 'down'
deg[which(abs(deg$logFC) <= 1 | deg$adj.P.Val  >= 0.05),'direction'] <- 'none'

colnames(deg)[1]<- "log2FoldChange" 
colnames(deg)[5]<- "padj" 

#画火山图看差异基因
library(tidyverse)
deg<-rownames_to_column(deg,var = "ID")
top_de <- filter(deg,abs(log2FoldChange) > 5)
#
library(ggrepel)
my_palette <- c('#4DBBD5FF','#999999','#E64B35FF')
#创建一个画布
library(ggplot2)
#添加几何对象  geom_point散点图，将direction映射给点颜色  aes映射颜色
ggplot(data = deg,aes(x = log2FoldChange, y = -log10(padj))) +
  geom_point(aes(color = direction,size = abs(log2FoldChange))) + 
  geom_label_repel(data = top_de,aes(label = ID))+
  geom_hline(yintercept = -log10(0.05),linetype = 'dashed') +
  geom_vline(xintercept = c(-1,1),linetype = 'dashed') +
  scale_color_manual(values = my_palette) +
  scale_size(range = c(0.1,2)) +
  labs(x = 'log2 fold change',
       y = '-log10(padj)',
       title = 'Volcano plot',
       size = 'log2 fold change') +
  guides(size = FALSE) +
  theme_classic() + 
  theme(plot.title = element_text(size = 18,hjust = 0.5),
        legend.background = element_blank(),
        legend.position = c(0.93,0.85)) 

