#FPKM calculation
library(AnnotationHub)
library(GenomicFeatures)
ah <- AnnotationHub()
ath <- query(ah,'Arabidopsis thaliana')
ath_tx <- ath[['AH52247']]
athTAIR10TXDB <- makeTxDbFromGFF("TAIR10_GFF3_genes_transposons.gff")
columns(athTAIR10TXDB)
k=keys(athTAIR10TXDB,keytype = "GENEID")
df=select(athTAIR10TXDB, keys=k, keytype = "GENEID",columns = "TXNAME")
head(df)
tx2gene=df[,2:1]
head(tx2gene)
gene_txdb <- genes(athTAIR10TXDB)
exon_txdb <- exons(athTAIR10TXDB)
overlapinfo = findOverlaps(exon_txdb, gene_txdb)
overlapinfo
t1 = exon_txdb[queryHits(overlapinfo)]
t2 = gene_txdb[subjectHits(overlapinfo)]
t1$gene_id=mcols(t2)[,1]
t1=as.data.frame(t1)
g_l = lapply(split(t1, t1$gene_id), function(x){
  tmp=apply(x, 1, function(y){
    y[2]:y[3]
  })
  length(unique(unlist(tmp)))
})
g_l = data.frame(gene_id = names(g_l), lenth = as.numeric(g_l))
head(g_l)

pistilcounts <- read.csv("pistil_counts.csv", header = TRUE, row.names=1)
pistilcounts[1:6, 1:12]
ng=intersect(rownames(pistilcounts), g_l$gene_id)
exprSet=pistilcounts[ng,]
lengths=g_l[match(ng, g_l$gene_id), 2]
head(lengths)
head(rownames(exprSet))
exprSet[1:6, 1:12]
dim(pistilcounts)
dim(exprSet)
total_count <- colSums(exprSet)
total_count[1:12]
rpkm <- t(do.call(rbind, 
                  lapply(1:length(total_count),
                         function(i){
                           10^9*exprSet[,i]/lengths/total_count[i]
                         })))
rpkm[1:6, 1:12]
rownames(rpkm)=rownames(exprSet)
colnames(rpkm)=colnames(exprSet)
rpkm[1:6, 1:12]
write.table(rpkm,"pistil_rpkmNormalized-TAIR10.csv", sep=",")

#Differential expression analysis
library(DESeq2)
raw.count=read.csv("Counts_WT22_tcp22.csv")
count.data=raw.count[,2:7]
row.names(count.data)=raw.count[,1]
condition=factor(c("WT_22","WT_22","WT_22","tcp_22","tcp_22","tcp_22"),levels = c("WT_22","tcp_22"))
col_data=data.frame(row.names = colnames(count.data),condition)
dds=DESeqDataSetFromMatrix(countData = count.data,colData = col_data,design = ~ condition)
dds_out=DESeq(dds)
res=results(dds_out)
write.csv(res,file = "FC_WT22_tcp22.csv")
rm(list=ls()) 
---------
raw.count=read.csv("Counts_WT28_tcp28.csv")
count.data=raw.count[,2:7]
row.names(count.data)=raw.count[,1]
condition=factor(c("WT_28","WT_28","WT_28","tcp_28","tcp_28","tcp_28"),levels = c("WT_28","tcp_28")) 
col_data=data.frame(row.names = colnames(count.data),condition)
dds=DESeqDataSetFromMatrix(countData = count.data,colData = col_data,design = ~ condition)
dds_out=DESeq(dds)
res=results(dds_out)
write.csv(res,file = "FC_WT28_tcp28.csv")
rm(list=ls()) 
----------
raw.count=read.csv("Counts_WT22_WT28.csv")
count.data=raw.count[,2:7]
row.names(count.data)=raw.count[,1]
condition=factor(c("WT_22","WT_22","WT_22","WT_28","WT_28","WT_28"),levels = c("WT_22","WT_28")) 
col_data=data.frame(row.names = colnames(count.data),condition)
dds=DESeqDataSetFromMatrix(countData = count.data,colData = col_data,design = ~ condition)
dds_out=DESeq(dds)
res=results(dds_out)
write.csv(res,file = "FC_WT22_WT28.csv")
rm(list=ls()) 
-------
raw.count=read.csv("Counts_tcp22_tcp28.csv")
count.data=raw.count[,2:7]
row.names(count.data)=raw.count[,1]
condition=factor(c("tcp_22","tcp_22","tcp_22","tcp_28","tcp_28","tcp_28"),levels = c("tcp_22","tcp_28")) 
col_data=data.frame(row.names = colnames(count.data),condition)
dds=DESeqDataSetFromMatrix(countData = count.data,colData = col_data,design = ~ condition)
dds_out=DESeq(dds)
res=results(dds_out)
write.csv(res,file = "FC_tcp22_tcp28.csv")
------
raw.count=read.csv("Counts_WT22_tcp28.csv")
count.data=raw.count[,2:7]
row.names(count.data)=raw.count[,1]
condition=factor(c("WT_22","WT_22","WT_22","tcp_28","tcp_28","tcp_28"),levels = c("WT_22","tcp_28")) 
col_data=data.frame(row.names = colnames(count.data),condition)
dds=DESeqDataSetFromMatrix(countData = count.data,colData = col_data,design = ~ condition)
dds_out=DESeq(dds)
res=results(dds_out)
write.csv(res,file = "FC_WT22_tcp28.csv")

#volcano plot
  
library("ggpubr")
library("ggthemes")
library("ggplot2")
library("patchwork")
setwd("D:/rna-seq-project/210705_pistil/noneTE")
WT22_tcp22 <- read.csv("FC_WT22_tcp22.csv", header = T)
head(WT22_tcp22)
WT22_tcp22$Group = "not-significant"
WT22_tcp22$Group[which((WT22_tcp22$padj < 0.01) & (WT22_tcp22$log2FoldChange > 1))] = "up-regulated"
WT22_tcp22$Group[which((WT22_tcp22$padj < 0.01) & (WT22_tcp22$log2FoldChange < -1))] = "down-regulated"
WT22_tcp22$log10Padj <- -log10(WT22_tcp22$padj)
WT22_tcp22_volcano <- ggscatter(WT22_tcp22, x = "log2FoldChange", y = "log10Padj", color = "Group", 
                palette = c("#2f5688", "#BBBBBB", "#CC0000"), size = 1, 
                xlab = "log2FoldChange", ylab = "-log10Padj", xlim = c(-10, 10)) +
  theme_base()+
  theme(axis.title = element_text(size = 18)) +
  geom_hline(yintercept = 2.00, linetype = "dashed") +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed") 
WT22_tcp22_volcano
--
WT28_tcp28 <- read.csv("FC_WT28_tcp28.csv", header = T)
head(WT28_tcp28)
WT28_tcp28$Group = "not-significant"
WT28_tcp28$Group[which((WT28_tcp28$padj < 0.01) & (WT28_tcp28$log2FoldChange > 1))] = "up-regulated"
WT28_tcp28$Group[which((WT28_tcp28$padj < 0.01) & (WT28_tcp28$log2FoldChange < -1))] = "down-regulated"
WT28_tcp28$log10Padj <- -log10(WT28_tcp28$padj)
WT28_tcp28_volcano <- ggscatter(WT28_tcp28, x = "log2FoldChange", y = "log10Padj", color = "Group", 
                palette = c("#2f5688", "#BBBBBB", "#CC0000"), size = 1, 
                xlab = "log2FoldChange", ylab = "-log10Padj", xlim = c(-10, 10)) +
  theme_base()+
  theme(axis.title = element_text(size = 18)) +
  geom_hline(yintercept = 2.00, linetype = "dashed") +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed") 
WT28_tcp28_volcano
--
WT22_WT28 <- read.csv("FC_WT22_WT28.csv", header = T)
head(WT22_WT28)
WT22_WT28$Group = "not-significant"
WT22_WT28$Group[which((WT22_WT28$padj < 0.01) & (WT22_WT28$log2FoldChange > 1))] = "up-regulated"
WT22_WT28$Group[which((WT22_WT28$padj < 0.01) & (WT22_WT28$log2FoldChange < -1))] = "down-regulated"
WT22_WT28$log10Padj <- -log10(WT22_WT28$padj)
WT22_WT28_volcano <- ggscatter(WT22_WT28, x = "log2FoldChange", y = "log10Padj", color = "Group", 
                                palette = c("#2f5688", "#BBBBBB", "#CC0000"), size = 1, 
                                xlab = "log2FoldChange", ylab = "-log10Padj", xlim = c(-10, 10)) +
  theme_base()+
  theme(axis.title = element_text(size = 18)) +
  geom_hline(yintercept = 2.00, linetype = "dashed") +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed") 
WT22_WT28_volcano
--
tcp22_tcp28 <- read.csv("FC_tcp22_tcp28.csv", header = T)
head(tcp22_tcp28)
tcp22_tcp28$Group = "not-significant"
tcp22_tcp28$Group[which((tcp22_tcp28$padj < 0.01) & (tcp22_tcp28$log2FoldChange > 1))] = "up-regulated"
tcp22_tcp28$Group[which((tcp22_tcp28$padj < 0.01) & (tcp22_tcp28$log2FoldChange < -1))] = "down-regulated"
tcp22_tcp28$log10Padj <- -log10(tcp22_tcp28$padj)
tcp22_tcp28_volcano <- ggscatter(tcp22_tcp28, x = "log2FoldChange", y = "log10Padj", color = "Group", 
                               palette = c("#2f5688", "#BBBBBB", "#CC0000"), size = 1, 
                               xlab = "log2FoldChange", ylab = "-log10Padj", xlim = c(-10, 10)) +
  theme_base()+
  theme(axis.title = element_text(size = 18)) +
  geom_hline(yintercept = 2.00, linetype = "dashed") +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed") 
tcp22_tcp28_volcano
WT22_tcp22_volcano + WT28_tcp28_volcano + WT22_WT28_volcano + tcp22_tcp28_volcano


#tcp_related gene co-expression cluster analysis
setwd("D:/rna-seq-project/210705_pistil/noneTE/tcp_related")
zscore <- read.csv("Geneid_tcprelated_FC1Padj001_zscore_mean_up.csv",header=TRUE,row.names=1)

tested_cluster <- 20
wss <- (nrow(zscore)-1) * sum(apply(zscore, 2, var))
for (i in 2:tested_cluster) {
  wss[i] <- kmeans(zscore, centers=i,iter.max=100,  nstart=25)$tot.withinss
}
plot(1:tested_cluster, wss, type="b", xlab="Number of Clusters", 
     ylab="Within groups sum of squares")
---
down_zscore <- read.csv("Geneid_tcprelated_FC1Padj001_zscore_mean_down.csv",header=TRUE,row.names=1)

tested_cluster <- 20
wss <- (nrow(down_zscore)-1) * sum(apply(down_zscore, 2, var))
for (i in 2:tested_cluster) {
  wss[i] <- kmeans(down_zscore, centers=i,iter.max=100,  nstart=25)$tot.withinss
}
plot(1:tested_cluster, wss, type="b", xlab="Number of Clusters", 
     ylab="Within groups sum of squares")

library(cluster)
library(fpc)

center = 3
fit <- kmeans(zscore, centers=center, iter.max=100, nstart=25)
withinss <- fit$tot.withinss
print(paste("Get withinss for the first run", withinss))
try_count = 10
for (i in 1:try_count) {
  tmpfit <- kmeans(zscore, centers=center, iter.max=100, nstart=25)
  tmpwithinss <- tmpfit$tot.withinss
  print(paste(("The additional "), i, 'run, withinss', tmpwithinss))
  if (tmpwithinss < withinss){
    withins <- tmpwithinss
    fit <- tmpfit
  }
}
fit_cluster = fit$cluster
clusplot(zscore, fit_cluster, shade=T, labels=5, lines=0, color=T,
         lty=4, main='K-means clusters')

cluster_mean <- aggregate(zscore, by=list(fit_cluster), FUN=mean)
write.table(t(cluster_mean), file="ehbio.pam.cluster.mean.xls", sep='\t',
            col.names=F, row.names=T, quote=F)

write.table(fit_cluster, file="fit_cluster_3.xls", sep='\t',
            col.names=F, row.names=T, quote=F)
----
  center = 3
fit <- kmeans(down_zscore, centers=center, iter.max=100, nstart=25)
withinss <- fit$tot.withinss
print(paste("Get withinss for the first run", withinss))
try_count = 10
for (i in 1:try_count) {
  tmpfit <- kmeans(down_zscore, centers=center, iter.max=100, nstart=25)
  tmpwithinss <- tmpfit$tot.withinss
  print(paste(("The additional "), i, 'run, withinss', tmpwithinss))
  if (tmpwithinss < withinss){
    withins <- tmpwithinss
    fit <- tmpfit
  }
}
fit_cluster = fit$cluster
clusplot(zscore, fit_cluster, shade=T, labels=5, lines=0, color=T,
         lty=4, main='K-means clusters')

write.table(fit_cluster, file="fit_cluster_3.xls", sep='\t',
            col.names=F, row.names=T, quote=F)


setwd("D:/rna-seq-project/210705_pistil/noneTE/tcp_related")
library(reshape2)
library(ggplot2)
up_zscore_mean <- read.csv("Geneid_tcprelated_FC1Padj001_zscore_mean_up.csv",header=TRUE,row.names=1)
ggparcoord(up_zscore_mean, columns = 1:4,showPoints = TRUE,groupColumn=5, scale = "globalminmax",
           title = "", alphaLines = 0.3)+
 facet_wrap(~group) +
  theme_bw() +
  theme(plot.title = element_text(size=10))

down_zscore_mean <- read.csv("Geneid_tcprelated_FC1Padj001_zscore_mean_down.csv",header=TRUE,row.names=1)
ggparcoord(down_zscore_mean, columns = 1:4,showPoints = TRUE,groupColumn=5, scale = "globalminmax",
           title = "", alphaLines = 0.3)+
  facet_wrap(~group) +
  theme_bw() +
  scale_color_brewer(palette = "Blues", direction = -1) +
  theme(plot.title = element_text(size=10))


#TCP-heatmap
setwd("D:/rna-seq-project/210705_pistil/noneTE/FPKM")
library(pheatmap)
TCP <- read.csv(file = "FPKM_TCP.csv", header = T, row.names = 1)
pheatmap(TCP, cluster_row = FALSE, cluster_cols = FALSE, display_numbers = TRUE, cellwidth = 40, cellheight = 15, scale = "row", color = colorRampPalette(c("#2f5688", "#FFFFFF", "#CC0000"))(20))


##GO bubble
library("ggplot2")
library("forcats")
library("ggsci")
library(RColorBrewer)
setwd("D:/rna-seq-project/210705_pistil/noneTE/tcp_related")
GO3 <- read.csv("cluster3_GO.csv")
p1 <- ggplot(GO3, aes(x=Enrichment,y=Term,size=lgFDR))  + 
  geom_point(color ="#DC0000")+
  scale_size(rang=c(5,15))+
  theme_bw(base_size = 12)+
  theme(text = element_text(size = 15))
p1
GO4 <- read.csv("cluster4_GO.csv")
p2 <- ggplot(GO4, aes(x=Enrichment,y=Term,size=lgFDR)) + 
  geom_point(color = "#3c5488")+
  xlim(5, 20)+
  scale_size(rang=c(5,15))+
  theme_bw(base_size = 12)+
  theme(text = element_text(size = 15))
p2

##venn plot of tcpII-related-upregulated genes
setwd("D:/rna-seq-project/210705_pistil/noneTE/tcp_related")
tcprelated <- read.csv("Geneid_tcprelated_FC1Padj001.csv")
tcp_inter <- get.venn.partitions(tcprelated)
for (i in 1:nrow(tcp_inter)) tcp_inter[i,'values'] <- paste(tcp_inter[[i,'..values..']], collapse = ', ')
write.table(tcp_inter[-c(5, 6)], 'tcp_inter.txt', row.names = FALSE, sep = '\t', quote = FALSE)

tcprelated_venn <- venn.diagram(
  x = list(WT22_tcp22_up  = tcprelated$WT22_tcp22_up , 
           WT28_tcp28_up  = tcprelated$WT28_tcp28_up ), 
  filename = "tcprelated_up_venn.png", 
  imagetype = "png", lwd = 3, fill = c("#F39B7F", "#DC0000"), alpha = 0.6,  
  label.col = "white", cex =1.5, fontfamily = "sans", fontface = "bold",  
  cat.col = c("#F39B7F", "#DC0000"),  cat.cex = 0.7, cat.fontfamily = "sans", 
  cat.fontface = "bold", margin = 0.05)

##venn plot of tcpII-related-downregulated genes
tcprelated_venn <- venn.diagram(
  x = list(WT22_tcp22_down   = tcprelated$WT22_tcp22_down , 
           WT28_tcp28_down  = tcprelated$WT28_tcp28_down ), 
  filename = "tcprelated_down_venn.png", 
  imagetype = "png", lwd = 3, fill = c("#3C5488", "#4DBBD5"), alpha = 0.6,  
  label.col = "white", cex =1.5, fontfamily = "sans", fontface = "bold",  
  cat.col = c("#3C5488", "#4DBBD5"),  cat.cex = 0.7, cat.fontfamily = "sans", 
  cat.fontface = "bold", margin = 0.05)
  
  