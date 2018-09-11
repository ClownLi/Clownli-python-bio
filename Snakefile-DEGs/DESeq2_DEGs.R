#!/usr/bin/Rscript
##############################################################################################
#Author:limingzhu
#Date:2018-07-22
#DESeq2 find DEGs
##############################################################################################
WorkDir = getwd()
SetDir = paste(WorkDir,"Results",sep = '/')
dir.create(SetDir)

rm(list=ls())
args<-commandArgs(T)


library("DESeq2")
raw_count <- read.table(args[1],sep = "\t",header = TRUE) 
count_data <- raw_count[,2:7]
row.names(count_data) <- raw_count[,1]
condition <- factor(c("control","control","control","treat","treat","treat"))
col_data <- data.frame(row.names = colnames(count_data), condition) 

##############################################################################################

dds <- DESeqDataSetFromMatrix(countData = count_data,colData = col_data,design = ~ condition)
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]
dds_out <- DESeq(dds)
res <- results(dds_out,alpha = 0.05,na.rm = TRUE)
summary(res)

##############################################################################################
#boxplot
exprSet_new=assay(dds_out)
par(cex = 0.7)
n.sample=ncol(count_data)
cols <- rainbow(n.sample*1.2)
pdf(file="Results/expression_boxplot.pdf")
boxplot(log2(exprSet_new+1), col = cols,main="expression value",las=2)
dev.off()

##############################################################################################

library( "gplots" )
library( "RColorBrewer" )
hmcol<- colorRampPalette(brewer.pal(9, 'GnBu'))(100)
rld<- rlogTransformation(dds, blind=TRUE)
vsd<-varianceStabilizingTransformation(dds, blind=TRUE)
distsRL <- dist(t(assay(rld)))
mat<- as.matrix(distsRL)
sampleFiles <- as.character(c("1","2","3","1","2","3"))
rownames(mat) <- colnames(mat) <- with(colData(dds),
                                       paste(condition,sampleFiles , sep='-'))
##############################################################################################
#heatmaps_samples
hc <- hclust(distsRL)
pdf(file="Results/heatmaps.pdf")
heatmap.2(mat, Rowv=as.dendrogram(hc),
          symm=TRUE, trace='none',
          col = rev(hmcol), margin=c(13, 13))
dev.off()

#PCA
pdf(file="Results/PCA.pdf")
plotPCA(rld, intgroup=c('condition'))
dev.off()

##############################################################################################

resLFC <- lfcShrink(dds=dds_out,coef=2,res=res)
write.table(resLFC,"Results/tbb.deseq2.resLFC.csv",quote = F,sep = ",")

table(res$padj<0.05) 
res_deseq <- res[order(res$padj),] 
diff_gene_deseq2 <- subset(res_deseq, padj<0.05 & (log2FoldChange > 1 | log2FoldChange < -1)) 

diff_gene_deseq2 <- row.names(diff_gene_deseq2) 
res_diff_data <- merge(as.data.frame(res),as.data.frame(counts(dds_out,normalize=TRUE)),by="row.names",sort=FALSE) 
write.csv(res_diff_data,file = "Results/Control_Treat.csv",row.names = F)

##############################################################################################
#MA
library(ggplot2)
pdf(file="Results/MA1.pdf")
plotMA(res,ylim=c(-5,5)) 
topGene <- rownames(res)[which.min(res$padj)] 
with(res[topGene, ], { 
  points(baseMean, log2FoldChange, col="dodgerblue", cex=2, lwd=2) 
  text(baseMean, log2FoldChange, topGene, pos=2, col="dodgerblue") 
})
dev.off()

pdf(file="Results/MA2.pdf")
resLFC <- lfcShrink(dds_out,coef = 2,res=res) 
plotMA(resLFC, ylim=c(-5,5)) 
topGene <- rownames(resLFC)[which.min(res$padj)] 
with(resLFC[topGene, ], { 
  points(baseMean, log2FoldChange, col="dodgerblue", cex=2, lwd=2) 
  text(baseMean, log2FoldChange, topGene, pos=2, col="dodgerblue") 
}) 
idx <- identify(res$baseMean, res$log2FoldChange)
dev.off()
##############################################################################################
#heatmap
library(genefilter) 
library(pheatmap) 

rld <- rlogTransformation(dds_out,blind = F) 
write.csv(assay(rld),file="Results/mm.DESeq2.pseudo.counts.csv") 
topVarGene <- head(order(rowVars(assay(rld)),decreasing = TRUE),1000) 
mat <- assay(rld)[ topVarGene, ] 

anno <- as.data.frame(colData(rld)[,c("condition","sizeFactor")])
pdf(file="Results/heatmap.pdf",height=12,width=8)
pheatmap(mat, annotation_col = anno,show_rownames = F)
dev.off()

#################################################################################
#up
resOrdered <- resLFC[order(res$padj),]
resOrdered=as.data.frame(resOrdered)
head(resOrdered)
filter_vector_p_value = resOrdered$padj < 0.05
filter_vector_up = resOrdered$log2FoldChange > 1
filter_vector_up_all = filter_vector_p_value & filter_vector_up

res.sign.up = resOrdered[filter_vector_up_all,]
write.table(res.sign.up,"Results/up.csv",quote = F,sep = ",")

#down
filter_vector_down = resOrdered$log2FoldChange < -1
filter_vector_down_all = filter_vector_p_value & filter_vector_down
res.sign.down = resOrdered[filter_vector_down_all,]
write.table(res.sign.down,"Results/down.csv",quote = F,sep = ",")
################################################################################
#Volcano plot
pdf(file="Results/Volcano_plot.pdf")
par(mar=c(5,5,5,5), cex=1.0, cex.main=1.4, cex.axis=1.4, cex.lab=1.4)

topT <- as.data.frame(res)
with(topT, plot(log2FoldChange, -log10(padj), pch=20, main="Volcano plot", cex=1.0, xlab=bquote(~Log[2]~fold~change), ylab=bquote(~-log[10]~Q~value)))

with(subset(topT, padj<0.05 & log2FoldChange > 2), points(log2FoldChange, -log10(padj), pch=20, col="red", cex=0.5))

with(subset(topT, padj<0.05 & log2FoldChange < -2), points(log2FoldChange, -log10(padj), pch=20, col="green", cex=0.5))

abline(v=0, col="black", lty=3, lwd=1.0)
abline(v=-2, col="black", lty=4, lwd=2.0)
abline(v=2, col="black", lty=4, lwd=2.0)
abline(h=-log10(max(topT$pvalue[topT$padj<0.05], na.rm=TRUE)), col="black", lty=4, lwd=2.0)
dev.off()
