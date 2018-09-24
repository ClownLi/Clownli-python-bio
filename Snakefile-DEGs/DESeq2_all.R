##############################################################################################
#Author:limingzhu
#Date:2018-07-22
#DESeq2 find DEGs
##############################################################################################
rm(list=ls())

#input data
setwd("～")
library("DESeq2")
raw_count <- read.table("result.txt",sep = "\t",header = TRUE) 
count_data <- raw_count[,2:7]
row.names(count_data) <- raw_count[,1]
condition <- factor(c("control","control","control","treat","treat","treat"))
col_data <- data.frame(row.names = colnames(count_data), condition) 

##############################################################################################
####构造dds的对象,统计差异基因

#构建dds矩阵的基本代码
dds <- DESeqDataSetFromMatrix(countData = count_data,colData = col_data,design = ~ condition)
#过滤低质量的低count数据
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]
#使用DESeq进行差异表达分析
dds_out <- DESeq(dds)
res <- results(dds_out,alpha = 0.05,na.rm = TRUE)
summary(res)

##############################################################################################
#表达量箱线图
exprSet_new=assay(dds_out)
par(cex = 0.7)
n.sample=ncol(count_data)
cols <- rainbow(n.sample*1.2)
pdf(file="expression_boxplot.pdf")
boxplot(log2(exprSet_new+1), col = cols,main="expression value",las=2)
dev.off()

##############################################################################################
# #clustering of samples图
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
hc <- hclust(distsRL)
pdf(file="heatmaps.pdf")
heatmap.2(mat, Rowv=as.dendrogram(hc),
          symm=TRUE, trace='none',
          col = rev(hmcol), margin=c(13, 13))
dev.off()
#dev.copy(png,'deseq2_heatmaps_samplebysample.png')
#dev.off()

#画PCA图
pdf(file="PCA.pdf")
plotPCA(rld, intgroup=c('condition'))
dev.off()

##############################################################################################
#统计差异基因
resLFC <- lfcShrink(dds=dds_out,coef=2,res=res)
write.table(resLFC,"tbb.deseq2.resLFC.csv",quote = F,sep = ",")

table(res$padj<0.05) 
res_deseq <- res[order(res$padj),] 
diff_gene_deseq2 <- subset(res_deseq, padj<0.05 & (log2FoldChange > 1 | log2FoldChange < -1)) 
## or diff_gene_deseq2 <- subset(res_desq, padj<0.05 & abs(log2FoldChange) >=1) 
diff_gene_deseq2 <- row.names(diff_gene_deseq2) 
res_diff_data <- merge(as.data.frame(res),as.data.frame(counts(dds_out,normalize=TRUE)),by="row.names",sort=FALSE) 
write.csv(res_diff_data,file = "Control_Treat.csv",row.names = F)

##############################################################################################
#绘制MA图
library(ggplot2)
pdf(file="MA1.pdf")
plotMA(res,ylim=c(-5,5)) 
topGene <- rownames(res)[which.min(res$padj)] 
with(res[topGene, ], { 
  points(baseMean, log2FoldChange, col="dodgerblue", cex=2, lwd=2) 
  text(baseMean, log2FoldChange, topGene, pos=2, col="dodgerblue") 
})
dev.off()

pdf(file="MA2.pdf")
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
#绘制热图
library(genefilter) 
library(pheatmap) 
#pdf(file="heatmap2.pdf")
rld <- rlogTransformation(dds_out,blind = F) 
write.csv(assay(rld),file="mm.DESeq2.pseudo.counts.csv") 
topVarGene <- head(order(rowVars(assay(rld)),decreasing = TRUE),1000) 
mat <- assay(rld)[ topVarGene, ] 
### mat <- mat - rowMeans(mat) 减去一个平均值，让数值更加集中。第二个图 
anno <- as.data.frame(colData(rld)[,c("condition","sizeFactor")])
pheatmap(mat, annotation_col = anno,show_rownames = F)
#dev.off()

#################################################################################
#统计up差异基因
resOrdered <- resLFC[order(res$padj),]
resOrdered=as.data.frame(resOrdered)
head(resOrdered)
filter_vector_p_value = resOrdered$padj < 0.05
filter_vector_up = resOrdered$log2FoldChange > 1
filter_vector_up_all = filter_vector_p_value & filter_vector_up
#res = na.omit(res)
res.sign.up = resOrdered[filter_vector_up_all,]
write.table(res.sign.up,"up.csv",quote = F,sep = ",")
#统计down差异基因
filter_vector_down = resOrdered$log2FoldChange < -1
filter_vector_down_all = filter_vector_p_value & filter_vector_down
res.sign.down = resOrdered[filter_vector_down_all,]
write.table(res.sign.down,"down.csv",quote = F,sep = ",")
################################################################################
##火山图
pdf(file="Volcano_plot.pdf")
par(mar=c(5,5,5,5), cex=1.0, cex.main=1.4, cex.axis=1.4, cex.lab=1.4)

topT <- as.data.frame(res)

#Adjusted P values (FDR Q values)
with(topT, plot(log2FoldChange, -log10(padj), pch=20, main="Volcano plot", cex=1.0, xlab=bquote(~Log[2]~fold~change), ylab=bquote(~-log[10]~Q~value)))

with(subset(topT, padj<0.05 & log2FoldChange > 2), points(log2FoldChange, -log10(padj), pch=20, col="red", cex=0.5))

with(subset(topT, padj<0.05 & log2FoldChange < -2), points(log2FoldChange, -log10(padj), pch=20, col="green", cex=0.5))

#with(subset(topT, padj<0.05 & abs(log2FoldChange)< 2), text(log2FoldChange, -log10(padj), labels=subset(rownames(topT), topT$padj<0.05 & abs(topT$log2FoldChange)>2), cex=0.8, pos=3))

#Add lines for absolute FC>2 and P-value cut-off at FDR Q<0.05
abline(v=0, col="black", lty=3, lwd=1.0)
abline(v=-2, col="black", lty=4, lwd=2.0)
abline(v=2, col="black", lty=4, lwd=2.0)
abline(h=-log10(max(topT$pvalue[topT$padj<0.05], na.rm=TRUE)), col="black", lty=4, lwd=2.0)
dev.off()

################################################################################
#GO和KEGG分析

table(res$padj<0.05)
res_deseq <- res[order(res$padj),]
diff_gene_deseq2 <- subset(res_deseq, padj<0.05 & (log2FoldChange > 1 | log2FoldChange < -1))


library(clusterProfiler)
require(DOSE)
library(DO.db)
library(AnnotationHub)
hub <- AnnotationHub()
query(hub, "Aspergillus flavus")
AFL.OrgDb <- hub[["AH59435"]]
#keytypes(AFL.OrgDb)
gene <- row.names(diff_gene_deseq2)
#tansid <- select(sl,keys = gene,columns = c("GENENAME","GO","ENTREZID"),keytype = "SYMBOL")
#head(tansid)

#display_number = c(15, 10, 15)
enrich.go.BP <- enrichGO(
gene = row.names(diff_gene_deseq2),
OrgDb = AFL.OrgDb,
keyType = "SYMBOL",
ont = "BP",
pvalueCutoff = 0.05
)

#Biological process 点图
pdf("Biological_process.pdf")
dotplot(enrich.go.BP, font.size=8)
dev.off()

# pdf(file="~/enrich.go.BP.tree.pdf",width = 10, height = 15)
# plotGOgraph(enrich.go.BP) #画出树型图
# dev.off()

enrich.go.MF <- enrichGO(
gene = row.names(diff_gene_deseq2),
OrgDb = AFL.OrgDb,
keyType = "SYMBOL",
ont = "MF",
pvalueCutoff = 0.05
)
#Biological process 点图
pdf("Molecular_function.pdf")
dotplot(enrich.go.MF, font.size=8)
dev.off()

enrich.go.CC <- enrichGO(
gene = row.names(diff_gene_deseq2),
OrgDb = AFL.OrgDb,
keyType = "SYMBOL",
ont = "CC",
pvalueCutoff = 0.05
)
#Biological process 点图
pdf("Cellular_component.pdf")
dotplot(enrich.go.CC , font.size=8)
dev.off()

ego_result_BP <- as.data.frame(enrich.go.BP)[1:nrow(enrich.go.BP), ]
ego_result_MF <- as.data.frame(enrich.go.MF)[1:nrow(enrich.go.MF), ]
ego_result_CC <- as.data.frame(enrich.go.CC)[1:nrow(enrich.go.CC), ]

go_enrich_df <- data.frame(ID=c(ego_result_BP$ID, ego_result_CC$ID, ego_result_MF$ID),
Description=c(ego_result_BP$Description, ego_result_CC$Description, ego_result_MF$Description),
GeneNumber=c(ego_result_BP$Count, ego_result_CC$Count, ego_result_MF$Count),
type=factor(c(rep("biological process", nrow(enrich.go.BP)), rep("cellular component", nrow(enrich.go.CC)),
rep("molecular function", nrow(enrich.go.MF))), levels=c("molecular function", "cellular component", "biological process")))

## numbers as data on x axis
go_enrich_df$number <- factor(rev(1:nrow(go_enrich_df)))
## shorten the names of GO terms
shorten_names <- function(x, n_word=4, n_char=40){
    if (length(strsplit(x, " ")[[1]]) > n_word || (nchar(x) > 40))
    {
        if (nchar(x) > 40) x <- substr(x, 1, 40)
        x <- paste(paste(strsplit(x, " ")[[1]][1:min(length(strsplit(x," ")[[1]]), n_word)],
        collapse=" "), "...", sep="")
        return(x)
    }
    else
    {
        return(x)
    }
}

labels=(sapply(
levels(go_enrich_df$Description)[as.numeric(go_enrich_df$Description)],
shorten_names))
names(labels) = rev(1:nrow(go_enrich_df))

## colors for bar // green, blue, orange
CPCOLS <- c("#8DA1CB", "#FD8D62", "#66C3A5")
library(ggplot2)
p <- ggplot(data=go_enrich_df, aes(x=number, y=GeneNumber, fill=type)) +
geom_bar(stat="identity", width=0.8) + coord_flip() +
scale_fill_manual(values = CPCOLS) + theme_bw() +
scale_x_discrete(labels=labels) +
xlab("GO term") +
theme(axis.text=element_text(face = "bold", color="gray50")) +
labs(title = "The Most Enriched GO Terms")

pdf("GO_Enrichment.pdf")
p
dev.off()

dotplot(enrich.go.BP, font.size=8)
enrichMap(ego, vertex.label.cex=1.2, layout=igraph::layout.kamada.kawai)

#KEGG分析
KEGG <- rownames(diff_gene_deseq2)
ekk <- enrichKEGG(KEGG, keyType = "kegg",organism = "afv", pvalueCutoff = 0.05, pAdjustMethod = "BH", qvalueCutoff = 0.1)
#summary(as.data.frame(ekk))
pdf("KEGG_enrichment.pdf")
dotplot(ekk, font.size=8)
dev.off()
# enrichMap(ekk, vertex.label.cex=1.2, layout=igraph::layout.kamada.kawai,)
# browseKEGG(ekk,'afv01100')




