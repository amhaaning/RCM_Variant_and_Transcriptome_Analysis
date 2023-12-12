library(edgeR)
library(ggplot2)
library(ggrepel)
library(org.Hs.eg.db)
library(AnnotationDbi)
library(pheatmap)
library(ggtext)

setwd("C:/Users/ahaaning/OneDrive - Indiana University/Manuscripts/RCM_Variants_and_RNAseq/RCM_RNAseq/")

#Read in expression data
data <- read.delim("HS_062_ILMN_204_Ware_mRNAseq18_Dec2017_STAR_featureCounts_edgeR_02202018.txt")

#Remove "CMM" from IDs
colnames(data) <- gsub("\\.CMM","",colnames(data))

#Get raw counts
counts <- data[,c(30:47)]

#Remove ".1" from column names
colnames(counts) <- gsub("\\.1","",colnames(counts))

#Set gene symbols as row names
rownames(counts) <- data$GeneSymbol

#Make data frame with sample group
info <- data.frame(Sample=colnames(counts),check.names = F)
info$Group <- "RCM"
info$Group[info$Sample %in% c("P0197","P0256","P0144","P0142")] <- "Control"

#Order info and counts by group and then sample
info <- info[order(info$Group,info$Sample),]
counts <- counts[,info$Sample]

annot <- data.frame(ENTREZID=data$ENTREZID)

#Make a DGEList from counts
MyDGE <- DGEList(counts,genes=annot)

#Calculate normalization factors
MyDGE <- calcNormFactors(MyDGE)

#Estimate dispersions
MyDGE <- estimateDisp(MyDGE)

mds <- plotMDS(MyDGE, cex=0.5, main="edgeR MDS Plot")
mds_data <- data.frame(Sample=rownames(mds$distance.matrix.squared),Dim1=mds$x,Dim2=mds$y,check.names = F)
mds_data <- merge(mds_data,info,by="Sample")

dev.off()

#Plot Dim 1 and Dim 2 - Points have sample IDs
pdf("mRNA_MDS_leadingFC_EdgeR_with_Sample_Labels.pdf",width=5,height=4.5)
ggplot(mds_data,aes(x=Dim1,y=Dim2,label=Sample,color=Group,group=Group)) + geom_text_repel(size=4,key_glyph = draw_key_blank) + theme_bw() + theme(legend.text=element_text(size=14),axis.text.x = element_text(size=14,color="black"),axis.text.y = element_text(size=12,color="black"),axis.title.x = element_text(size=14,face="bold"),axis.title.y = element_text(size=14,color="black",face="bold"),legend.title = element_blank(),legend.position = c(0.05,0.9),legend.background = element_blank(),legend.key = element_blank()) + xlab("Leading logFC dim 1 (18%)") + ylab("Leading logFC dim 2 (15%)") + theme(legend.text=element_markdown(size=12)) + scale_color_manual(labels = paste("<span style='color:",c("#006100","#0000FF"),"'>",unique(mds_data$Group)[2:1],"</span>"),values = c("#006100","#0000FF"))
dev.off()

pdf("BCV_edgeR.pdf",width=6.5,height=5)
plotBCV(MyDGE, cex=0.4, main="edgeR: Biological coefficient of variation (BCV) vs abundance")
dev.off()

#Create a model matrix
design <- model.matrix(~info$Group)
rownames(design) <- info$Sample

#Fit a generalized log-linear model to count data
fit_glmQL <- glmQLFit(MyDGE, design)

#Conduct genewise statistical tests between groups
qlf <- glmQLFTest(fit_glmQL)
qlf_DE_Genes <- topTags(qlf, n=Inf, p=0.05)
qlf_DE_Genes <- qlf_DE_Genes$table

qlf_DE_Genes_Up <- qlf_DE_Genes[qlf_DE_Genes$logFC > 0,]
qlf_DE_Genes_Down <- qlf_DE_Genes[qlf_DE_Genes$logFC < 0,]

write.table(qlf_DE_Genes[qlf_DE_Genes$FDR < 0.01,],"DE_Genes_FDR_0.01.txt",row.names=T,col.names=T,sep="\t",quote=F)

#Get log2 CPM, raw counts, and edgeR results for all genes
logCPM <- data.frame(GeneSymbol=rownames(MyDGE),cpm(MyDGE,log=TRUE,prior.count=1),check.names = F)
counts <- data.frame(GeneSymbol=rownames(MyDGE),counts,check.names = F)
edgeR_all_genes <- topTags(qlf, n=Inf)
edgeR_all_genes <- data.frame(GeneSymbol=rownames(edgeR_all_genes),edgeR_all_genes$table,check.names = F)

#Get gene info for all genes
gene_info <- data[,c(2,1,3:8)]

#Merge gene info with log2 CPM and raw counts
logCPM <- merge(gene_info, logCPM, by="GeneSymbol")
counts <- merge(gene_info, counts, by="GeneSymbol")
edgeR_all_genes <- merge(gene_info, edgeR_all_genes, by="GeneSymbol")

#Write all results to files
write.table(logCPM,"EdgeR_logCPM_all_genes.txt",sep="\t",row.names=F,col.names=T,quote=F)
write.table(counts,"Raw_counts_all_genes.txt",sep="\t",row.names=F,col.names=T,quote=F)
write.table(edgeR_all_genes,"EdgeR_results_qlf_test_all_genes.txt",sep="\t",row.names=F,col.names=T,quote=F)

is.de <- decideTestsDGE(qlf)

#Make MD plot
pdf("MD_Plot.pdf",height=5,width=5)
plotMD(qlf,status=is.de,main = "Mean-Difference Plot\n(RCM vs Control)")
dev.off()

#Create a volcano plot
volcano_data <- data.frame(GeneSymbol=rownames(qlf),logFC=qlf$table$logFC,Minuslog10pval=-1*log10(qlf$table$PValue),is.de@.Data,check.names = F)

volcano_data$Legend <- "Not Sig"
volcano_data$Legend[volcano_data$`info$GroupRCM`<0] <- "Down"
volcano_data$Legend[volcano_data$`info$GroupRCM`>0] <- "Up"

pdf("Volcano_Plot.pdf",height=4,width=4)
ggplot(volcano_data,aes(x=logFC,y=Minuslog10pval,color=Legend)) + geom_point() + theme_bw() + scale_color_manual(values=c("blue","black","red")) + theme(legend.title = element_blank(),legend.position = c(0.12,0.12),legend.background = element_blank(),legend.key = element_blank()) + xlab("log2(fold-change)") + ylab("-log10(p-value)")
dev.off()

#Get table of DE genes with FDR p-val < 0.05
DE_Genes <- topTags(qlf, n=Inf, p=0.05)
DE_Genes <- DE_Genes$table

#Make a heatmap with all DE genes (FDR adjusted p-val < 0.05)####

#Get log2 CPM
logCPM <- cpm(MyDGE,log=TRUE,prior.count=1)

#Keep DE genes
logCPM <- logCPM[rownames(logCPM) %in% rownames(qlf_DE_Genes),]

#Transpose and convert to data frame
logCPM <- data.frame(t(logCPM),check.names = F)

#Add a group column to logCPM
logCPM$Group <- "RCM"
logCPM$Group[rownames(logCPM) %in% c("P0197","P0256","P0144","P0142")] <- "Control"

#Split off RCM and control groups
control <- logCPM[logCPM$Group == "Control",]
RCM <- logCPM[logCPM$Group == "RCM",]

#Remove group column and transpose data frames
control$Group <- NULL
control <- data.frame(t(control),check.names = F)
RCM$Group <- NULL
RCM <- data.frame(t(RCM),check.names = F)

#Calculate mean log2 logCPM
control$Mean <- rowMeans(control)
RCM$Mean <- rowMeans(RCM)

#Get mean columns and remove them from data frames
control_Mean <- control$Mean
RCM_Mean <- RCM$Mean

control$Mean <- NULL
RCM$Mean <- NULL

#Convert data frames to matrices
control <- as.matrix(control)
RCM <- as.matrix(RCM)

#Calculate log2 FC for each RCM sample over the mean of control samples, and vice versa
control_FC <- sweep(control, 1, RCM_Mean)
RCM_FC <- sweep(RCM, 1, control_Mean)

#Convert matrices back to data frames, add a group column, and combine them
control_FC <- t(control_FC)
control_FC <- data.frame(Group="Control",control_FC,check.names = F)

RCM_FC <- t(RCM_FC)
RCM_FC <- data.frame(Group="RCM",RCM_FC,check.names = F)

All_FC <- rbind(RCM_FC,control_FC)

col_annotation <- data.frame(Group=as.factor(All_FC$Group))
rownames(col_annotation) <- rownames(All_FC)

All_FC$Group <- NULL
All_FC <- data.frame(t(All_FC),check.names = F)

write.table(data.frame(GeneSymbol=rownames(All_FC),All_FC),"RCM_and_Control_FC_over_Mean_FDR_0.05.txt",sep="\t",row.names=F,col.names = T,quote=F)

my_colors = list(Group=c(RCM="#0000FF",Control="#006100"))

heatmap <- pheatmap(All_FC,cluster_cols = F,color=colorRampPalette(c("navy", "white", "red"))(50),annotation_col = col_annotation,annotation_colors = my_colors)

dev.off()

order <- heatmap$tree_row[["order"]]
All_FC <- All_FC[order,]

write.table(data.frame(GeneSymbol=rownames(All_FC),All_FC,check.names=F),"Heatmap_1931_Genes_FDR_0.05.txt",sep="\t",row.names=T,col.names=T,quote=F)

pdf("Heatmap_log2FC_1931_DE_Genes_FDR_0.05.pdf",width=10,height=10)
pheatmap(All_FC,cluster_cols = F,color=colorRampPalette(c("navy", "white", "red"))(50),annotation_col = col_annotation,annotation_colors = my_colors,show_rownames=FALSE,annotation_names_col = F)
dev.off()

#Make a heatmap with a subset of DE genes (|log2 FC| >= 1.5 and FDR adjusted p-val < 0.05)####

#Keep genes with |log2 FC| >= 1.5
qlf_DE_Genes_Sub <- qlf_DE_Genes[abs(qlf_DE_Genes$logFC) >= 1.5,]

Sub_FC <- All_FC[rownames(All_FC) %in% rownames(qlf_DE_Genes_Sub),]

heatmap <- pheatmap(Sub_FC,cluster_cols = F,color=colorRampPalette(c("navy", "white", "red"))(50),annotation_col = col_annotation,annotation_colors = my_colors)

dev.off()

order <- heatmap$tree_row[["order"]]
Sub_FC <- Sub_FC[order,]

write.table(data.frame(GeneSymbol=rownames(Sub_FC),Sub_FC,check.names = F),"Heatmap_413_Genes_logFC_1.5.txt",sep="\t",row.names=F,col.names=T,quote=F)

pdf("Heatmap_log2FC_413_DE_Genes_logFC_1.5.pdf",width=10,height=10)
pheatmap(Sub_FC,cluster_cols = F,color=colorRampPalette(c("navy", "white", "red"))(50),annotation_col = col_annotation,annotation_colors = my_colors,show_rownames=FALSE,annotation_names_col = F)
dev.off()

#Make a heatmap with the top 60 genes, based on p-value####

#Read in 60 genes with most divergent expression pattern between RCM and control
Sub_FC <- read.delim("Heatmap_top_60_Genes_pval.txt")
rownames(Sub_FC) <- Sub_FC$GeneSymbol
Sub_FC$GeneSymbol <- NULL

pdf("Heatmap_log2FC_Top_60_Genes_pval.pdf",width=10,height=10)
pheatmap(Sub_FC,cluster_cols = F,color=colorRampPalette(c("navy", "white", "red"))(50),annotation_col = col_annotation,annotation_colors = my_colors,show_rownames=TRUE,annotation_names_col = F)
dev.off()

#Make a heatmap with the top 60 genes, based on fold-change####

#Read in 60 genes with most divergent expression pattern between RCM and control
Sub_FC <- read.delim("Heatmap_top_60_Genes_FC.txt")
rownames(Sub_FC) <- Sub_FC$GeneSymbol
Sub_FC$GeneSymbol <- NULL

pdf("Heatmap_log2FC_Top_60_Genes_FC.pdf",width=10,height=10)
pheatmap(Sub_FC,cluster_cols = F,color=colorRampPalette(c("navy", "white", "red"))(50),annotation_col = col_annotation,annotation_colors = my_colors,show_rownames=TRUE,annotation_names_col = F)
dev.off()

#Make a heatmap with the top 55 lncRNAs, based on p-value####

#Read in fold-change values for 55 lncRNAs
Sub_FC <- read.delim("lncRNA_FC_Heatmap.txt")
rownames(Sub_FC) <- Sub_FC$GeneSymbol

Sub_FC <- Sub_FC[,c(15:32)]

pdf("Heatmap_log2FC_Top_55_lncRNA_pval.pdf",width=10,height=10)
pheatmap(Sub_FC,cluster_cols = F,color=colorRampPalette(c("navy", "white", "red"))(50),annotation_col = col_annotation,annotation_colors = my_colors,show_rownames=TRUE,annotation_names_col = F)
dev.off()

