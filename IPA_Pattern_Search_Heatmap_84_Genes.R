library(pheatmap)

setwd("C:/Users/ahaaning/OneDrive - Indiana University/Manuscripts/RCM_Variants_and_RNAseq/IPA_Analysis/")

#Read in data from IPA from heatmap
Heatmap_data <- read.delim("Heatmap_Data.txt",check.names = F)

#Format column names for visualization
colnames(Heatmap_data) <- gsub(" vs.*\\("," \\(",colnames(Heatmap_data))
colnames(Heatmap_data) <- gsub(", idiopathic","",colnames(Heatmap_data))

#Read in log2 FC values for all genes and keep the genes from the IPA data
RCM_DE_data <- read.delim("../RCM_RNAseq/RCM_and_Control_FC_over_Mean_FDR_0.05.txt")
RCM_DE_data <- RCM_DE_data[RCM_DE_data$GeneSymbol %in% Heatmap_data$GeneSymbol,]

#Keep RCM samples only
RCM_DE_data <- RCM_DE_data[,c(1:15)]

#Combine RCM data with IPA data
All_data <- merge(RCM_DE_data,Heatmap_data,by="GeneSymbol")

#Set gene symbols as rownames
rownames(All_data) <- All_data$GeneSymbol

paletteLength <- 50
myColor <- colorRampPalette(c("navy", "white", "darkred"))(paletteLength)
#length(breaks) == length(paletteLength) + 1
#use floor and ceiling to deal with even/odd length pallettelengths

myBreaks <- c(seq(min(All_data[,-1],na.rm=T), 0, length.out=ceiling(paletteLength/2+ 1)),seq(max(All_data[,-1],na.rm=T)/paletteLength, max(All_data[,-1],na.rm=T), length.out=floor(paletteLength/2)))


#Make heatmap showing log2 FC values for RCM samples and studies identified through IPA pattern search
pdf("Heatmap_84_Genes_IPA_Pattern_Search.pdf",height=11,width=9)
heatmap <- pheatmap(All_data[,-1],color=myColor,breaks=myBreaks,na_col = "gray40",treeheight_col = 0,cluster_cols = F)
dev.off()

#Read in log2 FC from EdgeR
logFC_edgeR <- read.delim("../RCM_RNAseq/EdgeR_results_qlf_test_all_genes.txt")

#Keep 84 genes
logFC_edgeR <- logFC_edgeR[logFC_edgeR$GeneSymbol %in% Heatmap_data$GeneSymbol,]

#Keep necessary columns
logFC_edgeR <- logFC_edgeR[,c("GeneSymbol","logFC")]

#Combine RCM data with IPA data
All_data <- merge(logFC_edgeR,Heatmap_data,by="GeneSymbol")

#Set gene symbols as rownames
rownames(All_data) <- All_data$GeneSymbol

#Change column name for RCM log2 FC
colnames(All_data)[2] <- "RCM (current study)"

myBreaks <- c(seq(min(All_data[,-1],na.rm=T), 0, length.out=ceiling(paletteLength/2+ 1)),seq(max(All_data[,-1],na.rm=T)/paletteLength, max(All_data[,-1],na.rm=T), length.out=floor(paletteLength/2)))

#Make heatmap showing log2 FC values for RCM compared to controls and studies identified through IPA pattern search
pdf("Heatmap_84_Genes_IPA_Pattern_Search_EdgeR_log2FC.pdf",height=11,width=7)
heatmap <- pheatmap(All_data[,-1],color=myColor,breaks=myBreaks,na_col = "gray40",treeheight_col = 0,cluster_cols = F)
dev.off()

#Write all results for final heatmap to a file
write.table(All_data,"Final_Heatmap_Data_Fig_S2.txt",row.names = F,col.names = T,sep="\t",quote=F)
