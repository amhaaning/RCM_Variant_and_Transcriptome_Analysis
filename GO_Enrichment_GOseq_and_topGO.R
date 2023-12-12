library(topGO)
library(dplyr)
library(org.Hs.eg.db)
library(AnnotationDbi)
library(reshape2)
library(pheatmap)

setwd("")

#Read in expression data
data <- read.delim("EdgeR_results_qlf_test_all_genes.txt")

#Keep significant DE genes (FDR < 0.05)
DE_Genes <- data[data$FDR < 0.05,]

#Make a vector to indicate whether gene is differentially expressed
genes <- rep(0,length(data$GeneSymbol))
names(genes) <- data$GeneSymbol
genes[names(genes) %in% DE_Genes$GeneSymbol] <- 1

#topGO analysis

#Make character vector with all GO categories
GO_Cats <- c("MF","BP","CC")

#Down-regulated genes
down_genes <- data[data$FDR < 0.05 & data$logFC < 0,]

#Up-regulated genes
up_genes <- data[data$FDR < 0.05 & data$logFC > 0,]

#GO Term enrichment for all down-regulated genes####

for (i in 1:length(GO_Cats)){
    
  #Set a GO category (BP, CC, or MF)
  GO_Cat <- GO_Cats[i]
    
  #Get GO terms for specific category
  geneList <- annFUN.org(GO_Cat,mapping="org.Hs.eg.db",ID="SYMBOL")
  allGenes <- unique(unlist(geneList))
  geneList <- factor(as.integer(allGenes %in% down_genes$GeneSymbol))
  names(geneList) <- allGenes
    
  #Create a new topGO data container for enrichment analysis
  topGO <- new("topGOdata", ontology=GO_Cat, allGenes=geneList, annot=annFUN.org, mapping="org.Hs.eg.db", ID="SYMBOL",nodeSize=5)
    
  #Perform statistical tests for enrichment
  classic_fisher_result=runTest(topGO, algorithm='classic', statistic='fisher')
  weight_fisher_result=runTest(topGO, algorithm='weight01', statistic='fisher')
    
  # generate a table of results: we can use the GenTable function to generate a summary table with the results from tests applied to the topGOdata object.
  allGO=usedGO(topGO)
  all_res=GenTable(topGO, weightFisher=weight_fisher_result, weightClassicFisher=classic_fisher_result, orderBy='weightFisher', topNodes=length(allGO),numChar=1000)
    all_res$weightFisher_adj <- p.adjust(all_res$weightFisher,method="fdr")
    
  #Get significant results and write to file
  results.table.p= all_res[which(all_res$weightFisher_adj<=0.05),]
  
  if (nrow(results.table.p) > 0){
    #Get genes associated with each GO term
    genes_in_go <- genesInTerm(topGO)
    #Keep significant GO terms
    genes_in_go <- genes_in_go[names(genes_in_go) %in% results.table.p$GO.ID]
    
    #convert the list of genes associated with GO terms to a data frame
    n.obs <- sapply(genes_in_go, length)
    seq.max <- seq_len(max(n.obs))
    mat <- t(sapply(genes_in_go, "[", i = seq.max))
    genes_in_go <- data.frame(GO.ID=rownames(mat),mat,check.names = F)
    
    #Convert to long format
    genes_in_go <- melt(genes_in_go,id="GO.ID")
    genes_in_go$variable <- NULL
    colnames(genes_in_go)[2] <- "Associated_DE_Gene"
    
    #Keep genes associated with GO terms if they are in the list of DE genes
    genes_in_go <- genes_in_go[genes_in_go$Associated_DE_Gene %in% names(geneList)[geneList==1],]
    
    #Add descriptions of GO terms
    genes_in_go <- merge(results.table.p[,c("GO.ID","Term")],genes_in_go,by="GO.ID")
    
    write.table(results.table.p,file=paste0("./topGO_Enrichment_Results/Downregulated_Genes/","RCM_Downregulated_Genes_",GO_Cat,"_topGO_enrichment.txt"),sep="\t",row.names=F,col.names=T,quote=F)
    write.table(genes_in_go,file=paste0("./topGO_Enrichment_Results/Downregulated_Genes/","RCM_Downregulated_Genes_",GO_Cat,"_associated_DE_Genes.txt"),sep="\t",row.names=F,col.names=T,quote=F)
  }
}

#GO Term enrichment for all up-regulated genes####

for (i in 1:length(GO_Cats)){
  
  #Set a GO category (BP, CC, or MF)
  GO_Cat <- GO_Cats[i]
  
  #Get GO terms for specific category
  geneList <- annFUN.org(GO_Cat,mapping="org.Hs.eg.db",ID="SYMBOL")
  allGenes <- unique(unlist(geneList))
  geneList <- factor(as.integer(allGenes %in% up_genes$GeneSymbol))
  names(geneList) <- allGenes
  
  #Create a new topGO data container for enrichment analysis
  topGO <- new("topGOdata", ontology=GO_Cat, allGenes=geneList, annot=annFUN.org, mapping="org.Hs.eg.db", ID="SYMBOL",nodeSize=5)
  
  #Perform statistical tests for enrichment
  classic_fisher_result=runTest(topGO, algorithm='classic', statistic='fisher')
  weight_fisher_result=runTest(topGO, algorithm='weight01', statistic='fisher')
  
  # generate a table of results: we can use the GenTable function to generate a summary table with the results from tests applied to the topGOdata object.
  allGO=usedGO(topGO)
  all_res=GenTable(topGO, weightFisher=weight_fisher_result, weightClassicFisher=classic_fisher_result, orderBy='weightFisher', topNodes=length(allGO),numChar=1000)
  all_res$weightFisher_adj <- p.adjust(all_res$weightFisher,method="fdr")
  
  #Get significant results and write to file
  results.table.p= all_res[which(all_res$weightFisher_adj<=0.05),]
  
  if (nrow(results.table.p) > 0){
    #Get genes associated with each GO term
    genes_in_go <- genesInTerm(topGO)
    #Keep significant GO terms
    genes_in_go <- genes_in_go[names(genes_in_go) %in% results.table.p$GO.ID]
    
    #convert the list of genes associated with GO terms to a data frame
    n.obs <- sapply(genes_in_go, length)
    seq.max <- seq_len(max(n.obs))
    mat <- t(sapply(genes_in_go, "[", i = seq.max))
    genes_in_go <- data.frame(GO.ID=rownames(mat),mat,check.names = F)
    
    #Convert to long format
    genes_in_go <- melt(genes_in_go,id="GO.ID")
    genes_in_go$variable <- NULL
    colnames(genes_in_go)[2] <- "Associated_DE_Gene"
    
    #Keep genes associated with GO terms if they are in the list of DE genes
    genes_in_go <- genes_in_go[genes_in_go$Associated_DE_Gene %in% names(geneList)[geneList==1],]
    
    #Add descriptions of GO terms
    genes_in_go <- merge(results.table.p[,c("GO.ID","Term")],genes_in_go,by="GO.ID")
    
    write.table(results.table.p,file=paste0("./topGO_Enrichment_Results/Upregulated_Genes/","RCM_Upregulated_Genes_",GO_Cat,"_topGO_enrichment.txt"),sep="\t",row.names=F,col.names=T,quote=F)
    write.table(genes_in_go,file=paste0("./topGO_Enrichment_Results/Upregulated_Genes/","RCM_Upregulated_Genes_",GO_Cat,"_associated_DE_Genes.txt"),sep="\t",row.names=F,col.names=T,quote=F)
  }
}

#GO Term enrichment for all DE genes####

for (i in 1:length(GO_Cats)){
  
  #Set a GO category (BP, CC, or MF)
  GO_Cat <- GO_Cats[i]
  
  #Get GO terms for specific category
  geneList <- annFUN.org(GO_Cat,mapping="org.Hs.eg.db",ID="SYMBOL")
  allGenes <- unique(unlist(geneList))
  geneList <- factor(as.integer(allGenes %in% DE_Genes$GeneSymbol))
  names(geneList) <- allGenes
  
  #Create a new topGO data container for enrichment analysis
  topGO <- new("topGOdata", ontology=GO_Cat, allGenes=geneList, annot=annFUN.org, mapping="org.Hs.eg.db", ID="SYMBOL",nodeSize=5)
  
  #Perform statistical tests for enrichment
  classic_fisher_result=runTest(topGO, algorithm='classic', statistic='fisher')
  weight_fisher_result=runTest(topGO, algorithm='weight01', statistic='fisher')
  
  # generate a table of results: we can use the GenTable function to generate a summary table with the results from tests applied to the topGOdata object.
  allGO=usedGO(topGO)
  all_res=GenTable(topGO, weightFisher=weight_fisher_result, weightClassicFisher=classic_fisher_result, orderBy='weightFisher', topNodes=length(allGO),numChar=1000)
  all_res$weightFisher_adj <- p.adjust(all_res$weightFisher,method="fdr")
  
  #Get significant results and write to file
  results.table.p= all_res[which(all_res$weightFisher_adj<=0.05),]
  
  if (nrow(results.table.p) > 0){
    #Get genes associated with each GO term
    genes_in_go <- genesInTerm(topGO)
    #Keep significant GO terms
    genes_in_go <- genes_in_go[names(genes_in_go) %in% results.table.p$GO.ID]
    
    #convert the list of genes associated with GO terms to a data frame
    n.obs <- sapply(genes_in_go, length)
    seq.max <- seq_len(max(n.obs))
    mat <- t(sapply(genes_in_go, "[", i = seq.max))
    genes_in_go <- data.frame(GO.ID=rownames(mat),mat,check.names = F)
    
    #Convert to long format
    genes_in_go <- melt(genes_in_go,id="GO.ID")
    genes_in_go$variable <- NULL
    colnames(genes_in_go)[2] <- "Associated_DE_Gene"
    
    #Keep genes associated with GO terms if they are in the list of DE genes
    genes_in_go <- genes_in_go[genes_in_go$Associated_DE_Gene %in% names(geneList)[geneList==1],]
    
    #Add descriptions of GO terms
    genes_in_go <- merge(results.table.p[,c("GO.ID","Term")],genes_in_go,by="GO.ID")
    
    write.table(results.table.p,file=paste0("./topGO_Enrichment_Results/All_DE_Genes/","RCM_All_DE_Genes_",GO_Cat,"_topGO_enrichment.txt"),sep="\t",row.names=F,col.names=T,quote=F)
    write.table(genes_in_go,file=paste0("./topGO_Enrichment_Results/All_DE_Genes/","RCM_All_DE_Genes_",GO_Cat,"_associated_DE_Genes.txt"),sep="\t",row.names=F,col.names=T,quote=F)
  }
}

#Combine results from topGO and add DE genes that are associated with significant GO terms####

#Read in results from EdgeR
edgeR_results <- read.delim("EdgeR_results_qlf_test_all_genes.txt")

#Keep necessary columns
colnames(edgeR_results)[2] <- "ENTREZID"
edgeR_results$Chr <- NULL
edgeR_results$Start <- NULL
edgeR_results$End <- NULL
edgeR_results$Strand <- NULL
edgeR_results$Length <- NULL

#Read in results from topGO analysis
topGO_MF <- read.delim("./topGO_Enrichment_Results/All_DE_Genes/RCM_All_DE_Genes_MF_topGO_enrichment.txt")

topGO_CC <- read.delim("./topGO_Enrichment_Results/All_DE_Genes/RCM_All_DE_Genes_CC_topGO_enrichment.txt")

#Read in genes associated with significantly overrepresented GO terms
topGO_MF_Genes <- read.delim("./topGO_Enrichment_Results/All_DE_Genes/RCM_All_DE_Genes_MF_associated_DE_Genes.txt")

topGO_MF_Genes$Term <- NULL

topGO_CC_Genes <- read.delim("./topGO_Enrichment_Results/All_DE_Genes/RCM_All_DE_Genes_CC_associated_DE_Genes.txt")

topGO_CC_Genes$Term <- NULL

#Merge all results together
topGO_MF <- merge(topGO_MF_Genes,topGO_MF,by="GO.ID")
colnames(topGO_MF)[2] <- "GeneSymbol"
topGO_MF <- merge(topGO_MF,edgeR_results,by="GeneSymbol")

topGO_CC <- merge(topGO_CC_Genes,topGO_CC,by="GO.ID")
colnames(topGO_CC)[2] <- "GeneSymbol"
topGO_CC <- merge(topGO_CC,edgeR_results,by="GeneSymbol")

#Order results by edgeR p-value
topGO_MF <- topGO_MF[order(topGO_MF$FDR),]
topGO_CC <- topGO_CC[order(topGO_CC$FDR),]

#Add a column stating whether gene is up- or down-regulated
topGO_MF$Direction <- NA
topGO_MF$Direction[topGO_MF$logFC < 0] <- "Downregulated"
topGO_MF$Direction[topGO_MF$logFC > 0] <- "Upregulated"

topGO_CC$Direction <- NA
topGO_CC$Direction[topGO_CC$logFC < 0] <- "Downregulated"
topGO_CC$Direction[topGO_CC$logFC > 0] <- "Upregulated"

#Add a column for the type of GO term
topGO_MF$GO.TYPE <- "Molecular Function"
topGO_CC$GO.TYPE <- "Cellular Compartment"

#Combine all results
topGO <- rbind(topGO_MF,topGO_CC)

#Add a column indicating the number of times a gene is associated with a GO term
Count <- data.frame(table(topGO$GeneSymbol))
colnames(Count) <- c("GeneSymbol","Count")
topGO <- merge(topGO,Count,by="GeneSymbol")

#Write results to table
write.table(topGO,"./topGO_Enrichment_Results/topGO_edgeR_results_merged.txt",sep="\t",quote=F,row.names=F,col.names=T)

