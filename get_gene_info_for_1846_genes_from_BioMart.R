library(biomaRt)

ensembl <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl")
attributes <- listAttributes(ensembl)

attributes[grepl("hgnc",attributes$name,ignore.case = T),]

searchDatasets(mart = ensembl, pattern = "hsapiens")

setwd("/N/project/WareLab_ARP_NGS/RCM_VCFs/")

genes <- read.delim("1846_genes.txt",header=F)
genes <- genes[!duplicated(genes$V1),]

human_genes <- getBM(attributes=c('ensembl_gene_id','external_gene_name','external_synonym','chromosome_name','start_position','end_position','strand','gene_biotype','description','hgnc_id'),mart = ensembl)

colnames(human_genes) <- c('ensembl_gene_id','hgnc_gene_symbol','gene_symbol_synonym','chromosome','start','end','strand','gene_biotype','description','hgnc_id')

write.table(human_genes,"Human_Genes_BioMart.txt",sep="\t",row.names = F,col.names=T,quote = F)

ensembl_IDs_1846_genes <- human_genes[human_genes$hgnc_gene_symbol %in% genes | human_genes$ensembl_gene_id %in% genes,]

ensembl_IDs_1846_genes <- ensembl_IDs_1846_genes[ensembl_IDs_1846_genes$chromosome %in% c("1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","X","Y","MT"),]

ensembl_IDs_1846_genes <- ensembl_IDs_1846_genes[!duplicated(ensembl_IDs_1846_genes$ensembl_gene_id),]

#More than 1846 genes - check to see if hgnc gene symbols are duplicated
duplicated_symbols <- ensembl_IDs_1846_genes$hgnc_gene_symbol[duplicated(ensembl_IDs_1846_genes$hgnc_gene_symbol)]

#Two genes with duplicated symbols - MKKS and MATR3. I'm going to keep all Ensemble gene IDs for these symbols.

#Write all info for 1846 genes to a file
write.table(ensembl_IDs_1846_genes,file="All_Info_1846_Genes.txt",sep="\t",quote=F,row.names=F,col.names=T)

#Get ensembl IDs only and write them to a file
ensembl_IDs_1846_genes <- ensembl_IDs_1846_genes$ensembl_gene_id

write.table(ensembl_IDs_1846_genes,"ensembl_IDs_1846_genes.txt",sep="\t",quote=F,row.names=F,col.names=F)
