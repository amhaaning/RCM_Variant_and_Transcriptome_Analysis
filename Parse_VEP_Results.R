#Set working directory
setwd("/N/project/WareLab_ARP_NGS/RCM_VCFs/")

#Read in Ensembl gene IDs
genes <- read.delim("ensembl_IDs_1846_genes.txt",header=F)

#INDELs####

#Read in file with variant effects from VEP
variant_effects <- read.delim("RCM_WES_Apr18_Joint_INDELs_FINAL_1846_Genes.txt",skip = 107,header=T)
variants <- read.delim("RCM_WES_Apr18_Joint_INDELs_FINAL_1846_Genes_Variants.txt")

variant_effects[variant_effects=="-"] <- NA
variant_effects$Allele[is.na(variant_effects$Allele)] <- "-"

#Create columns for merging and merge
chrom <- variant_effects$X.Uploaded_variation
chrom <- gsub("_.*","",chrom)
pos <- variant_effects$X.Uploaded_variation
pos <- gsub("chr.._","",pos)
pos <- gsub("chr._","",pos)
pos <- gsub("_.*","",pos)
pos <- (as.numeric(pos))-1
variant_effects$merge_col <- paste(chrom,pos,sep="_")
variants$merge_col <- paste(variants$CHROM,variants$POS,sep="_")

variant_effects <- merge(variant_effects,variants,by="merge_col")

#Make sure gene is in the cadidate gene list
variant_effects <- variant_effects[variant_effects$Gene %in% genes$V1,]

#Change first column name
colnames(variant_effects)[2] <- "Variant"

#Parse SpliceAI Delta score so that it can be used for sorting
variant_effects$SpliceAI_Delta <- sub(".*?\\|","",variant_effects$SpliceAI_pred)
variant_effects$SpliceAI_Delta <- gsub("\\|.*","",variant_effects$SpliceAI_Delta)
variant_effects$SpliceAI_Delta <- as.numeric(variant_effects$SpliceAI_Delta)

#Add single apostraphe to prevent excel from converting values to dates
variant_effects$EXON <- paste0("'",variant_effects$EXON)
variant_effects$INTRON <- paste0("'",variant_effects$INTRON)

#Write to a file
write.table(variant_effects,"RCM_WES_Apr18_Joint_INDELs_1846_CM_Genes_Final.txt",sep = "\t",quote=F,row.names = F,col.names = T)

#Number of unique variants 
length(unique(variant_effects$Variant))

#Get moderate impact variants
moderate_variants <- variant_effects[variant_effects$IMPACT=="MODERATE",]

#Number of moderate impact variants
length(unique(moderate_variants$Variant))

#Get high impact variants
high_variants <- variant_effects[variant_effects$IMPACT=="HIGH",]

#Number of high impact variants
length(unique(high_variants$Variant))

#write to files
write.table(moderate_variants,"RCM_WES_Apr18_Joint_INDELs_1846_CM_Genes_Moderate_Impact.txt",sep="\t",quote=F,row.names=F,col.names=T)

write.table(high_variants,"RCM_WES_Apr18_Joint_INDELs_1846_CM_Genes_High_Impact.txt",sep="\t",quote=F,row.names=F,col.names=T)

#Create a short list of unique variants for quick reference
moderate_variants <- moderate_variants[,c("Variant","Location","Allele","Gene","Consequence","Existing_variation","IMPACT","STRAND","VARIANT_CLASS","SYMBOL","HGNC_ID","MANE_SELECT","TSL","SWISSPROT","AF","gnomADe_AF","gnomADg_AF","MAX_AF","CLIN_SIG","CADD_PHRED","CADD_RAW","SpliceAI_pred","SpliceAI_Delta","IU.0105","IU.0106","P0142","P0144","P0177","P0178","P0197","P0199","P0203","P0215","P0220","P0223","P0233","P0247","P0256","P0288","P0328","P0339")]

high_variants <- high_variants[,c("Variant","Location","Allele","Gene","Consequence","Existing_variation","IMPACT","STRAND","VARIANT_CLASS","SYMBOL","HGNC_ID","MANE_SELECT","TSL","SWISSPROT","AF","gnomADe_AF","gnomADg_AF","MAX_AF","CLIN_SIG","CADD_PHRED","CADD_RAW","SpliceAI_pred","SpliceAI_Delta","IU.0105","IU.0106","P0142","P0144","P0177","P0178","P0197","P0199","P0203","P0215","P0220","P0223","P0233","P0247","P0256","P0288","P0328","P0339")]

#Remove redundancies
moderate_variants <- moderate_variants[!duplicated(moderate_variants$Variant),]
high_variants <- high_variants[!duplicated(high_variants$Variant),]

#Write to files
write.table(moderate_variants,"RCM_WES_Apr18_Joint_INDELs_1846_CM_Genes_Moderate_Impact_Short.txt",sep="\t",quote=F,row.names=F,col.names=T)

write.table(high_variants,"RCM_WES_Apr18_Joint_INDELs_1846_CM_Genes_High_Impact_Short.txt",sep="\t",quote=F,row.names=F,col.names=T)

#SNVs####

#Read in file with variant effects from VEP
variant_effects <- read.delim("RCM_WES_Apr18_Joint_SNVs_FINAL_1846_Genes.txt",skip = 110,header=T)
variants <- read.delim("RCM_WES_Apr18_Joint_SNVs_FINAL_1846_Genes_Variants.txt")

variant_effects[variant_effects=="-"] <- NA

#Create columns for merging and merge
chrom <- variant_effects$X.Uploaded_variation
chrom <- gsub("_.*","",chrom)
pos <- variant_effects$X.Uploaded_variation
pos <- gsub("chr.._","",pos)
pos <- gsub("chr._","",pos)
pos <- gsub("_.*","",pos)
pos <- (as.numeric(pos))
variant_effects$merge_col <- paste(chrom,pos,sep="_")
variants$merge_col <- paste(variants$CHROM,variants$POS,sep="_")

variant_effects <- merge(variant_effects,variants,by="merge_col")

#Make sure gene is in the cadidate gene list
variant_effects <- variant_effects[variant_effects$Gene %in% genes$V1,]

#Change first column name
colnames(variant_effects)[2] <- "Variant"

#Parse SpliceAI Delta score so that it can be used for sorting
variant_effects$SpliceAI_Delta <- sub(".*?\\|","",variant_effects$SpliceAI_pred)
variant_effects$SpliceAI_Delta <- gsub("\\|.*","",variant_effects$SpliceAI_Delta)
variant_effects$SpliceAI_Delta <- as.numeric(variant_effects$SpliceAI_Delta)

#Add single apostraphe to prevent excel from converting values to dates
variant_effects$EXON <- paste0("'",variant_effects$EXON)
variant_effects$INTRON <- paste0("'",variant_effects$INTRON)

#Write to a file
write.table(variant_effects,"RCM_WES_Apr18_Joint_SNVs_1846_CM_Genes_Final.txt",sep = "\t",quote=F,row.names = F,col.names = T)

#Get moderate impact variants
moderate_variants <- variant_effects[variant_effects$IMPACT=="MODERATE",]

#Number of moderate impact variants
length(unique(moderate_variants$Variant))

#Get high impact variants
high_variants <- variant_effects[variant_effects$IMPACT=="HIGH",]

#Number of high impact variants
length(unique(high_variants$Variant))

#write to files
write.table(moderate_variants,"RCM_WES_Apr18_Joint_SNVs_1846_CM_Genes_Moderate_Impact.txt",sep="\t",quote=F,row.names=F,col.names=T)

write.table(high_variants,"RCM_WES_Apr18_Joint_SNVs_1846_CM_Genes_High_Impact.txt",sep="\t",quote=F,row.names=F,col.names=T)

#Create a short list of unique variants for quick reference
moderate_variants <- moderate_variants[,c("Variant","Location","Allele","Gene","Consequence","Existing_variation","IMPACT","STRAND","VARIANT_CLASS","SYMBOL","HGNC_ID","MANE_SELECT","TSL","SWISSPROT","AF","gnomADe_AF","gnomADg_AF","MAX_AF","CLIN_SIG","CADD_PHRED","CADD_RAW","SpliceAI_pred","SpliceAI_Delta","PrimateAI","am_class","am_pathogenicity","SIFT","PolyPhen","IU.0105","IU.0106","P0142","P0144","P0177","P0178","P0197","P0199","P0203","P0215","P0220","P0223","P0233","P0247","P0256","P0288","P0328","P0339")]

high_variants <- high_variants[,c("Variant","Location","Allele","Gene","Consequence","Existing_variation","IMPACT","STRAND","VARIANT_CLASS","SYMBOL","HGNC_ID","MANE_SELECT","TSL","SWISSPROT","AF","gnomADe_AF","gnomADg_AF","MAX_AF","CLIN_SIG","CADD_PHRED","CADD_RAW","SpliceAI_pred","SpliceAI_Delta","PrimateAI","am_class","am_pathogenicity","SIFT","PolyPhen","IU.0105","IU.0106","P0142","P0144","P0177","P0178","P0197","P0199","P0203","P0215","P0220","P0223","P0233","P0247","P0256","P0288","P0328","P0339")]

#Remove redundancies
moderate_variants <- moderate_variants[!duplicated(moderate_variants$Variant),]
high_variants <- high_variants[!duplicated(high_variants$Variant),]

#Write to files
write.table(moderate_variants,"RCM_WES_Apr18_Joint_SNVs_1846_CM_Genes_Moderate_Impact_Short.txt",sep="\t",quote=F,row.names=F,col.names=T)

write.table(high_variants,"RCM_WES_Apr18_Joint_SNVs_1846_CM_Genes_High_Impact_Short.txt",sep="\t",quote=F,row.names=F,col.names=T)

#INDELs - All####

#Read in file with variant effects from VEP
variant_effects <- read.delim("RCM_WES_Apr18_Joint_INDELs_FINAL_Anno.txt",skip = 107,header=T)
variants <- read.delim("RCM_WES_Apr18_Joint_INDELs_FINAL_Anno_Variants.txt")

variant_effects[variant_effects=="-"] <- NA
variant_effects$Allele[is.na(variant_effects$Allele)] <- "-"

#Create columns for merging and merge
chrom <- variant_effects$X.Uploaded_variation
chrom <- gsub("_.*","",chrom)
pos <- variant_effects$X.Uploaded_variation
pos <- gsub("chr.._","",pos)
pos <- gsub("chr._","",pos)
pos <- gsub("_.*","",pos)
pos <- (as.numeric(pos))-1
variant_effects$merge_col <- paste(chrom,pos,sep="_")
variants$merge_col <- paste(variants$CHROM,variants$POS,sep="_")

variant_effects <- merge(variant_effects,variants,by="merge_col")

#Change first column name
colnames(variant_effects)[2] <- "Variant"

#Parse SpliceAI Delta score so that it can be used for sorting
variant_effects$SpliceAI_Delta <- sub(".*?\\|","",variant_effects$SpliceAI_pred)
variant_effects$SpliceAI_Delta <- gsub("\\|.*","",variant_effects$SpliceAI_Delta)
variant_effects$SpliceAI_Delta <- as.numeric(variant_effects$SpliceAI_Delta)

#Add single apostraphe to prevent excel from converting values to dates
variant_effects$EXON <- paste0("'",variant_effects$EXON)
variant_effects$INTRON <- paste0("'",variant_effects$INTRON)

#Write to a file
write.table(variant_effects,"RCM_WES_Apr18_Joint_INDELs_All_Genes_Final.txt",sep = "\t",quote=F,row.names = F,col.names = T)

#Number of unique variants 
length(unique(variant_effects$Variant))

#Get moderate impact variants
moderate_variants <- variant_effects[variant_effects$IMPACT=="MODERATE",]

#Number of moderate impact variants
length(unique(moderate_variants$Variant))

#Get high impact variants
high_variants <- variant_effects[variant_effects$IMPACT=="HIGH",]

#Number of high impact variants
length(unique(high_variants$Variant))

#write to files
write.table(moderate_variants,"RCM_WES_Apr18_Joint_INDELs_All_Genes_Moderate_Impact.txt",sep="\t",quote=F,row.names=F,col.names=T)

write.table(high_variants,"RCM_WES_Apr18_Joint_INDELs_All_High_Impact.txt",sep="\t",quote=F,row.names=F,col.names=T)

#Create a short list of unique variants for quick reference
moderate_variants <- moderate_variants[,c("Variant","Location","Allele","Gene","Consequence","Existing_variation","IMPACT","STRAND","VARIANT_CLASS","SYMBOL","HGNC_ID","MANE_SELECT","TSL","SWISSPROT","AF","gnomADe_AF","gnomADg_AF","MAX_AF","CLIN_SIG","CADD_PHRED","CADD_RAW","SpliceAI_pred","SpliceAI_Delta","IU.0105","IU.0106","P0142","P0144","P0177","P0178","P0197","P0199","P0203","P0215","P0220","P0223","P0233","P0247","P0256","P0288","P0328","P0339")]

high_variants <- high_variants[,c("Variant","Location","Allele","Gene","Consequence","Existing_variation","IMPACT","STRAND","VARIANT_CLASS","SYMBOL","HGNC_ID","MANE_SELECT","TSL","SWISSPROT","AF","gnomADe_AF","gnomADg_AF","MAX_AF","CLIN_SIG","CADD_PHRED","CADD_RAW","SpliceAI_pred","SpliceAI_Delta","IU.0105","IU.0106","P0142","P0144","P0177","P0178","P0197","P0199","P0203","P0215","P0220","P0223","P0233","P0247","P0256","P0288","P0328","P0339")]

#Remove redundancies
moderate_variants <- moderate_variants[!duplicated(moderate_variants$Variant),]
high_variants <- high_variants[!duplicated(high_variants$Variant),]

#Write to files
write.table(moderate_variants,"RCM_WES_Apr18_Joint_INDELs_All_Genes_Moderate_Impact_Short.txt",sep="\t",quote=F,row.names=F,col.names=T)

write.table(high_variants,"RCM_WES_Apr18_Joint_INDELs_All_Genes_High_Impact_Short.txt",sep="\t",quote=F,row.names=F,col.names=T)

#SNVs - All####

#Read in file with variant effects from VEP
variant_effects <- read.delim("RCM_WES_Apr18_Joint_SNVs_FINAL_Anno.txt",skip = 110,header=T)
variants <- read.delim("RCM_WES_Apr18_Joint_SNVs_FINAL_Anno_Variants.txt")

variant_effects[variant_effects=="-"] <- NA

#Create columns for merging and merge
chrom <- variant_effects$X.Uploaded_variation
chrom <- gsub("_.*","",chrom)
pos <- variant_effects$X.Uploaded_variation
pos <- gsub("chr.._","",pos)
pos <- gsub("chr._","",pos)
pos <- gsub("_.*","",pos)
pos <- (as.numeric(pos))
variant_effects$merge_col <- paste(chrom,pos,sep="_")
variants$merge_col <- paste(variants$CHROM,variants$POS,sep="_")

variant_effects <- merge(variant_effects,variants,by="merge_col")

#Change first column name
colnames(variant_effects)[2] <- "Variant"

#Parse SpliceAI Delta score so that it can be used for sorting
variant_effects$SpliceAI_Delta <- sub(".*?\\|","",variant_effects$SpliceAI_pred)
variant_effects$SpliceAI_Delta <- gsub("\\|.*","",variant_effects$SpliceAI_Delta)
variant_effects$SpliceAI_Delta <- as.numeric(variant_effects$SpliceAI_Delta)

#Add single apostraphe to prevent excel from converting values to dates
variant_effects$EXON <- paste0("'",variant_effects$EXON)
variant_effects$INTRON <- paste0("'",variant_effects$INTRON)

#Write to a file
write.table(variant_effects,"RCM_WES_Apr18_Joint_SNVs_All_Genes_Final.txt",sep = "\t",quote=F,row.names = F,col.names = T)

#Number of unique variants 
length(unique(variant_effects$Variant))

#Get moderate impact variants
moderate_variants <- variant_effects[variant_effects$IMPACT=="MODERATE",]

#Number of moderate impact variants
length(unique(moderate_variants$Variant))

#Get high impact variants
high_variants <- variant_effects[variant_effects$IMPACT=="HIGH",]

#Number of high impact variants
length(unique(high_variants$Variant))

#write to files
write.table(moderate_variants,"RCM_WES_Apr18_Joint_SNVs_All_Genes_Moderate_Impact.txt",sep="\t",quote=F,row.names=F,col.names=T)

write.table(high_variants,"RCM_WES_Apr18_Joint_SNVs_All_High_Impact.txt",sep="\t",quote=F,row.names=F,col.names=T)

#Create a short list of unique variants for quick reference
moderate_variants <- moderate_variants[,c("Variant","Location","Allele","Gene","Consequence","Existing_variation","IMPACT","STRAND","VARIANT_CLASS","SYMBOL","HGNC_ID","MANE_SELECT","TSL","SWISSPROT","AF","gnomADe_AF","gnomADg_AF","MAX_AF","CLIN_SIG","CADD_PHRED","CADD_RAW","SpliceAI_pred","SpliceAI_Delta","PrimateAI","am_class","am_pathogenicity","SIFT","PolyPhen","IU.0105","IU.0106","P0142","P0144","P0177","P0178","P0197","P0199","P0203","P0215","P0220","P0223","P0233","P0247","P0256","P0288","P0328","P0339")]

high_variants <- high_variants[,c("Variant","Location","Allele","Gene","Consequence","Existing_variation","IMPACT","STRAND","VARIANT_CLASS","SYMBOL","HGNC_ID","MANE_SELECT","TSL","SWISSPROT","AF","gnomADe_AF","gnomADg_AF","MAX_AF","CLIN_SIG","CADD_PHRED","CADD_RAW","SpliceAI_pred","SpliceAI_Delta","PrimateAI","am_class","am_pathogenicity","SIFT","PolyPhen","IU.0105","IU.0106","P0142","P0144","P0177","P0178","P0197","P0199","P0203","P0215","P0220","P0223","P0233","P0247","P0256","P0288","P0328","P0339")]

#Remove redundancies
moderate_variants <- moderate_variants[!duplicated(moderate_variants$Variant),]
high_variants <- high_variants[!duplicated(high_variants$Variant),]

#Write to files
write.table(moderate_variants,"RCM_WES_Apr18_Joint_SNVs_All_Genes_Moderate_Impact_Short.txt",sep="\t",quote=F,row.names=F,col.names=T)

write.table(high_variants,"RCM_WES_Apr18_Joint_SNVs_All_Genes_High_Impact_Short.txt",sep="\t",quote=F,row.names=F,col.names=T)

#Parse Variants by proband####
INDELs <- read.delim("RCM_WES_Apr18_Joint_INDELs_All_Genes_Final.txt")
SNVs <- read.delim("RCM_WES_Apr18_Joint_SNVs_All_Genes_Final.txt")

Probands <- c("IU.0105","IU.0106","P0177","P0178","P0199","P0203","P0215","P0220","P0223","P0233","P0247","P0288","P0328","P0339")
Controls <- c("P0142","P0144","P0197","P0256")

for (i in 1:length(Probands)){
  #Get individual proband
  Proband <- Probands[i]
  #Get IDs for other probands
  Other_Probands <- Probands[Probands != Proband]
  #Get INDELs for individual proband
  INDELs_Proband <- INDELs[,!(colnames(INDELs)%in% Other_Probands) & !(colnames(INDELs) %in% Controls)]
  #Change column name containing genotype
  colnames(INDELs_Proband)[colnames(INDELs_Proband)==Proband] <- "Geno"
  #Keep values that aren't missing or ref/ref
  INDELs_Proband <- INDELs_Proband[INDELs_Proband$Geno != "0/0" & INDELs_Proband$Geno != "0|0" & INDELs_Proband$Geno != "./." & INDELs_Proband$Geno != ".|.",]
  #Keep INDELs with high or moderate impact
  INDELs_Proband <- INDELs_Proband[INDELs_Proband$IMPACT=="HIGH" | INDELs_Proband$IMPACT=="MODERATE",]
  #Add a column to indicate if the gene affected is in the list of candidate genes
  INDELs_Proband$Candidate_Gene <- FALSE
  INDELs_Proband$Candidate_Gene[INDELs_Proband$Gene %in% genes$V1] <- TRUE
  #Change Geno column name back to proband
  colnames(INDELs_Proband)[colnames(INDELs_Proband)=="Geno"] <- Proband
  #Keep some columns and remove redundancies
  cols <- c("Variant","Location","Allele","Gene","Candidate_Gene","Consequence","Existing_variation","IMPACT","STRAND","VARIANT_CLASS","SYMBOL","HGNC_ID","MANE_SELECT","TSL","SWISSPROT","AF","gnomADe_AF","gnomADg_AF","MAX_AF","CLIN_SIG","CADD_PHRED","CADD_RAW","SpliceAI_pred","SpliceAI_Delta",Proband)
  INDELs_Proband <- INDELs_Proband[,cols]
  INDELs_Proband <- INDELs_Proband[!duplicated(INDELs_Proband$Variant),]
  
  #Get SNVs for individual proband
  SNVs_Proband <- SNVs[,!(colnames(SNVs)%in% Other_Probands) & !(colnames(SNVs) %in% Controls)]
  #Change column name containing genotype
  colnames(SNVs_Proband)[colnames(SNVs_Proband)==Proband] <- "Geno"
  #Keep values that aren't missing or ref/ref
  SNVs_Proband <- SNVs_Proband[SNVs_Proband$Geno != "0/0" & SNVs_Proband$Geno != "0|0" & SNVs_Proband$Geno != "./." & SNVs_Proband$Geno != ".|.",]
  #Keep SNVs with high or moderate impact
  SNVs_Proband <- SNVs_Proband[SNVs_Proband$IMPACT=="HIGH" | SNVs_Proband$IMPACT=="MODERATE",]
  #Add a column to indicate if the gene affected is in the list of candidate genes
  SNVs_Proband$Candidate_Gene <- FALSE
  SNVs_Proband$Candidate_Gene[SNVs_Proband$Gene %in% genes$V1] <- TRUE
  #Change Geno column name back to proband
  colnames(SNVs_Proband)[colnames(SNVs_Proband)=="Geno"] <- Proband
  #Keep some columns and remove redundancies
  cols <- c("Variant","Location","Allele","Gene","Candidate_Gene","Consequence","Existing_variation","IMPACT","STRAND","VARIANT_CLASS","SYMBOL","HGNC_ID","MANE_SELECT","TSL","SWISSPROT","AF","gnomADe_AF","gnomADg_AF","MAX_AF","CLIN_SIG","CADD_PHRED","CADD_RAW","SpliceAI_pred","SpliceAI_Delta","PrimateAI","am_class","am_pathogenicity","SIFT","PolyPhen",Proband)
  SNVs_Proband <- SNVs_Proband[,cols]
  SNVs_Proband <- SNVs_Proband[!duplicated(SNVs_Proband$Variant),]
  
  #Write results to files
  write.table(INDELs_Proband,file=paste0(Proband,"_INDELs_All_Genes_Moderate_and_High_Impact.txt"),sep="\t",row.names=F,col.names=T,quote=F)
  write.table(SNVs_Proband,file=paste0(Proband,"_SNVs_All_Genes_Moderate_and_High_Impact.txt"),sep="\t",row.names=F,col.names=T,quote=F)
}
