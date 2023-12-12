library(ggplot2)
library(gridExtra)
library(dplyr)
library(stringr)

setwd("C:/Users/ahaaning/OneDrive - Indiana University/Manuscripts/RCM_Variants_and_RNAseq/New_Variants/")

#SNV metrics distributions - variant level metrics####

#Read in SNV qualities
SNV_Qual <- read.delim("RCM_WES_Apr18_Joint_QC_Metrics_SNVs.txt",header=F)
colnames(SNV_Qual) <- c("CHROM","POS","REF","ALT","FS","MQ","MQRankSum","QD","ReadPosRankSum","SOR")

#Replace "." with NA
SNV_Qual[SNV_Qual=="."] <- NA

#Make sure metrics are numeric
SNV_Qual$QD <- as.numeric(SNV_Qual$QD)
SNV_Qual$SOR <- as.numeric(SNV_Qual$SOR)
SNV_Qual$FS <- as.numeric(SNV_Qual$FS)
SNV_Qual$MQ <- as.numeric(SNV_Qual$MQ)
SNV_Qual$MQRankSum <- as.numeric(SNV_Qual$MQRankSum)
SNV_Qual$ReadPosRankSum <- as.numeric(SNV_Qual$ReadPosRankSum)

#Count the number of variants that would be removed if GATK starting thresholds were used
nrow(SNV_Qual[SNV_Qual$FS > 60 & !is.na(SNV_Qual$FS),])
nrow(SNV_Qual[SNV_Qual$MQ < 40 & !is.na(SNV_Qual$MQ),])
nrow(SNV_Qual[SNV_Qual$MQRankSum < -12.5 & !is.na(SNV_Qual$MQRankSum),])
nrow(SNV_Qual[SNV_Qual$QD < 2 & !is.na(SNV_Qual$QD),])
nrow(SNV_Qual[SNV_Qual$ReadPosRankSum < -8 & !is.na(SNV_Qual$ReadPosRankSum),])
nrow(SNV_Qual[SNV_Qual$SOR > 3 & !is.na(SNV_Qual$SOR),])

#Make density plots for SNV qual scores
tiff("SNV_QD_Density.tiff",height=2,width=3,units = "in",res=600)
ggplot(SNV_Qual,aes(x=QD)) + geom_density() + geom_vline(aes(xintercept=2),color="blue",linetype="dashed",linewidth=1) + xlim(0,100) + geom_rect(aes(xmin=0,xmax=2,ymin=0,ymax=Inf),fill="blue",alpha=0.005) + theme_bw()
dev.off()

tiff("SNV_SOR_Density.tiff",height=2,width=3,units = "in",res=600)
ggplot(SNV_Qual,aes(x=SOR)) + geom_density() + geom_vline(aes(xintercept=3),color="blue",linetype="dashed",linewidth=1) + xlim(0,10) + geom_rect(aes(xmin=3,xmax=Inf,ymin=0,ymax=Inf),fill="blue",alpha=0.005) + theme_bw()
dev.off()

tiff("SNV_FS_Density.tiff",height=2,width=3,units = "in",res=600)
ggplot(SNV_Qual,aes(x=FS)) + geom_density() + geom_vline(aes(xintercept=60),color="blue",linetype="dashed",linewidth=1) + xlim(0,100) + geom_rect(aes(xmin=60,xmax=Inf,ymin=0,ymax=Inf),fill="blue",alpha=0.005) + theme_bw()
dev.off()

tiff("SNV_MQ_Density.tiff",height=2,width=3,units = "in",res=600)
ggplot(SNV_Qual,aes(x=MQ)) + geom_density() + geom_vline(aes(xintercept=40),color="blue",linetype="dashed",linewidth=1) + xlim(0,100) + geom_rect(aes(xmin=0,xmax=40,ymin=0,ymax=Inf),fill="blue",alpha=0.005) + theme_bw()
dev.off()

tiff("SNV_MQRankSum_Density.tiff",height=2,width=3,units = "in",res=600)
ggplot(SNV_Qual,aes(x=MQRankSum)) + geom_density() + geom_vline(aes(xintercept=-12.5),color="blue",linetype="dashed",linewidth=1) + xlim(-13,13) + geom_rect(aes(xmin=-13,xmax=-12.5,ymin=0,ymax=Inf),fill="blue",alpha=0.005) + theme_bw()
dev.off()

tiff("SNV_ReadPosRankSum_Density.tiff",height=2,width=3,units = "in",res=600)
ggplot(SNV_Qual,aes(x=ReadPosRankSum)) + geom_density() + geom_vline(aes(xintercept=-8),color="blue",linetype="dashed",linewidth=1) + xlim(-10,10) + geom_rect(aes(xmin=-10,xmax=-8,ymin=0,ymax=Inf),fill="blue",alpha=0.005) + theme_bw()
dev.off()

tiff("SNV_SOR_Density.tiff",height=2,width=3,units = "in",res=600)
ggplot(SNV_Qual,aes(x=SOR)) + geom_density() + geom_vline(aes(xintercept=3),color="blue",linetype="dashed",linewidth=1) + xlim(0,10) + geom_rect(aes(xmin=3,xmax=Inf,ymin=0,ymax=Inf),fill="blue",alpha=0.005) + theme_bw()
dev.off()

#INDEL metrics distributions - variant level metrics####

#Read in INDEL qualities
INDEL_Qual <- read.delim("RCM_WES_Apr18_Joint_QC_Metrics_INDELs.txt",header=F)
colnames(INDEL_Qual) <- c("CHROM","POS","REF","ALT","FS","MQ","MQRankSum","QD","ReadPosRankSum","SOR")

#Replace "." with NA
INDEL_Qual[INDEL_Qual=="."] <- NA

#Make sure metrics are numeric
INDEL_Qual$QD <- as.numeric(INDEL_Qual$QD)
INDEL_Qual$FS <- as.numeric(INDEL_Qual$FS)
INDEL_Qual$ReadPosRankSum <- as.numeric(INDEL_Qual$ReadPosRankSum)

#Count the number of variants that would be removed if GATK starting thresholds were used
nrow(INDEL_Qual[INDEL_Qual$FS > 60 & !is.na(INDEL_Qual$FS),])
nrow(INDEL_Qual[INDEL_Qual$QD < 2 & !is.na(INDEL_Qual$QD),])
nrow(INDEL_Qual[INDEL_Qual$ReadPosRankSum < -20 & !is.na(INDEL_Qual$ReadPosRankSum),])

#Make density plots for INDEL qual scores
tiff("INDEL_QD_Density.tiff",height=2,width=3,units = "in",res=600)
ggplot(INDEL_Qual,aes(x=QD)) + geom_density() + geom_vline(aes(xintercept=2),color="blue",linetype="dashed",size=1) + xlim(0,100) + geom_rect(aes(xmin=0,xmax=2,ymin=0,ymax=Inf),fill="blue",alpha=0.005) + theme_bw()
dev.off()

tiff("INDEL_FS.tiff",height=2,width=3,units = "in",res=600)
ggplot(INDEL_Qual,aes(x=FS)) + geom_density() + geom_vline(aes(xintercept=60),color="blue",linetype="dashed",size=1) + xlim(0,100) + geom_rect(aes(xmin=60,xmax=Inf,ymin=0,ymax=Inf),fill="blue",alpha=0.005) + theme_bw()
dev.off()

tiff("INDEL_ReadPosRankSum_Density.tiff",height=2,width=3,units = "in",res=600)
ggplot(INDEL_Qual,aes(x=ReadPosRankSum)) + geom_density() + geom_vline(aes(xintercept=-20),color="blue",linetype="dashed",size=1) + xlim(-30,10) + geom_rect(aes(xmin=-30,xmax=-20,ymin=0,ymax=Inf),fill="blue",alpha=0.005) + theme_bw()
dev.off()

#SNV metrics distributions - sample level metrics####

#Read in SNV qualities
SNV_Qual_Ind <- read.delim("RCM_WES_Apr18_Joint_QC_Metrics_Sample_Level_SNVs.txt",header=F)

#Get genotype quality for all samples and convert to numeric
SNV_GQ <- SNV_Qual_Ind[,c(22:39)]
SNV_GQ <- mutate_all(SNV_GQ, function(x) as.numeric(x))

#Get the number of reads supporting ref and alt alleles
SNV_Ref_Alt_Reads <- SNV_Qual_Ind[,c(40:57)]

#Get the number of total reads for all samples and convert to numeric
SNV_Total_Reads <- mutate_all(SNV_Qual_Ind[,c(58:75)], function(x) as.numeric(x))

#Get individual genotypes and assign categories
SNV_GT <- SNV_Qual_Ind[,c(4:21)]

#Convert genotypes, alt read ratios, and genotype qualities to vectors and combine in dataframe
SNV_Qual_Ind <- data.frame(cbind(as.vector(as.matrix(SNV_GT)),as.vector(as.matrix(SNV_Total_Reads)),as.vector(as.matrix(SNV_Ref_Alt_Reads)),as.vector(as.matrix(SNV_GQ))))

#Set column names
colnames(SNV_Qual_Ind) <- c("GT","Total_Reads","Ref_Alt_Reads","GQ")

#Get the total number of possible alleles and remove multiallelics
SNV_Qual_Ind$Number_of_Alleles <- (str_count(SNV_Qual_Ind$Ref_Alt_Reads,",")) + 1
SNV_Qual_Ind <- SNV_Qual_Ind[SNV_Qual_Ind$Number_of_Alleles < 3,]

#Get the number of ref and alt alleles
SNV_Qual_Ind$Ref <- gsub(",.*","",SNV_Qual_Ind$Ref_Alt_Reads)
SNV_Qual_Ind$Alt <- gsub(".*,","",SNV_Qual_Ind$Ref_Alt_Reads)

#Make sure values are numeric
SNV_Qual_Ind$GQ <- as.numeric(SNV_Qual_Ind$GQ)
SNV_Qual_Ind$Ref <- as.numeric(SNV_Qual_Ind$Ref)
SNV_Qual_Ind$Alt <- as.numeric(SNV_Qual_Ind$Alt)
SNV_Qual_Ind$Alt_Read_Ratio <- SNV_Qual_Ind$Alt/(SNV_Qual_Ind$Alt + SNV_Qual_Ind$Ref)

#Remove rows with NA values
SNV_Qual_Ind <- SNV_Qual_Ind[complete.cases(SNV_Qual_Ind),]

#Recalculate total reads based on number of ref and alt alleles
SNV_Qual_Ind$Total_Reads <- SNV_Qual_Ind$Alt + SNV_Qual_Ind$Ref

#Distribution of total reads
tiff("SNV_Total_Reads_Distribution.tiff",height=2,width=3,units = "in",res=600)
ggplot(SNV_Qual_Ind,aes(x=Total_Reads)) + geom_density() + theme_bw() + xlab("Read Depth of Sample Calls")
dev.off()

#Mean read depth
mean(SNV_Qual_Ind$Total_Reads[!is.na(SNV_Qual_Ind$Total_Reads)])

#Distribution of genotype quality
tiff("SNV_GQ_Distribution.tiff",height=2,width=3,units = "in",res=600)
ggplot(SNV_Qual_Ind,aes(x=GQ)) + geom_density() + theme_bw() + xlab("Genotype Quality (GQ)")
dev.off()

#Remove rows with fewer than 8 total reads
SNV_Qual_Ind <- SNV_Qual_Ind[SNV_Qual_Ind$Total_Reads >= 5,]

#Individual genotype scores for hets####

#Pull data for hets
Het_SNV_Qual_Ind <- SNV_Qual_Ind[SNV_Qual_Ind$GT=="0/1" | SNV_Qual_Ind$GT=="0|1",]

#Make density plot for alt read ratio in hets
tiff("SNV_Het_Alt_Read_Ratio.tiff",height=2,width=3,units = "in",res=600)
ggplot(Het_SNV_Qual_Ind,aes(x=Alt_Read_Ratio)) + geom_density() + theme_bw() + xlab("Alt Read Ratio - Hets")
dev.off()

#Look at relationship between alt read ratio and genotype quality in hets
tiff("SNV_Hets_Alt_Read_Ratio_vs_GQ.tiff",height=2,width=3,units = "in",res=600)
ggplot(Het_SNV_Qual_Ind,aes(x=Alt_Read_Ratio,y=GQ)) + geom_point() + theme_bw() + xlab("Alt Read Ratio") + ylab("Genotype Quality")
dev.off()

#Individual genotype scores for hom ref####

#Pull data for hom refs
Hom_Ref_SNV_Qual_Ind <- SNV_Qual_Ind[SNV_Qual_Ind$GT=="0/0" | SNV_Qual_Ind$GT=="0|0",]

#Make density plot for alt read ratio in hom refs
tiff("SNV_Hom_Ref_Alt_Read_Ratio.tiff",height=2,width=3,units = "in",res=600)
ggplot(Hom_Ref_SNV_Qual_Ind,aes(x=Alt_Read_Ratio)) + geom_density() + theme_bw() + xlab("Alt Read Ratio - Hom Ref")
dev.off()

#Look at relationship between alt read ratio and genotype quality in hom refs
tiff("SNV_Hom_Ref_Alt_Read_Ratio_vs_GQ.tiff",height=2,width=3,units = "in",res=600)
ggplot(Hom_Ref_SNV_Qual_Ind,aes(x=Alt_Read_Ratio,y=GQ)) + geom_point() + theme_bw() + xlab("Alt Read Ratio") + ylab("Genotype Quality")
dev.off()

#Individual genotype scores for hom alt####

#Pull data for hom alts
Hom_Alt_SNV_Qual_Ind <- SNV_Qual_Ind[SNV_Qual_Ind$GT=="1/1" | SNV_Qual_Ind$GT=="1|1",]

#Make density plot for alt read ratio in hom alts
tiff("SNV_Hom_Alt_Alt_Read_Ratio.tiff",height=2,width=3,units = "in",res=600)
ggplot(Hom_Alt_SNV_Qual_Ind,aes(x=Alt_Read_Ratio)) + geom_density() + theme_bw() + xlab("Alt Read Ratio - Hom Alt")
dev.off()

#Look at relationship between alt read ratio and genotype quality in hom alts
tiff("SNV_Hom_Alt_Alt_Read_Ratio_vs_GQ.tiff",height=2,width=3,units = "in",res=600)
ggplot(Hom_Alt_SNV_Qual_Ind,aes(x=Alt_Read_Ratio,y=GQ)) + geom_point() + theme_bw() + xlab("Alt Read Ratio") + ylab("Genotype Quality")
dev.off()

#Look at relationship between total reads and genotype quality

#Calculate the number of missing genotype calls for each variant
SNV_GT$Number_Missing <- apply(SNV_GT, 1, function(x) length(which(x=="./.")) + length(which(x==".|.")))

tiff("SNV_Missing_Distribution.tiff",height=2,width=3,units = "in",res=600)
ggplot(SNV_GT,aes(x=Number_Missing)) + geom_density() + xlab("Number of missing calls/variant") + theme_bw()
dev.off()

#Number of SNVs missing more than 75% of calls
nrow(SNV_GT[SNV_GT$Number_Missing > 4,])

#Number of missing calls 
sum(SNV_GT$Number_Missing)
max(SNV_GT$Number_Missing)
