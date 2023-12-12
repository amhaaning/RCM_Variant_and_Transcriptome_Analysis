library(ggplot2)
library(gridExtra)
library(dplyr)

setwd("/N/project/WareLab_ARP_NGS/RCM_VCFs/")

SNV_Qual <- read.delim("SNV_Qual_Metrics.txt",sep=" ")
Indel_Qual <- read.delim("Indel_Qual_Metrics.txt",sep=" ")

#Replace "." with NA
SNV_Qual[SNV_Qual=="."] <- NA
Indel_Qual[Indel_Qual=="."] <- NA

SNV_Qual$QD <- as.numeric(SNV_Qual$QD)
SNV_Qual$SOR <- as.numeric(SNV_Qual$SOR)
SNV_Qual$FS <- as.numeric(SNV_Qual$FS)
SNV_Qual$MQ <- as.numeric(SNV_Qual$MQ)
SNV_Qual$MQRankSum <- as.numeric(SNV_Qual$MQRankSum)
SNV_Qual$ReadPosRankSum <- as.numeric(SNV_Qual$ReadPosRankSum)

Indel_Qual$QD <- as.numeric(Indel_Qual$QD)
Indel_Qual$FS <- as.numeric(Indel_Qual$FS)
Indel_Qual$ReadPosRankSum <- as.numeric(Indel_Qual$ReadPosRankSum)

nrow(SNV_Qual[SNV_Qual$FS > 60 & !is.na(SNV_Qual$FS),])
nrow(SNV_Qual[SNV_Qual$MQ < 40 & !is.na(SNV_Qual$MQ),])
nrow(SNV_Qual[SNV_Qual$MQRankSum < -12.5 & !is.na(SNV_Qual$MQRankSum),])
nrow(SNV_Qual[SNV_Qual$QD < 2 & !is.na(SNV_Qual$QD),])
nrow(SNV_Qual[SNV_Qual$QUAL < 30 & !is.na(SNV_Qual$QUAL),])
nrow(SNV_Qual[SNV_Qual$ReadPosRankSum < -8 & !is.na(SNV_Qual$ReadPosRankSum),])
nrow(SNV_Qual[SNV_Qual$SOR > 3 & !is.na(SNV_Qual$SOR),])

#Make density plots for SNV qual scores
SNV_QD <- ggplot(SNV_Qual,aes(x=QD)) + geom_density() + geom_vline(aes(xintercept=2),color="blue",linetype="dashed",size=1) + xlim(0,100) + geom_rect(aes(xmin=0,xmax=2,ymin=0,ymax=Inf),fill="blue",alpha=0.005)
SNV_QUAL <- ggplot(SNV_Qual,aes(x=QUAL)) + geom_density() + geom_vline(aes(xintercept=30),color="blue",linetype="dashed",size=1) + xlim(0,5000) + geom_rect(aes(xmin=0,xmax=30,ymin=0,ymax=Inf),fill="blue",alpha=0.005)
SNV_SOR <- ggplot(SNV_Qual,aes(x=SOR)) + geom_density() + geom_vline(aes(xintercept=3),color="blue",linetype="dashed",size=1) + xlim(0,10) + geom_rect(aes(xmin=3,xmax=Inf,ymin=0,ymax=Inf),fill="blue",alpha=0.005)
SNV_FS <- ggplot(SNV_Qual,aes(x=FS)) + geom_density() + geom_vline(aes(xintercept=60),color="blue",linetype="dashed",size=1) + xlim(0,100) + geom_rect(aes(xmin=60,xmax=Inf,ymin=0,ymax=Inf),fill="blue",alpha=0.005)
SNV_MQ <- ggplot(SNV_Qual,aes(x=MQ)) + geom_density() + geom_vline(aes(xintercept=40),color="blue",linetype="dashed",size=1) + xlim(0,100) + geom_rect(aes(xmin=0,xmax=40,ymin=0,ymax=Inf),fill="blue",alpha=0.005)
SNV_MQRS <- ggplot(SNV_Qual,aes(x=MQRankSum)) + geom_density() + geom_vline(aes(xintercept=-12.5),color="blue",linetype="dashed",size=1) + xlim(-13,13) + geom_rect(aes(xmin=-13,xmax=-12.5,ymin=0,ymax=Inf),fill="blue",alpha=0.005)
SNV_RPRS <- ggplot(SNV_Qual,aes(x=ReadPosRankSum)) + geom_density() + geom_vline(aes(xintercept=-8),color="blue",linetype="dashed",size=1) + xlim(-10,10) + geom_rect(aes(xmin=-10,xmax=-8,ymin=0,ymax=Inf),fill="blue",alpha=0.005)
SNV_SOR <- ggplot(SNV_Qual,aes(x=SOR)) + geom_density() + geom_vline(aes(xintercept=3),color="blue",linetype="dashed",size=1) + xlim(0,10) + geom_rect(aes(xmin=3,xmax=Inf,ymin=0,ymax=Inf),fill="blue",alpha=0.005)

pdf("SNV_Qual_Density_Plots.pdf",width=6,height=8)
grid.arrange(SNV_FS,SNV_MQ,SNV_MQRS,SNV_QD,SNV_QUAL,SNV_RPRS,SNV_SOR,nrow=4)
dev.off()

nrow(Indel_Qual[Indel_Qual$FS > 60 & !is.na(Indel_Qual$FS),])
nrow(Indel_Qual[Indel_Qual$QD < 2 & !is.na(Indel_Qual$QD),])
nrow(Indel_Qual[Indel_Qual$QUAL < 30 & !is.na(Indel_Qual$QUAL),])
nrow(Indel_Qual[Indel_Qual$ReadPosRankSum < -20 & !is.na(Indel_Qual$ReadPosRankSum),])

#Make density plots for Indel qual scores
Indel_QD <- ggplot(Indel_Qual,aes(x=QD)) + geom_density() + geom_vline(aes(xintercept=2),color="blue",linetype="dashed",size=1) + xlim(0,100) + geom_rect(aes(xmin=0,xmax=2,ymin=0,ymax=Inf),fill="blue",alpha=0.005)
Indel_QUAL <- ggplot(Indel_Qual,aes(x=QUAL)) + geom_density() + geom_vline(aes(xintercept=30),color="blue",linetype="dashed",size=1) + xlim(0,5000) + geom_rect(aes(xmin=0,xmax=30,ymin=0,ymax=Inf),fill="blue",alpha=0.005)
Indel_FS <- ggplot(Indel_Qual,aes(x=FS)) + geom_density() + geom_vline(aes(xintercept=60),color="blue",linetype="dashed",size=1) + xlim(0,100) + geom_rect(aes(xmin=60,xmax=Inf,ymin=0,ymax=Inf),fill="blue",alpha=0.005)
Indel_RPRS <- ggplot(Indel_Qual,aes(x=ReadPosRankSum)) + geom_density() + geom_vline(aes(xintercept=-20),color="blue",linetype="dashed",size=1) + xlim(-30,10) + geom_rect(aes(xmin=-30,xmax=-20,ymin=0,ymax=Inf),fill="blue",alpha=0.005)

pdf("Indel_Qual_Density_Plots.pdf",width=6,height=4)
grid.arrange(Indel_FS,Indel_QD,Indel_QUAL,Indel_RPRS,nrow=2)
dev.off()

SNV_Qual <- read.delim("SNVs_Individual_Variant_Quality_Parameters.txt",sep=" ",header=F)
Indel_Qual <- read.delim("Indels_Individual_Variant_Quality_Parameters.txt",sep=" ",header=F)

SNV_GQ <- SNV_Qual[,c(2:19)]
SNV_GQ[SNV_GQ=="."] <- NA
SNV_GQ <- mutate_all(SNV_GQ, function(x) as.numeric(x))
SNV_GQ_Mean <- data.frame(rowMeans(SNV_GQ,T))
colnames(SNV_GQ_Mean) <- "SNV_GQ_Mean"

#Density plot of mean GQ (genotyping quality) for all variants
SNV_GQ_plot <- ggplot(SNV_GQ_Mean,aes(x=SNV_GQ_Mean)) + geom_density() + geom_vline(aes(xintercept=40),color="blue",linetype="dashed",size=1) + geom_rect(aes(xmin=0,xmax=40,ymin=0,ymax=Inf),fill="blue",alpha=0.005)

length(SNV_GQ_Mean[SNV_GQ_Mean$SNV_GQ_Mean < 40 & !is.na(SNV_GQ_Mean$SNV_GQ_Mean),])

#Get all GT scores for all variants
SNV_GT <- SNV_Qual[,c(22:39)]

#Calculate the number of missing genotype calls for each variant
SNV_GT$Number_Missing <- apply(SNV_GT, 1, function(x) length(which(x=="./.")))

SNV_Missing <- ggplot(SNV_GT,aes(x=Number_Missing)) + geom_density() + geom_vline(aes(xintercept=4),color="blue",linetype="dashed",size=1) + geom_rect(aes(xmin=4,xmax=18,ymin=0,ymax=Inf),fill="blue",alpha=0.005) + xlab("Number of missing calls")

#Number of SNVs missing more than 75% of calls
nrow(SNV_GT[SNV_GT$Number_Missing < 1,])

pdf("Individual_SNV_Qual_Density_Plots.pdf",width=6,height=2)
grid.arrange(SNV_GQ_plot,SNV_Missing,nrow=1)
dev.off()
