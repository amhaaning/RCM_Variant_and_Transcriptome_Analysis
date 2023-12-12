#!/bin/bash -l 
#SBATCH --time=24:00:00 
#SBATCH --ntasks=1 
#SBATCH --cpus-per-task=5 
#SBATCH --mem=50gb
#SBATCH -A r00126
#SBATCH --mail-user=ahaaning@iu.edu 
#SBATCH --mail-type=ALL

cd /N/project/WareLab_ARP_NGS/RCM_VCFs

vcf_in=RCM_WES_Apr18_Joint.vcf.gz
vcf_SNV_out=RCM_WES_Apr18_Joint_SNVs.vcf
vcf_INDELS_out=RCM_WES_Apr18_Joint_INDELs.vcf

#Separate variants into SNVs and Indels
gatk --java-options "-Xms50g -Xmx50g -XX:ParallelGCThreads=5" SelectVariants -V $vcf_in -select-type SNP -O $vcf_SNV_out
gatk --java-options "-Xms50g -Xmx50g -XX:ParallelGCThreads=5" SelectVariants -V $vcf_in -select-type INDEL -O $vcf_INDELS_out

vcf_in=RCM_WES_Apr18_Joint_INDELs.vcf
vcf_out=RCM_WES_Apr18_Joint_INDELs_Norm.vcf

#Normalize indels
gatk --java-options "-Xms50g -Xmx50g -XX:ParallelGCThreads=5" LeftAlignAndTrimVariants -R /N/project/Ware-lab_NGS/References/gatk_resource_bundle/Homo_sapiens_assembly38.fasta -V $vcf_in -O $vcf_out --max-indel-length 1000 --max-leading-bases 2000


