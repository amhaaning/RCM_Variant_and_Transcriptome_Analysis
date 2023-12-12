#!/bin/bash -l 
#SBATCH --time=12:00:00 
#SBATCH --ntasks=1 
#SBATCH --cpus-per-task=10 
#SBATCH --mem=100gb
#SBATCH -A r00126
#SBATCH --mail-user=ahaaning@iu.edu 
#SBATCH --mail-type=ALL

#Directory containing recalibrated, mapped bam files for all individual samples
BAM_Dir=./ILMN_249_Ware_WESeq18_Apr2018_BAM/

#Directory for output GVCFs
GVCF_Dir=./GVCF/

#Populate array with directory info
IFS=$'\n' read -d '' -r -a samples < $BAM_Dir"sample_list.txt"

#Get individual sample ID
sample=${samples[${SLURM_ARRAY_TASK_ID}]}

#BAM file input for HaplotypeCaller
BAM=$BAM_Dir$sample".UMI_consensus.hg38.mapped.recal.bam"

#GVCF file output for HaplotypeCaller
GVCF=$GVCF_Dir$sample".g.vcf.gz"

#Run HaplotypeCaller to generate GVCFs from Bam files
gatk --java-options "-Xms100g -Xmx100g -XX:ParallelGCThreads=10" HaplotypeCaller  \
     -R ./References/gatk_resource_bundle/Homo_sapiens_assembly38.fasta \
     -I $BAM \
     -O $GVCF \
     -L ./WESeq18_Apr2018_Bait_and_Target_Files/S07604514_Regions_target.interval_list \
     -ip 100 \
     -ERC GVCF


