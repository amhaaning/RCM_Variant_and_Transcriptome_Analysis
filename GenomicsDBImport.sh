#!/bin/bash -l 
#SBATCH --time=72:00:00 
#SBATCH --ntasks=1 
#SBATCH --cpus-per-task=20 
#SBATCH --mem=400gb
#SBATCH -A r00126
#SBATCH --mail-user=ahaaning@iu.edu 
#SBATCH --mail-type=ALL
#SBATCH --partition=largememory

ulimit -c unlimited

#Change to directory containing GVCFs
cd /N/project/WareLab_ARP_NGS/ILMN_249_Ware_WESseq18_Apr2018_GVCF/

#Path to cohort sample map
samples=/N/project/WareLab_ARP_NGS/RCM_VCFs/scripts/cohort.sample_map

#Path to genomicsdb workspace
genomicsdb=/N/project/WareLab_ARP_NGS/ILMN_249_Ware_WESseq18_Apr2018_GVCF/gatk_genomicsdb

#Path to file containing intervals
intervals=/N/project/WareLab_ARP_NGS/WESeq18_Apr2018_Bait_and_Target_Files/S07604514_Regions_target.interval_list

gatk --java-options "-Xms400g -Xmx400g -XX:ParallelGCThreads=20" GenomicsDBImport \
     --genomicsdb-workspace-path $genomicsdb \
     --interval-padding 100 -L $intervals \
     --sample-name-map $samples --tmp-dir /N/scratch/ahaaning/tmp_genomicsdb_RCM \
     --overwrite-existing-genomicsdb-workspace --reader-threads 20 \
     --genomicsdb-shared-posixfs-optimizations true \
     --bypass-feature-reader \
     --merge-input-intervals true
