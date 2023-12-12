#!/bin/bash -l
#SBATCH --time=12:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=40gb
#SBATCH -A r00126
#SBATCH --mail-user=ahaaning@iu.edu
#SBATCH --mail-type=ALL

genomicsdb=/N/project/WareLab_ARP_NGS/ILMN_249_Ware_WESseq18_Apr2018_GVCF/gatk_genomicsdb/

#Output directory
outdir=/N/project/WareLab_ARP_NGS/RCM_VCFs/

#Create vcf output file name
outfile=$outdir"RCM_WES_Apr18_Joint.vcf.gz"

#Reference genome
Ref=/N/project/Ware-lab_NGS/References/gatk_resource_bundle/Homo_sapiens_assembly38.fasta

gatk --java-options "-Xms35g -Xmx35g -XX:ParallelGCThreads=4" GenotypeGVCFs -R $Ref -V gendb://$genomicsdb -O $outfile --allow-old-rms-mapping-quality-annotation-data true
