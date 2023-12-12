#!/bin/bash -l
#SBATCH --time=12:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=40gb
#SBATCH -A r00126
#SBATCH --mail-user=ahaaning@iu.edu
#SBATCH --mail-type=ALL

genomicsdb=./GVCF/gatk_genomicsdb/

#Output directory
outdir=./RCM_VCFs/

#Create vcf output file name
outfile=$outdir"RCM_WES_Apr18_Joint.vcf.gz"

#Reference genome
Ref=./References/gatk_resource_bundle/Homo_sapiens_assembly38.fasta

gatk --java-options "-Xms35g -Xmx35g -XX:ParallelGCThreads=4" GenotypeGVCFs -R $Ref -V gendb://$genomicsdb -O $outfile --allow-old-rms-mapping-quality-annotation-data true
