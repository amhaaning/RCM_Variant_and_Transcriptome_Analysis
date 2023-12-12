#!/bin/bash -l 
#SBATCH --time=24:00:00 
#SBATCH --ntasks=1 
#SBATCH --cpus-per-task=1 
#SBATCH --mem=10gb
#SBATCH -A r00126
#SBATCH --mail-user=ahaaning@iu.edu 
#SBATCH --mail-type=ALL

cd /N/project/WareLab_ARP_NGS/RCM_VCFs

vcf_in=RCM_WES_Apr18_Joint_INDELs.vcf
vcf_out=RCM_WES_Apr18_Joint_INDELs_VariantFiltration.vcf

#Count the number of variants in the unfiltered file
echo $vcf_in" Number of Variants:"
gatk CountVariants \
     -V $vcf_in

Ref=/N/project/Ware-lab_NGS/References/gatk_resource_bundle/Homo_sapiens_assembly38.fasta

#Filter variants based on cutoffs for parameters recommended by GATK
#Also set genotype to missing if it doesn't pass the genotype filter
gatk VariantFiltration \
     -R $Ref -V $vcf_in \
     -O $vcf_out \
     --filter-expression "FS > 200.0" --filter-name "highFS"\
     --filter-expression "QD < 2.0" --filter-name "lowQD" \
     --filter-expression "ReadPosRankSum < -20.0" --filter-name "lowRPRS" \
     --genotype-filter-expression "GQ < 25.0" --genotype-filter-name "lowGQ" \
     --genotype-filter-expression "DP < 8.0" --genotype-filter-name "lowDP"

echo $vcf_out" Number of Variants:"
gatk CountVariants \
      -V $vcf_out

vcf_in=RCM_WES_Apr18_Joint_INDELs_VariantFiltration.vcf
vcf_out=RCM_WES_Apr18_Joint_INDELs_SelectVariants.vcf

#Remove variants based on filters above or set to missing based on
#individual genotype filters
gatk SelectVariants \
     -R $Ref \
     -V $vcf_in \
     -O $vcf_out \
     --set-filtered-gt-to-nocall true \
     --exclude-filtered true \
     --max-nocall-number 4 \
     --remove-unused-alternates true
     
echo $vcf_out" Number of Variants:"
gatk CountVariants \
      -V $vcf_out

vcf_in=RCM_WES_Apr18_Joint_INDELs_SelectVariants.vcf
vcf_out=RCM_WES_Apr18_Joint_INDELs_SelectVariants_Final.vcf

#Remove variants if they have more than 4 missing calls
#It's possible that this won't remove any variants
#I'm not sure of the order of operations with the above filters
#For example, does it filter out variants with too many missing
#before it sets them to missing based on gt filters?
gatk SelectVariants \
     -R $Ref -V $vcf_in \
     -O $vcf_out \
     --max-nocall-number 4

echo $vcf_out" Number of Variants:"
gatk CountVariants \
      -V $vcf_out
