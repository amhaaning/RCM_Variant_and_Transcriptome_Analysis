#!/bin/bash -l 
#SBATCH --time=24:00:00 
#SBATCH --ntasks=1 
#SBATCH --cpus-per-task=1
#SBATCH --mem=50gb
#SBATCH -A r00126
#SBATCH --mail-user=ahaaning@iu.edu 
#SBATCH --mail-type=ALL

module load perl

cd .

vcf_in=RCM_WES_Apr18_Joint_SNVs_SelectVariants_Final.vcf
vcf_out=RCM_WES_Apr18_Joint_SNVs_SelectVariants_Final_Anno.vcf

#Use vep to annotate quality filtered VCF
vep --cache --dir_cache ./Software_Packages/ensembl_vep_cache \
    -i $vcf_in --vcf -o $vcf_out --everything

vcf_in=RCM_WES_Apr18_Joint_SNVs_SelectVariants_Final_Anno.vcf
vcf_out=RCM_WES_Apr18_Joint_SNVs_SelectVariants_Final_Anno_Filt_Max_AF.vcf

#Filter VCF based on population allele frequencies
filter_vep -i $vcf_in \
	   -f "((gnomADe_AF <= 0.001 or not gnomADe_AF) and (gnomADg_AF <= 0.001 or not gnomADg_AF))" \
	   -o $vcf_out

#Use SnpSift to count the number of hom alt and hets in cases (+) and controls (-)
vcf_in=RCM_WES_Apr18_Joint_SNVs_SelectVariants_Final_Anno_Filt_Max_AF.vcf
vcf_out=RCM_WES_Apr18_Joint_SNVs_SelectVariants_Final_Filt_Anno_caseControl.vcf

java -jar /N/project/Ware-lab_NGS/Software_Packages/snpEff/SnpSift.jar \
    caseControl "++--++-+++++++-+++" \
    $vcf_in > $vcf_out

#Remove varaints that are hom alt or het in controls
vcf_in=RCM_WES_Apr18_Joint_SNVs_SelectVariants_Final_Filt_Anno_caseControl.vcf
vcf_out=RCM_WES_Apr18_Joint_SNVs_SelectVariants_Final_Filt_Anno_caseControl_Filt.vcf

java -jar ./Software_Packages/snpEff/SnpSift.jar filter \
    "Controls[2]=0" $vcf_in > $vcf_out
    
#Run vep again on final variants
vcf_in=RCM_WES_Apr18_Joint_SNVs_SelectVariants_Final_Filt_Anno_caseControl_Filt.vcf

vep --cache --dir_cache ./Software_Packages/ensembl_vep_cache \
    -i $vcf_in \
    --dir_plugins ./Software_Packages/ensembl_vep_plugins/ \
    --plugin CADD,snvs=./Software_Packages/ensembl_vep_plugins/whole_genome_SNVs_inclAnno.tsv.gz \
    --plugin AlphaMissense,file=./Software_Packages/ensembl_vep_plugins/AlphaMissense_hg38.tsv.gz \
    --plugin SpliceAI,snv=./Software_Packages/ensembl_vep_plugins/spliceai_scores.raw.snv.hg38.vcf.gz,indel=./Software_Packages/ensembl_vep_plugins/spliceai_scores.raw.indel.hg38.vcf.gz \
    --plugin PrimateAI,./Software_Packages/ensembl_vep_plugins/PrimateAI_scores_v0.2_GRCh38_sorted.tsv.bgz \
    -o RCM_WES_Apr18_Joint_SNVs_FINAL_Anno.vcf --force --everything --vcf --buffer_size 1000

vep --cache --dir_cache ./Software_Packages/ensembl_vep_cache \
    -i $vcf_in \
    --dir_plugins ./Software_Packages/ensembl_vep_plugins/ \
    --plugin CADD,snvs=./Software_Packages/ensembl_vep_plugins/whole_genome_SNVs_inclAnno.tsv.gz \
    --plugin AlphaMissense,file=./Software_Packages/ensembl_vep_plugins/AlphaMissense_hg38.tsv.gz \
    --plugin SpliceAI,snv=./Software_Packages/ensembl_vep_plugins/spliceai_scores.raw.snv.hg38.vcf.gz,indel=./Software_Packages/ensembl_vep_plugins/spliceai_scores.raw.indel.hg38.vcf.gz \
    --plugin PrimateAI,./Software_Packages/ensembl_vep_plugins/PrimateAI_scores_v0.2_GRCh38_sorted.tsv.bgz \
    -o RCM_WES_Apr18_Joint_SNVs_FINAL_Anno.txt --force --tab --everything --buffer_size 1000

module load vcftools

vcftools --vcf RCM_WES_Apr18_Joint_SNVs_FINAL_Anno.vcf --extract-FORMAT-info GT
mv out.GT.FORMAT RCM_WES_Apr18_Joint_SNVs_FINAL_Anno_Variants.txt
