#!/bin/sh -l
#SBATCH -A fnrtowhee
#SBATCH -N 1
#SBATCH -n 10
#SBATCH -t 14-00:00:00
#SBATCH --job-name=norelMAFsyn
#SBATCH -e %x
#SBATCH -o %x

##########################################################################
###                             Samarth Mathur                          ###
###                        PhD Candidate, DeWoody Lab                   ###
###                           Purdue University                         ###
###                                                                     ###
###     Date Created: 03/23/22                  Last Modified: 03/23/22 ###
###########################################################################
###########################################################################
###                     norelMAF.sh                   		        	###
###########################################################################

cd $SLURM_SUBMIT_DIR
module purge
module load bioinfo
module load vcftools


# Recalculating the MAF by removing related WTX samples
#For each combination of the 24 unrelated WTX samples, we randomly sampled 24 individuals 
# from AZ (with replacement) and compared mean heterozygosity and DAF between the two populations using Student’s t-test


# To create null distribution

#To test for significance for the observed statistic (t*), we first created a null distribution of the t-statistic by 
# randomly choosing individuals from the pool of combined samples (N = 48)

cd /scratch/bell/mathur20/MQU/ch3/revise/frequency/norel/daf/

for comb in $(seq 9 1 64)
do
	#mkdir comb${comb}
	cd comb${comb}
	
	#shuf -n 24 ../../allAZ28.list > az.rand.comb${comb}
	#less ../../Norel_WTX_combinations.txt | cut -f${comb} > wtx.unrel.comb${comb}

	#vcftools --vcf /scratch/bell/mathur20/MQU/ch3/revise/vcfs/best66.old.auto.poly.nomiss.recode.vcf \
	#--freq \
	#--keep az.rand.comb${comb} \
	#--out az.rand.comb${comb}

	#vcftools --vcf /scratch/bell/mathur20/MQU/ch3/revise/vcfs/best66.old.auto.poly.nomiss.recode.vcf \
	#--freq \
	#--keep wtx.unrel.comb${comb} \
	#--out wtx.unrel.comb${comb}

	vcftools --vcf /scratch/bell/mathur20/MQU/ch3/revise/vcfs/best66.old.auto.poly.nomiss.recode.vcf \
	--freq \
	--positions /scratch/bell/mathur20/MQU/ch3/revise/annotation/snps/onlyAZ.synonymous_variant.SNPs \
	--keep az.rand.comb${comb} \
	--out az.rand.syn.comb${comb}

	vcftools --vcf /scratch/bell/mathur20/MQU/ch3/revise/vcfs/best66.old.auto.poly.nomiss.recode.vcf \
	--freq \
	--positions /scratch/bell/mathur20/MQU/ch3/revise/annotation/snps/onlyWTX.synonymous_variant.SNPs \
	--keep wtx.unrel.comb${comb} \
	--out wtx.unrel.syn.comb${comb}

	#less az.rand.comb${comb}.frq | cut -f6 | cut -f2 -d ":" > az.rand.comb${comb}.MAF
	#less wtx.unrel.comb${comb}.frq | cut -f6 | cut -f2 -d ":" > wtx.unrel.comb${comb}.MAF
	less az.rand.syn.comb${comb}.frq | cut -f6 | cut -f2 -d ":" > az.rand.syn.comb${comb}.MAF
	less  wtx.unrel.syn.comb${comb}.frq | cut -f6 | cut -f2 -d ":" >  wtx.unrel.syn.comb${comb}.MAF

	cd ../
done


