#!/bin/sh -l
#SBATCH -A fnrquail
#SBATCH -N 1
#SBATCH -n 5
#SBATCH -t 14-00:00:00
#SBATCH --job-name=slimFst
#SBATCH -e slimFst
#SBATCH -o slimFst

##########################################################################
###                             Samarth Mathur                          ###
###                        PhD Candidate, DeWoody Lab                   ###
###                           Purdue University                         ###
###                                                                     ###
###     Date Created: 06/06/20                  Last Modified: 03/15/22 ###
###########################################################################
###########################################################################
###                     slimFst.sh                   		        	###
###########################################################################

cd $SLURM_SUBMIT_DIR
cd /scratch/bell/mathur20/MQU/ch3/revise/slim/1k/final_output/

module load bioinfo
module load vcftools
#module load plink/1.90b6.4

# Get Fst between pop1 and pop2 over time


#for gen in 125002 125252 125502 125752 126002 126302 126502
#do
#	for i in {41..100}
#	do 
#		vcftools --vcf /scratch/bell/mathur20/MQU/ch3/revise/slim/1k/final_output/vcf/1k.gen${gen}.run${i}.vcf \
#		--keep /scratch/bell/mathur20/MQU/ch3/revise/slim/1k/final_output/p1.list \
#		--freq \
#		--out freq/gen${gen}.run${i}.p1
#
#		vcftools --vcf /scratch/bell/mathur20/MQU/ch3/revise/slim/1k/final_output/vcf/1k.gen${gen}.run${i}.vcf \
#		--keep /scratch/bell/mathur20/MQU/ch3/revise/slim/1k/final_output/p2.list \
#		--freq \
#		--out freq/gen${gen}.run${i}.p2
#
#		vcftools --vcf /scratch/bell/mathur20/MQU/ch3/revise/slim/1k/final_output/vcf/1k.gen${gen}.run${i}.vcf \
#		--weir-fst-pop /scratch/bell/mathur20/MQU/ch3/revise/slim/1k/final_output/p1.list \
#		--weir-fst-pop /scratch/bell/mathur20/MQU/ch3/revise/slim/1k/final_output/p2.list\
#		--out fst/gen${gen}.run${i} &> fst/gen${gen}.run${i}.log
#	done
#done

#for gen in 125002 125252 125502 125752 126002 126302 126502
#do
#	for i in {1..40}
#	do 
		#grep -v "#" vcf/1k.gen${gen}.run${i}.vcf > geno/gen${gen}.run${i}.txt

		#grep -v "#" vcf/1k.gen${gen}.run${i}.vcf | cut -f8 | cut -d ";" -f2 | cut -f2 -d "=" \
		#> geno/gen${gen}.run${i}.selcoef.txt

		#less freq/gen${gen}.run${i}.p1.frq | cut -f6 | cut -f2 -d ":" > geno/gen${gen}.run${i}.p1.freq
		#less freq/gen${gen}.run${i}.p2.frq | cut -f6 | cut -f2 -d ":" > geno/gen${gen}.run${i}.p2.freq
#	done
#done

# Deleterious,neutral,and beneficial Fst

for gen in 125002 125252 125502 125752 126002 126302 126502
do
	for i in {1..100}
	do 
		for mut in neutral ben del
		do
			vcftools --vcf /scratch/bell/mathur20/MQU/ch3/revise/slim/1k/final_output/vcf/1k.gen${gen}.run${i}.vcf \
			--positions sites/gen${gen}.run${i}.${mut}.sites \
			--weir-fst-pop /scratch/bell/mathur20/MQU/ch3/revise/slim/1k/final_output/p1.list \
			--weir-fst-pop /scratch/bell/mathur20/MQU/ch3/revise/slim/1k/final_output/p2.list \
			--out fst/gen${gen}.run${i}.${mut} &> fst/gen${gen}.run${i}.${mut}.log
		done
	done
done

