#!/bin/sh -l
#SBATCH -A fnrdewoody
#SBATCH -N 1
#SBATCH -n 5
#SBATCH -t 14-00:00:00
#SBATCH --job-name=geno_age
#SBATCH -e geno_age
#SBATCH -o geno_age

##########################################################################
###                             Samarth Mathur                          ###
###                        PhD Candidate, DeWoody Lab                   ###
###                           Purdue University                         ###
###                                                                     ###
###     Date Created: 06/06/20                  Last Modified: 01/18/22 ###
###########################################################################
###########################################################################
###                     geno_age.sh                   		        	###
###########################################################################

# Get genotypes for mutations with ages

cd $SLURM_SUBMIT_DIR
module purge
module load bioinfo
module load vcftools

cd /scratch/bell/mathur20/MQU/ch3/revise/slim/1k/repeats

for i in {1..100}
do
	grep -v "initialize" run${i}/1k.postbot.run${i}.R2 \
	| grep -v "Working"  | sed '1,8d' > ../final_output/csv/run${i}.postbot.csv

done

for i in {1..100}
do
	

	

