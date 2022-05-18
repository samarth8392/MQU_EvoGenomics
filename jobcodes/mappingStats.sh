#!/bin/sh -l
#SBATCH -A fnrdewoody
#SBATCH -N 1
#SBATCH -n 32
#SBATCH -t 14-00:00:00
#SBATCH --job-name=mapstats
#SBATCH -e mapstats
#SBATCH -o mapstats

##########################################################################
###                             Samarth Mathur                          ###
###                        PhD Candidate, DeWoody Lab                   ###
###                           Purdue University                         ###
###                                                                     ###
###     Date Created: 08/04/20                 Last Modified: 12/31/21 	###
###########################################################################
###########################################################################
###                     mapstats.sh                   					###
###########################################################################

cd $SLURM_SUBMIT_DIR
module purge
module load bioinfo
module load samtools


cd /scratch/bell/mathur20/MQU/ch3/revise/mapstats/

while read -a line
do
	samtools depth -a /scratch/bell/mathur20/MQU/ch3/align/chick/final_bam/recal/${line[0]}.chick.recal.final.bam \
	| awk '{c++;s+=$3}END{print s/c}' \
	> ${line[0]}_chick_meandepth.txt

	samtools depth /scratch/bell/mathur20/MQU/ch3/align/chick/final_bam/recal/${line[0]}.chick.recal.final.bam \
	| awk '{c++;s+=$3}END{print s/c}' \
	> ${line[0]}_chick_meandepth2.txt

	samtools depth -a /scratch/bell/mathur20/MQU/ch3/align/chick/final_bam/recal/${line[0]}.chick.recal.final.bam \
	| awk '{c++; if($3>0) total+=1}END{print (total/c)*100}' \
	> ${line[0]}_chick_1x_breadth.txt

	samtools depth -a /scratch/bell/mathur20/MQU/ch3/align/chick/final_bam/recal/${line[0]}.chick.recal.final.bam \
	| awk '{c++; if($3>5) total+=1}END{print (total/c)*100}' \
	> ${line[0]}_chick_5x_breadth.txt
done < /scratch/bell/mathur20/MQU/ch3/revise/lists/best66.samples.list
