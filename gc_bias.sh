#!/bin/sh -l
#SBATCH -A fnrdewoody
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -t 300:00:00
#SBATCH --job-name=gc_bias_2
#SBATCH -e gc_bias
#SBATCH -o gc_bias

##########################################################################
###                             Samarth Mathur                          ###
###                        PhD Candidate, DeWoody Lab                   ###
###                           Purdue University                         ###
###                                                                     ###
###     Date Created: 07/22/20                  Last Modified: 07/23/20 ###
###########################################################################
###########################################################################
###                     ggc_bias.sh                        				###
###########################################################################

cd $SLURM_SUBMIT_DIR
module purge
module load bioinfo
module load bwa
module load picard-tools
module load bedops
module load GATK/3.6.0
module load samtools

while read -a line
do
	PicardCommandLine CollectGcBiasMetrics \
     I=/scratch/snyder/m/mathur20/MQU/ch3/align/chick/r2/step1/sorted_${line[0]}_chick.bam \
     O=/scratch/snyder/m/mathur20/MQU/ch3/align/chick/GC_bias/${line[0]}_gc_bias_metrics.txt \
     CHART=/scratch/snyder/m/mathur20/MQU/ch3/align/chick/GC_bias/${line[0]}_gc_bias_metrics.pdf \
     S=/scratch/snyder/m/mathur20/MQU/ch3/align/chick/GC_bias/${line[0]}_summary_metrics.txt \
     R=/scratch/snyder/m/mathur20/MQU/ch3/align/chick/ref/chicken_genome.fa
done < /scratch/snyder/m/mathur20/MQU/ch3/reads/sample.list2
