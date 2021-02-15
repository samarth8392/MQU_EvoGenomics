##########################################################################
###                             Samarth Mathur                          ###
###                        PhD Candidate, DeWoody Lab                   ###
###                           Purdue University                         ###
###                                                                     ###
###     Date Created: 06/06/20                  Last Modified: 07/20/20 ###
###########################################################################
###########################################################################
###                     vartiant_calling.sh                   			###
###########################################################################

#Run Hapltoype Caller for each sample


while read -a line
do 
	echo "#!/bin/sh -l
#SBATCH -A gcore
#SBATCH -N 1
#SBATCH -n 10
#SBATCH -t 14-00:00:00
#SBATCH --mem=20G
#SBATCH --job-name=${line[0]}_variants
#SBATCH -e ${line[0]}_variants
#SBATCH -o ${line[0]}_variants

cd $SLURM_SUBMIT_DIR
module purge
module load bioinfo
module load bwa
module load picard-tools
module load bedops
module load GATK/4.1.6.0
module load samtools


cd /scratch/snyder/m/mathur20/MQU/ch3/variants/no_recal/
mkdir ${line[0]}
cd ${line[0]}
mkdir log

# STEP1: Run Haplotypecaller in GVCF mode #

gatk --java-options \"-Xmx20g -XX:+UseParallelGC -XX:ParallelGCThreads=10\" HaplotypeCaller  \
-R /scratch/snyder/m/mathur20/MQU/ch3/align/chick/ref/chicken_genome.fa \
-I /scratch/snyder/m/mathur20/MQU/ch3/align/chick/final_bam/dedup/${line[0]}.chick.sorted.dedup.realigned.fixmate.bam \
-ERC GVCF \
-O ${line[0]}.raw.norecal.variants.g.vcf &> log/${line[0]}.norecal.logfile.txt" > /scratch/snyder/m/mathur20/MQU/ch3/jobcodes/variants/calling/${line[0]}_chick_variantcall.sh

done < /scratch/snyder/m/mathur20/MQU/ch3/reads/sample.list


