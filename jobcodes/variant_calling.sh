#!/bin/sh -l
#SBATCH -A fnrdewoody
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -t 300:00:00
#SBATCH --job-name=variant_call
#SBATCH -e variant_call
#SBATCH -o variant_call

##########################################################################
###                             Samarth Mathur                          ###
###                        PhD Candidate, DeWoody Lab                   ###
###                           Purdue University                         ###
###                                                                     ###
###     Date Created: 06/06/20                  Last Modified: 11/04/21 ###
###########################################################################
###########################################################################
###                     vartiant_calling.sh                   			###
###########################################################################

cd $SLURM_SUBMIT_DIR
#module purge
#module load bioinfo
#module load bwa
#module load picard-tools
#module load bedops
#module load GATK/4.1.6.0
#module load samtools


#while read -a line
#for line in E9569 E9570
while read -a line
do 
	echo "#!/bin/sh -l
#SBATCH -A fnrquail
#SBATCH -N 1
#SBATCH -n 8
#SBATCH -t 14-00:00:00
#SBATCH --mem=14G
#SBATCH --job-name=${line[0]}_HapCall
#SBATCH -e ${line[0]}_HapCall
#SBATCH -o ${line[0]}_HapCall

cd $SLURM_SUBMIT_DIR
module purge
module load bioinfo
module load bwa
module load picard-tools
module load bedops
module load GATK/4.1.6.0
module load samtools


cd /scratch/bell/mathur20/MQU/ch3/revise/gatk/gvcf/

# STEP1: Run Haplotypecaller in GVCF mode #

gatk --java-options \"-Xmx14g -XX:+UseParallelGC -XX:ParallelGCThreads=8\" HaplotypeCaller  \
-R /scratch/bell/mathur20/MQU/ch3/align/chick/ref/chicken_genome.fa \
-I /scratch/bell/mathur20/MQU/ch3/align/chick/final_bam/recal/${line}.chick.recal.final.bam \
-ERC GVCF \
-O ${line}.raw.variants.g.vcf" > /scratch/bell/mathur20/MQU/ch3/revise/jobcodes/per_ind/${line}/${line}.chick_variantcall.sh
done < /scratch/bell/mathur20/MQU/ch3/revise/lists/first64.samples.list
#done

while read -a line
do
	cd /scratch/bell/mathur20/MQU/ch3/revise/errors/per_ind/${line}/
	sbatch /scratch/bell/mathur20/MQU/ch3/revise/jobcodes/per_ind/${line}/${line}.chick_variantcall.sh
done < /scratch/bell/mathur20/MQU/ch3/revise/lists/first64.samples.list

