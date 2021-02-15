##########################################################################
###                             Samarth Mathur                          ###
###                        PhD Candidate, DeWoody Lab                   ###
###                           Purdue University                         ###
###                                                                     ###
###     Date Created: 06/08/20                  Last Modified: 06/23/20 ###
###########################################################################
###########################################################################
###                     base_recalibration.sh                  			###
###########################################################################

# This script will use a pre-determined set of known variants and run base-recalibration 
# step twice to get the recalibrated bases

cd $SLURM_SUBMIT_DIR


cd /scratch/snyder/m/mathur20/MQU/ch3/jobcodes/mapping/base_recal/

while read -a line
do 
	echo "#!/bin/sh -l
#SBATCH -A fnrdewoody
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -t 300:00:00
#SBATCH --job-name=bsqr_chick_${line[0]}
#SBATCH -e bsqr_chick_${line[0]}
#SBATCH -o bsqr_chick_${line[0]}


cd $SLURM_SUBMIT_DIR
module purge
module load bioinfo
module load bwa
module load picard-tools
module load bedops
module load GATK/3.6.0
module load samtools

# Step1 : Call variants

cd /scratch/snyder/m/mathur20/MQU/ch3/align/chick/r2/base_rec/
mkdir ${line[0]}
cd ${line[0]}
#rm *
mkdir plots
mkdir log

#cd /scratch/snyder/m/mathur20/MQU/ch3/align/chick/base_rec/no_recal/
#samtools index ${line[0]}.chick.norecal.bam

cd /scratch/snyder/m/mathur20/MQU/ch3/align/chick/base_rec/${line[0]}/

GenomeAnalysisTK -T BaseRecalibrator \
-R /scratch/snyder/m/mathur20/MQU/ch3/align/chick/ref/chicken_genome.fa \
-knownSites /scratch/snyder/m/mathur20/MQU/ch3/deep_cov/variants/hardfiltered_deep.chick.recode.vcf \
-I /scratch/snyder/m/mathur20/MQU/ch3/align/chick/r2/step3/${line[0]}.chick.sorted.dedup.realigned.fixmate.bam \
-o ${line[0]}.recal.1.table \
-cov ReadGroupCovariate \
-cov QualityScoreCovariate \
-cov CycleCovariate &> log/${line[0]}.recal.1.log.txt

GenomeAnalysisTK -T BaseRecalibrator \
-R /scratch/snyder/m/mathur20/MQU/ch3/align/chick/ref/chicken_genome.fa \
-knownSites /scratch/snyder/m/mathur20/MQU/ch3/deep_cov/variants/hardfiltered_deep.chick.recode.vcf \
-I /scratch/snyder/m/mathur20/MQU/ch3/align/chick/r2/step3/${line[0]}.chick.sorted.dedup.realigned.fixmate.bam \
-BQSR ${line[0]}.recal.1.table \
-o ${line[0]}.recal.2.table &> log/${line[0]}.baserecal.log.txt

GenomeAnalysisTK -T AnalyzeCovariates \
-R /scratch/snyder/m/mathur20/MQU/ch3/align/chick/ref/chicken_genome.fa \
-before ${line[0]}.recal.1.table \
-after ${line[0]}.recal.2.table \
-plots plots/${line[0]}_12_bsqr_plot.pdf \
-csv plots/${line[0]}_12_bsqr_plot.csv

GenomeAnalysisTK -T PrintReads \
-R /scratch/snyder/m/mathur20/MQU/ch3/align/chick/ref/chicken_genome.fa \
-BQSR ${line[0]}.recal.1.table \
-I /scratch/snyder/m/mathur20/MQU/ch3/align/chick/r2/step3/${line[0]}.chick.sorted.dedup.realigned.fixmate.bam \
-o ${line[0]}.chick.recal.1.bam &> log/${line[0]}.1.printreads.log.txt

# ROUND 2 # 

GenomeAnalysisTK -T BaseRecalibrator \
-R /scratch/snyder/m/mathur20/MQU/ch3/align/chick/ref/chicken_genome.fa \
-knownSites /scratch/snyder/m/mathur20/MQU/ch3/deep_cov/variants/hardfiltered_deep.chick.recode.vcf \
-I ${line[0]}.chick.recal.1.bam \
-o ${line[0]}.recal.3.table \
-cov ReadGroupCovariate \
-cov QualityScoreCovariate \
-cov CycleCovariate &> log/${line[0]}.recal.3.log.txt

GenomeAnalysisTK -T BaseRecalibrator \
-R /scratch/snyder/m/mathur20/MQU/ch3/align/chick/ref/chicken_genome.fa \
-knownSites /scratch/snyder/m/mathur20/MQU/ch3/deep_cov/variants/hardfiltered_deep.chick.recode.vcf \
-I ${line[0]}.chick.recal.1.bam \
-BQSR ${line[0]}.recal.3.table \
-o ${line[0]}.recal.4.table &> log/${line[0]}.baserecal2.log.txt

GenomeAnalysisTK -T AnalyzeCovariates \
-R /scratch/snyder/m/mathur20/MQU/ch3/align/chick/ref/chicken_genome.fa \
-before ${line[0]}.recal.3.table \
-after ${line[0]}.recal.4.table \
-plots plots/${line[0]}_34_bsqr_plot.pdf \
-csv plots/${line[0]}_34_bsqr_plot.csv

GenomeAnalysisTK -T PrintReads \
-R /scratch/snyder/m/mathur20/MQU/ch3/align/chick/ref/chicken_genome.fa \
-BQSR ${line[0]}.recal.3.table \
-I ${line[0]}.chick.recal.1.bam \
-o ${line[0]}.chick.recal.2.bam &> log/${line[0]}.2.printreads.log.tx" > /scratch/snyder/m/mathur20/MQU/ch3/jobcodes/mapping/base_recal/${line[0]}_chick_baserecal.sh

done < /scratch/snyder/m/mathur20/MQU/ch3/reads/sample.list2
