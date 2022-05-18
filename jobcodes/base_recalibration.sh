#!/bin/sh -l
#SBATCH -A fnrdewoody
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -t 300:00:00
#SBATCH --job-name=bsqr_chick
#SBATCH -e bsqr_chick
#SBATCH -o bsqr_chick

##########################################################################
###                             Samarth Mathur                          ###
###                        PhD Candidate, DeWoody Lab                   ###
###                           Purdue University                         ###
###                                                                     ###
###     Date Created: 06/08/20                  Last Modified: 01/24/22 ###
###########################################################################
###########################################################################
###                     base_recalibration.sh                  			###
###########################################################################

cd $SLURM_SUBMIT_DIR


#module purge
#module load bioinfo
#module load bwa
#module load picard-tools
#module load bedops
#module load GATK/3.6.0
#module load samtools

# Make before after plots for MQU Base recal #

#while read -a line
#do
#	GenomeAnalysisTK -T BaseRecalibrator \
#	-R /scratch/snyder/m/mathur20/MQU/ch3/align/MQU/ref/MQU_male.min500.fa \
#	-knownSites /scratch/snyder/m/mathur20/MQU/ch3/align/MQU/ref/raw_indels.vcf \
#	-knownSites /scratch/snyder/m/mathur20/MQU/ch3/align/MQU/ref/raw_snps.vcf \
#	-I /scratch/snyder/m/mathur20/MQU/ch3/align/MQU/step3/${line[0]}.mqu.sorted.dedup.realigned.fixmate.bam \
#	-BQSR /scratch/snyder/m/mathur20/MQU/ch3/align/MQU/step4/${line[0]}.mqu.sorted.dedup.realigned.fixmate.recal_data.table \
#	-o /scratch/snyder/m/mathur20/MQU/ch3/align/MQU/base_recal/${line[0]}.mqu.post.recal.table 
#done < /scratch/snyder/m/mathur20/MQU/ch3/reads/sample.list


#while read -a line
#do
#	GenomeAnalysisTK -T AnalyzeCovariates \
#	-R /scratch/snyder/m/mathur20/MQU/ch3/align/MQU/ref/MQU_male.min500.fa \
#	-before /scratch/snyder/m/mathur20/MQU/ch3/align/MQU/step4/${line[0]}.mqu.sorted.dedup.realigned.fixmate.recal_data.table \
#	-after /scratch/snyder/m/mathur20/MQU/ch3/align/MQU/base_recal/${line[0]}.mqu.post.recal.table \
#	-plots /scratch/snyder/m/mathur20/MQU/ch3/align/MQU/base_recal/plots/${line[0]}_bsqr_plot.pdf \
#	-csv /scratch/snyder/m/mathur20/MQU/ch3/align/MQU/base_recal/plots/${line[0]}_bsqr_plot.csv
#done < /scratch/snyder/m/mathur20/MQU/ch3/reads/sample.list

while read -a line
#for line in E9569 E9570
do 
	echo "#!/bin/sh -l
#SBATCH -A fnrquail
#SBATCH -N 1
#SBATCH -n 8
#SBATCH --mem=14G
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
module load GATK/4.1.6.0
module load samtools

# Step1 : Call variants

cd /scratch/bell/mathur20/MQU/ch3/revise/base_recal/
mkdir ${line[0]}
cd ${line[0]}
mkdir plots
mkdir log

cd /scratch/bell/mathur20/MQU/ch3/revise/base_recal/${line}/

gatk --java-options \"-Xmx14g\" BaseRecalibrator \
-R /scratch/bell/mathur20/MQU/ch3/align/chick/ref/chicken_genome.fa \
--known-sites /scratch/bell/mathur20/MQU/ch3/deep_cov/variants/deep_coverage_auto.recode.vcf \
-I /scratch/bell/mathur20/MQU/ch3/align/chick/final_bam/dedup/${line[0]}.chick.sorted.dedup.realigned.fixmate.bam \
-O ${line}.recal.1.table &> log/${line}.recal.1.log.txt

gatk --java-options \"-Xmx14g\" ApplyBQSR \
-R /scratch/bell/mathur20/MQU/ch3/align/chick/ref/chicken_genome.fa \
-I /scratch/bell/mathur20/MQU/ch3/align/chick/final_bam/dedup/${line[0]}.chick.sorted.dedup.realigned.fixmate.bam \
--bqsr-recal-file ${line}.recal.1.table \
-O ${line}.chick.recal.1.bam &> log/${line}.1.ApplyBQSR.log.txt


# ROUND 2 # 

gatk --java-options \"-Xmx14g\" BaseRecalibrator \
-R /scratch/bell/mathur20/MQU/ch3/align/chick/ref/chicken_genome.fa \
--known-sites /scratch/bell/mathur20/MQU/ch3/deep_cov/variants/deep_coverage_auto.recode.vcf \
-I ${line}.chick.recal.1.bam \
-O ${line}.recal.2.table &> log/${line}.recal.2.log.txt

gatk --java-options \"-Xmx14g\" ApplyBQSR \
-R /scratch/bell/mathur20/MQU/ch3/align/chick/ref/chicken_genome.fa \
-I ${line}.chick.recal.1.bam \
--bqsr-recal-file ${line}.recal.2.table \
-O ${line}.chick.recal.2.bam &> log/${line}.2.ApplyBQSR.log.txt" \
> /scratch/bell/mathur20/MQU/ch3/revise/jobcodes/per_ind/${line}/${line}.chick_baserecall.sh
done < /scratch/bell/mathur20/MQU/ch3/revise/lists/best66.samples.list

# Submit jobs

while read -a line
do
	cd /scratch/bell/mathur20/MQU/ch3/revise/errors/per_ind/${line}
	sbatch /scratch/bell/mathur20/MQU/ch3/revise/jobcodes/per_ind/${line}/${line}.chick_baserecall.sh
done < /scratch/bell/mathur20/MQU/ch3/revise/lists/best66.samples.list

