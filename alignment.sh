##########################################################################
###                             Samarth Mathur                          ###
###                        PhD Candidate, DeWoody Lab                   ###
###                           Purdue University                         ###
###                                                                     ###
###     Date Created: 05/31/20                  Last Modified: 08/22/20 ###
###########################################################################
###########################################################################
###                     alignment.sh                        			###
###########################################################################

cd $SLURM_SUBMIT_DIR
module purge
module load bioinfo
module load bwa
module load picard-tools
module load bedops
module load GATK/3.6.0
module load samtools

# STEP_0: Index fasta

bwa index /scratch/snyder/m/mathur20/MQU/ch3/align/chick/ref/chicken_genome.fa
samtools faidx /scratch/snyder/m/mathur20/MQU/ch3/align/chick/ref/chicken_genome.fa


# For each sample, create individual job file with each step #

while read -a line
do 
	echo "#!/bin/sh -l
#SBATCH -A gcore
#SBATCH -N 1
#SBATCH -n 5
#SBATCH -t 300:00:00
#SBATCH --job-name=${line[0]}_align
#SBATCH -e ${line[0]}_align
#SBATCH -o ${line[0]}_align

cd $SLURM_SUBMIT_DIR
module purge
module load bioinfo
module load bwa
module load picard-tools
module load bedops
module load GATK/3.6.0
module load samtools

cd /scratch/snyder/m/mathur20/MQU/ch3/align/chick/r4/step1/
bwa mem -t 50 -M -R \"@RG\tID:group1\tSM:${line[0]}\tPL:illumina\tLB:lib1\tPU:unit1\" \
/scratch/snyder/m/mathur20/MQU/ch3/align/chick/ref/chicken_genome.fa \
/scratch/snyder/m/mathur20/MQU/ch3/reads/r4/${line[0]}_R1.fastq \
/scratch/snyder/m/mathur20/MQU/ch3/reads/r4/${line[0]}_R2.fastq \
> ${line[0]}_chick.sam

PicardCommandLine ValidateSamFile I=${line[0]}_chick.sam MODE=SUMMARY O=${line[0]}_chick.sam.txt
PicardCommandLine SortSam INPUT=${line[0]}_chick.sam OUTPUT=sorted_${line[0]}_chick.bam SORT_ORDER=coordinate
PicardCommandLine MarkDuplicates INPUT=sorted_${line[0]}_chick.bam OUTPUT=dedup_${line[0]}_chick.bam METRICS_FILE=metrics_${line[0]}_chick.bam.txt
PicardCommandLine BuildBamIndex INPUT=dedup_${line[0]}_chick.bam


cd /scratch/snyder/m/mathur20/MQU/ch3/align/chick/r4/step2/

GenomeAnalysisTK -T RealignerTargetCreator \
-R /scratch/snyder/m/mathur20/MQU/ch3/align/chick/ref/chicken_genome.fa \
-I /scratch/snyder/m/mathur20/MQU/ch3/align/chick/r4/step1/dedup_${line[0]}_chick.bam \
-o forIndelRealigner.${line[0]}_chick.intervals

GenomeAnalysisTK -T IndelRealigner \
-R /scratch/snyder/m/mathur20/MQU/ch3/align/chick/ref/chicken_genome.fa \
-I /scratch/snyder/m/mathur20/MQU/ch3/align/chick/r4/step1/dedup_${line[0]}_chick.bam \
-targetIntervals forIndelRealigner.${line[0]}_chick.intervals \
-o realigned_${line[0]}_chick.bam \
&> log/indelrealign.${line[0]}.logfile.txt


cd /scratch/snyder/m/mathur20/MQU/ch3/align/chick/r4/step3/

PicardCommandLine FixMateInformation \
INPUT=/scratch/snyder/m/mathur20/MQU/ch3/align/chick/r4/step2/realigned_${line[0]}_chick.bam \
OUTPUT=${line[0]}.chick.sorted.dedup.realigned.fixmate.bam \
SO=coordinate \
CREATE_INDEX=true

cd /scratch/snyder/m/mathur20/MQU/ch3/align/chick/r4/step3/stats/

samtools depth -a /scratch/snyder/m/mathur20/MQU/ch3/align/chick/r4/step3/${line[0]}.chick.sorted.dedup.realigned.fixmate.bam \
| awk '{c++;s+=\$3}END{print s/c}' \
> ${line[0]}_chick_meandepth.txt

samtools depth -a /scratch/snyder/m/mathur20/MQU/ch3/align/chick/r4/step3/${line[0]}.chick.sorted.dedup.realigned.fixmate.bam \
| awk '{c++; if(\$3>0) total+=1}END{print (total/c)*100}' \
> ${line[0]}_chick_1x_breadth.txt
	
samtools flagstat /scratch/snyder/m/mathur20/MQU/ch3/align/chick/r4/step3/${line[0]}.chick.sorted.dedup.realigned.fixmate.bam \
> ${line[0]}_chick_mapped.txt

samtools stats /scratch/snyder/m/mathur20/MQU/ch3/align/chick/r4/step3/${line[0]}.chick.sorted.dedup.realigned.fixmate.bam \
> ${line[0]}_chick_stats.txt" > /scratch/snyder/m/mathur20/MQU/ch3/jobcodes/mapping/align/${line[0]}_chick_alignment.sh
done < /scratch/snyder/m/mathur20/MQU/ch3/reads/sample.list4




