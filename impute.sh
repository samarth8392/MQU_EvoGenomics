##########################################################################
###                             Samarth Mathur                          ###
###                        PhD Candidate, DeWoody Lab                   ###
###                           Purdue University                         ###
###                                                                     ###
###     Date Created: 08/04/20                 Last Modified: 09/04/20  ###
###########################################################################
###########################################################################
###                     impute.sh                                     	###
###########################################################################

cd $SLURM_SUBMIT_DIR
module load bioinfo
module load vcftools
module load samtools
module load beagle/3.1.2

# Get ref panel #
cd /scratch/snyder/m/mathur20/MQU/ch3/impute/beagle/

while read -a line
do
	wget http://gong_lab.hzau.edu.cn/static/imputeDB/download/species/chicken/panel/chr${line[1]}_chicken_impute.vcf.gz
done < /scratch/snyder/m/mathur20/MQU/ch3/variants/chr.list

gunzip *

# change chr name #

while read -a line
do
	awk '{if($0 !~ /^#/) $1=""; $2=""#; print}' chr${line[1]}_chicken_impute.vcf | awk -v VAR=${line[0]} '{if($0 !~ /^#/) print VAR, $0; else print $0}' > chr${line[1]}_chicken_impute.recode.vcf
done < /scratch/snyder/m/mathur20/MQU/ch3/variants/chr.list

cd ../
# Run beagle

while read -a line
do
	java -Xmx200g -jar /scratch/snyder/m/mathur20/MQU/ch3/softwares/beagle.18May20.d20.jar \
	nthreads=10 ne=50000 gp=true ap=true \
	gt=/scratch/snyder/m/mathur20/MQU/ch3/variants/angsd/best68.auto.longchr.angsd.recode.vcf \
	out=best68.auto.longchr.angsd
done < /scratch/snyder/m/mathur20/MQU/ch3/variants/chr.list

gunzip /scratch/snyder/m/mathur20/MQU/ch3/impute/beagle/best68.auto.longchr.angsd.vcf.gz

#convert vcf to impute file
vcftools --vcf /scratch/snyder/m/mathur20/MQU/ch3/impute/beagle/best68.auto.longchr.angsd.vcf \
--IMPUTE \
--out /scratch/snyder/m/mathur20/MQU/ch3/impute/beagle/best68.auto.longchr.angsd


# fastPHASE
module load fastPHASE/1.4.8
module load perl

#vcf to phase format

perl /scratch/snyder/m/mathur20/MQU/ch3/other_jobs/vcf-conversion-tools/vcf2fastPHASE.pl \
/scratch/snyder/m/mathur20/MQU/ch3/variants/angsd/best68.auto.longchr.angsd.recode.vcf \
best68.auto.longchr.angsd best68.auto.longchr.angsd 68


# IMPUTE FINAL
cd /scratch/snyder/m/mathur20/MQU/ch3/impute/beagle/

java -Xmx500g -jar /scratch/snyder/m/mathur20/MQU/ch3/softwares/beagle.18May20.d20.jar \
nthreads=20 ne=50000 gp=true ap=true burnin=100 iterations=200 \
gt=/scratch/snyder/m/mathur20/MQU/ch3/variants/angsd/best68.auto.longchr.angsd.final.recode.vcf \
out=best68.auto.longchr.beagle.phase



