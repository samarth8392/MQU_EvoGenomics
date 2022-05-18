#!/bin/sh -l
#SBATCH -A fnrdewoody
#SBATCH -N 1
#SBATCH -n 12
#SBATCH -t 14-00:00:00
#SBATCH --job-name=shapeit_step3
#SBATCH -e shapeit_step3
#SBATCH -o shapeit_step3


cd $SLURM_SUBMIT_DIR
module load bioinfo
module load vcftools
module load samtools

#cd /scratch/bell/mathur20/MQU/ch3/align/chick/final_bam/recal/
#while read -a name
#do
#	samtools index ${name}.chick.recal.final.bam  ${name}.chick.recal.final.bai
#done < /scratch/bell/mathur20/MQU/ch3/revise/lists/best66.samples.list

#Breaking VCF by chromosome
cd /scratch/bell/mathur20/MQU/ch3/revise/phase/r2/vcfBychr/

#while read -a chr
#do
#	vcftools --vcf /scratch/bell/mathur20/MQU/ch3/revise/vcfs/best66.old.auto.noSing.nomiss.recode.vcf \
#	--chr ${chr} \
#	--recode --recode-INFO-all \
#	--out best66.old.auto.noSing.nomiss.${chr}
#done < /scratch/bell/mathur20/MQU/ch3/revise/lists/chick_chr.list


#extracting information from sequence reads

#cd /scratch/bell/mathur20/MQU/ch3/revise/phase/r2/pirs/
#while read -a chr
#do
#	~/shapeit.v2.904/extractPIRs.v1.r68/extractPIRs \
#	--bam /scratch/bell/mathur20/MQU/ch3/revise/phase/r1/bamlist/best66.${chr}.bamlist \
#	--vcf /scratch/bell/mathur20/MQU/ch3/revise/phase/r2/vcfBychr/best66.old.auto.noSing.nomiss.${chr}.recode.vcf \
#	--base-quality 20 --read-quality 20 \
#	--out best66.${chr}.PIRsList
#done < /scratch/bell/mathur20/MQU/ch3/revise/lists/chick_chr.list

#running shape it using read information
cd /scratch/bell/mathur20/MQU/ch3/revise/phase/r2/assemble/
while read -a chr
do
	~/shapeit.v2.904/bin/shapeit -assemble --thread 30 \
	--input-vcf /scratch/bell/mathur20/MQU/ch3/revise/phase/r2/vcfBychr/best66.old.auto.noSing.nomiss.${chr}.recode.vcf \
	--input-pir /scratch/bell/mathur20/MQU/ch3/revise/phase/r2/pirs/best66.${chr}.PIRsList \
	--states 1000 --burn 200 --prune 210 --main 2000 --force \
	-O /scratch/bell/mathur20/MQU/ch3/revise/phase/r2/assemble/best66.phase.${chr}
done < /scratch/bell/mathur20/MQU/ch3/revise/lists/chick_chr.list
