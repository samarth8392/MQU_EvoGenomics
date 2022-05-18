#!/bin/sh -l
#SBATCH -A fnrdewoody
#SBATCH -N 1
#SBATCH -n 15
#SBATCH -t 14-00:00:00
#SBATCH --mem=90G
#SBATCH --job-name=beagle_error
#SBATCH -e beagle_error
#SBATCH -o beagle_error


cd $SLURM_SUBMIT_DIR
module load bioinfo
module load vcftools

cd /scratch/bell/mathur20/MQU/ch3/revise/phase/beagle/vcfBychr/

#Breaking VCF by chromosome (chromosome names =1,2,3..)

#for chr in $(seq 1 1 15; seq 17 1 28)
#do
#	vcftools --vcf /scratch/bell/mathur20/MQU/ch3/revise/vcfs/best66.old.auto.noSing.nomiss.rename \
#	--chr ${chr} \
#	--recode --recode-INFO-all \
#	--out best66.old.auto.noSing.nomiss.chr${chr}
#done

# Running phase by chr usping chicken recombination data (10 times)


cd /scratch/bell/mathur20/MQU/ch3/revise/phase/beagle/phase/

for i in $(seq 1 1 10)
do
	mkdir run${i}
	for chr in $(seq 1 1 15; seq 17 1 28)
	do
		java -Xmx90g -jar ~/beagle.28Jun21.220.jar \
		burnin=200 iterations=500 phase-states=500 \
		map=/scratch/bell/mathur20/MQU/ch3/revise/phase/beagle/chickMap/chick.chr${chr}.map \
		nthreads=52 ne=50000 \
		gt=/scratch/bell/mathur20/MQU/ch3/revise/phase/beagle/vcfBychr/best66.old.auto.noSing.nomiss.chr${chr}.recode.vcf \
		out=run${i}/best66.old.auto.noSing.nomiss.chr${chr}
	done
done

# Compare switch errors between runs

#count1=0
#count2=0
#for i in $(seq 1 1 10)
#do
#	((count1=count1+1))
#	count2=0
#	for j in $(seq 1 1 10)
#	do
#		((count2=count2+1))
#		if [[ ${count2} -gt ${count1} ]]
#		then
#			if [[ ${i} != ${j} ]]
#			then
#				for chr in $(seq 1 1 15; seq 17 1 28)
#				do	
#					vcftools --diff-switch-error \
#					--gzvcf /scratch/bell/mathur20/MQU/ch3/revise/phase/beagle/phase/run${i}/best66.old.auto.noSing.nomiss.chr${chr}.vcf.gz \
#					--gzdiff /scratch/bell/mathur20/MQU/ch3/revise/phase/beagle/phase/run${j}/best66.old.auto.noSing.nomiss.chr${chr}.vcf.gz \
#					--out /scratch/bell/mathur20/MQU/ch3/revise/phase/beagle/switch/chr${chr}.run${i}-run${j}
#				done
#			fi
#		fi
#	done
#done









