##########################################################################
###                             Samarth Mathur                          ###
###                        PhD Candidate, DeWoody Lab                   ###
###                           Purdue University                         ###
###                                                                     ###
###     Date Created: 08/17/20                  Last Modified: 08/17/20 ###
###########################################################################
###########################################################################
###                     finestructure.sh                        		###
###########################################################################



#see : https://people.maths.bris.ac.uk/~madjl/finestructure/finestructure_info.html
# see : https://people.maths.bris.ac.uk/~madjl/finestructure/manualse13.html


cd $SLURM_SUBMIT_DIR
#module load bioinfo
#module load parallel/20180222
#module load vcftools
#module load perl
#module load plink


#Using chicken recMap
cd /scratch/snyder/m/mathur20/MQU/ch3/structure/fineS/run2/

# Step1: Phase and impute by chr (using BEAGLE)

cd input/beagle/ 

for i in `seq 1 1 28`
do

	java -Xmx200g -jar /scratch/snyder/m/mathur20/MQU/ch3/softwares/beagle.18May20.d20.jar \
	nthreads=20 ne=100000 gp=true ap=true burnin=50 iterations=100 \
	chrom=$i \
	map=/scratch/snyder/m/mathur20/MQU/ch3/impute/chicken/chick.chr$i.map \
	gt=/scratch/snyder/m/mathur20/MQU/ch3/variants/angsd/best68.auto.longch.noSing.rename.vcf \
	out=best66.noSing.chr$i.beagle.phase
done

gunzip *phase.vcf.gz

# Step2: Convert vcf to impute (.haps) and then to right input format

cd input/
for i in `seq 1 1 28` 
do
	plink --chr-set 33 no-xy no-mt --recode12 \
	--vcf beagle/best66.noSing.chr$i.beagle.phase.vcf \
	--out hap/best66.noSing.chr$i.beagle.phase
done

for i in `seq 1 1 28` 
do
	perl /scratch/snyder/m/mathur20/MQU/ch3/softwares/fs_4.1.1/plink2chromopainter.pl \
	-p=hap/best66.noSing.chr$i.beagle.phase.ped \
	-m=hap/best66.noSing.chr$i.beagle.phase.map \
	-o=fineS/best66.noSing.chr$i.beagle.phase
done

# Step3: Get chicken recombination map by chr

for i in `seq 1 1 28` 
do
	cat /scratch/snyder/m/mathur20/MQU/ch3/impute/chicken/chick.chr$i.map \
	| awk '{ print $4, $3 }' | awk '{sub(/^M/,"")}1' > recMap/chicken.chr$i.recMap.txt
	
	perl /scratch/snyder/m/mathur20/MQU/ch3/softwares/fs_4.1.1/convertrecfile.pl \
	-D cdf -U c \
	fineS/best66.noSing.chr$i.beagle.phase recMap/chicken.chr$i.recMap.txt recMap/chicken.chr$i.recombfile
done


# Step4a: Run Chromopainter

cd /scratch/snyder/m/mathur20/MQU/ch3/structure/fineS/run2/painter/
 
for i in `seq 1 1 28`
for i in 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 17 18 19 20 21 22 23 24 25 26 27 28 
do
	/scratch/snyder/m/mathur20/MQU/ch3/softwares/fs_4.1.1/fs \
	best66.noSing.chr$i.cp \
	-n -hpc 1 \
	-phasefiles /scratch/snyder/m/mathur20/MQU/ch3/structure/fineS/run2/input/fineS/best66.noSing.chr$i.beagle.phase \
	-recombfiles /scratch/snyder/m/mathur20/MQU/ch3/structure/fineS/run2/input/recMap/chicken.chr$i.recombfile \
	-idfile /scratch/snyder/m/mathur20/MQU/ch3/structure/best66.list -go 
done

for i in `seq 1 1 28`
for i in 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 17 18 19 20 21 22 23 24 25 26 27 28
do
	cd /scratch/snyder/m/mathur20/MQU/ch3/structure/fine_structure/final/painter/best68.auto.long.chr$i/commandfiles/
	cat commandfile1.txt | awk '{ gsub("-in -iM --emfilesonly", "-in -iM --emfilesonly -n 100000 -k 20") ; print $0 }' > commandfile1a.txt
done

# Create jobfiles
for i in `seq 1 1 28`
#for i in 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 17 18 19 20 21 22 23 24 25 26 27 28 
do
	count=1
	cat best66.noSing.chr$i/commandfiles/commandfile1.txt | while read -r line
	do
		echo "#!/bin/sh -l
#SBATCH -A fnrfish
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -t 14-00:00:00
#SBATCH --job-name=fineS_chr$i.c1.j$count
#SBATCH -e fineS_chr$i.c1.j$count
#SBATCH -o fineS_chr$i.c1.j$count

cd /scratch/snyder/m/mathur20/MQU/ch3/structure/fineS/run2/painter/
/scratch/snyder/m/mathur20/MQU/ch3/softwares/fs_4.1.1/$line" > /scratch/snyder/m/mathur20/MQU/ch3/structure/fineS/run2/jobs/cmd1/chr$i.cmd1.j$count.job
count=$[$count+1]
	done 
done

# Submit jobs

cd /scratch/snyder/m/mathur20/MQU/ch3/structure/fineS/run2/jobs/cmd1/

for prefix in $(ls| uniq)
do
	cd /scratch/snyder/m/mathur20/MQU/ch3/structure/fineS/run2/errors/
	sbatch ../jobs/cmd1/$prefix
done

for i in `seq 1 1 28`
do
	/scratch/snyder/m/mathur20/MQU/ch3/softwares/fs_4.1.1/fs best66.noSing.chr$i.cp -go
done

#Ste4b: ROUND2
#Create jobfiles 

cd /scratch/snyder/m/mathur20/MQU/ch3/structure/fine_structure/final/painter/
for i in `seq 1 1 28`
do
	count=1
	cat best66.noSing.chr$i/commandfiles/commandfile2.txt | while read -r line
	do
		echo "#!/bin/sh -l
#SBATCH -A fnrfish
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -t 14-00:00:00
#SBATCH --job-name=finS_ch$i.c2.j$count
#SBATCH -e finS_ch$i.c2.j$count
#SBATCH -o finS_ch$i.c2.j$count

cd /scratch/snyder/m/mathur20/MQU/ch3/structure/fineS/run2/painter/
/scratch/snyder/m/mathur20/MQU/ch3/softwares/fs_4.1.1/$line" > /scratch/snyder/m/mathur20/MQU/ch3/structure/fineS/run2/jobs/cmd2/chr$i.cmd2.j$count.job
count=$[$count+1]
	done 
done

# Submit jobs

cd /scratch/snyder/m/mathur20/MQU/ch3/structure/fineS/run2/jobs/cmd2/

for prefix in $(ls| uniq)
do
	cd /scratch/snyder/m/mathur20/MQU/ch3/structure/fineS/run2/errors/cmd2
	sbatch ../../jobs/cmd2/$prefix
done

for i in `seq 1 1 28`
do
	cd /scratch/snyder/m/mathur20/MQU/ch3/structure/fineS/run2/painter/
	/scratch/snyder/m/mathur20/MQU/ch3/softwares/fs_4.1.1/fs best66.noSing.chr$i.cp -go
done

#Step4c: ROUND3
#Create jobfiles

cd /scratch/snyder/m/mathur20/MQU/ch3/structure/fineS/run2/painter/
for i in `seq 1 1 28`
do
	line=$(head -n 1 best66.noSing.chr$i/commandfiles/commandfile3.txt) 
	/scratch/snyder/m/mathur20/MQU/ch3/softwares/fs_4.1.1/$line
done

for i in `seq 1 1 28`
do
	cd /scratch/snyder/m/mathur20/MQU/ch3/structure/fineS/run2/painter/
	/scratch/snyder/m/mathur20/MQU/ch3/softwares/fs_4.1.1/fs best66.noSing.chr$i.cp -go
done

# Step5a" LAST ROUND

cd /scratch/snyder/m/mathur20/MQU/ch3/structure/fineS/run2/painter/
for i in `seq 1 1 28`
do
	count=1
	cat best66.noSing.chr$i/commandfiles/commandfile4.txt | while read -r line
	do
		echo "#!/bin/sh -l
#SBATCH -A fnrfish
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -t 14-00:00:00
#SBATCH --job-name=fineS_ch$i.c4.j$count
#SBATCH -e fineS_ch$i.c4.j$count
#SBATCH -o fineS_ch$i.c4.j$count

cd /scratch/snyder/m/mathur20/MQU/ch3/structure/fineS/run2/painter/
/scratch/snyder/m/mathur20/MQU/ch3/softwares/fs_4.1.1/$line" > /scratch/snyder/m/mathur20/MQU/ch3/structure/fineS/run2/jobs/cmd4/chr$i.cmd4.j$count.job
count=$[$count+1]
	done 
done


#Submit jobs

cd /scratch/snyder/m/mathur20/MQU/ch3/structure/fineS/run2/jobs/cmd4/

for prefix in $(ls| uniq)
do
	cd /scratch/snyder/m/mathur20/MQU/ch3/structure/fineS/run2/errors/cmd4/
	sbatch ../../jobs/cmd4/$prefix
done

for i in `seq 1 1 28`
do
	cd /scratch/snyder/m/mathur20/MQU/ch3/structure/fineS/run2/painter/
	/scratch/snyder/m/mathur20/MQU/ch3/softwares/fs_4.1.1/fs best66.noSing.chr$i.cp -go
done
