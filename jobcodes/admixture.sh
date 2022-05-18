#!/bin/sh -l
#SBATCH -A fnrquail
#SBATCH -N 1
#SBATCH -n 30
#SBATCH -t 14-00:00:00
#SBATCH --job-name=admix
#SBATCH -e admix
#SBATCH -o admix

##########################################################################
###                             Samarth Mathur                          ###
###                        PhD Candidate, DeWoody Lab                   ###
###                           Purdue University                         ###
###                                                                     ###
###     Date Created: 06/06/20                  Last Modified: 01/25/22 ###
###########################################################################
###########################################################################
###                     admixture.sh                   		        	###
###########################################################################

cd $SLURM_SUBMIT_DIR
module purge
module load bioinfo
module load plink/1.90b6.4


#Convert rename vcf (chr1, chr2,...) to PLINK

cd /scratch/bell/mathur20/MQU/ch3/revise/admix

plink --vcf /scratch/bell/mathur20/MQU/ch3/revise/vcfs/best66.old.auto.noSing.nomiss.rename \
--make-bed --allow-extra-chr --chr-set 33 no-xy no-mt \
--out best66.old.auto.noSing.nomiss.rename

#Run admixture

for k in $(seq 1 1 10)
do
	/scratch/bell/mathur20/MQU/ch3/softwares/admixture_linux-1.3.0/admixture \
	--cv best66.old.auto.noSing.nomiss.rename.bed $k -j30
done