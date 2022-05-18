#!/bin/sh -l
#SBATCH -A fnrdewoody
#SBATCH -N 1
#SBATCH -n 5
#SBATCH -t 14-00:00:00
#SBATCH --job-name=geno_age
#SBATCH -e geno_age
#SBATCH -o geno_age

##########################################################################
###                             Samarth Mathur                          ###
###                        PhD Candidate, DeWoody Lab                   ###
###                           Purdue University                         ###
###                                                                     ###
###     Date Created: 06/06/20                  Last Modified: 01/18/22 ###
###########################################################################
###########################################################################
###                     geno_age.sh                   		        	###
###########################################################################

# Get genotypes for mutations with ages

cd $SLURM_SUBMIT_DIR
module purge
module load bioinfo
module load vcftools


cd /scratch/bell/mathur20/MQU/ch3/revise/age/geno/

for pop in az wtx
do
	vcftools --vcf /scratch/bell/mathur20/MQU/ch3/revise/vcfs/best66.old.auto.noSing.nomiss.rename \
    --positions /scratch/bell/mathur20/MQU/ch3/revise/age/sites/${pop}.del.sites \
    --keep /scratch/bell/mathur20/MQU/ch3/revise/lists/${pop}.best66.list \
    --extract-FORMAT-info GT \
    --out ${pop}.del

    vcftools --vcf /scratch/bell/mathur20/MQU/ch3/revise/vcfs/best66.old.auto.noSing.nomiss.rename \
    --positions /scratch/bell/mathur20/MQU/ch3/revise/age/sites/priv.${pop}.del.sites \
    --keep /scratch/bell/mathur20/MQU/ch3/revise/lists/${pop}.best66.list \
    --extract-FORMAT-info GT \
    --out priv.${pop}.del
done


vcftools --vcf /scratch/bell/mathur20/MQU/ch3/revise/vcfs/best66.old.auto.noSing.nomiss.rename \
--positions /scratch/bell/mathur20/MQU/ch3/revise/age/sites/shared.del.sites \
--keep /scratch/bell/mathur20/MQU/ch3/revise/lists/az.best66.list \
--keep /scratch/bell/mathur20/MQU/ch3/revise/lists/wtx.best66.list \
--extract-FORMAT-info GT \
--out shared.del