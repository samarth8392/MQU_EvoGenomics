#!/bin/sh -l
#SBATCH -A fnrquail
#SBATCH -N 1
#SBATCH -n 10
#SBATCH -t 5-00:00:00
#SBATCH --job-name=load1
#SBATCH -e load1
#SBATCH -o load1

##########################################################################
###                             Samarth Mathur                          ###
###                        PhD Candidate, DeWoody Lab                   ###
###                           Purdue University                         ###
###                                                                     ###
###     Date Created: 06/06/20                  Last Modified: 03/02/22 ###
###########################################################################
###########################################################################
###                     genetic_load.sh                   		        ###
###########################################################################

cd $SLURM_SUBMIT_DIR
module purge
module load bioinfo
module load vcftools
module load R

# Get frequency and genotypes for deleterious (N= 19631),tolerated (N=37238), 
# synonymous (N=210505) and nonsynonymous (N=201311) mutations

cd /scratch/bell/mathur20/MQU/ch3/revise/load/freq/

#for i in deleterious tolerated synonymous nonsynonymous
#do
#    for j in az wtx mx etx
#    do
#        for p in freq counts
#        do
#            vcftools --vcf /scratch/bell/mathur20/MQU/ch3/revise/vcfs/best66.old.auto.poly.nomiss.recode.vcf \
#            --positions /scratch/bell/mathur20/MQU/ch3/revise/annotation/snps/best66.$i.SNPs \
#            --keep /scratch/bell/mathur20/MQU/ch3/revise/lists/$j.best66.list \
#            --$p \
#            --out $j.best66.$i
#        done
#    done
#done
##
#cd /scratch/bell/mathur20/MQU/ch3/revise/load/geno/
##
#for i in deleterious tolerated synonymous nonsynonymous
#do
##    for j in az wtx mx etx
#    do
#        vcftools --vcf /scratch/bell/mathur20/MQU/ch3/revise/vcfs/best66.old.auto.poly.nomiss.recode.vcf \
#        --positions /scratch/bell/mathur20/MQU/ch3/revise/annotation/snps/best66.$i.SNPs \
#        --keep /scratch/bell/mathur20/MQU/ch3/revise/lists/$j.best66.list \
#        --extract-FORMAT-info GT \
#        --out $j.best66.$i
#    done
#done
###

# Remove header and extract MAF

#for i in deleterious tolerated synonymous nonsynonymous
#do
#    for j in az wtx etx mx
#    do
#        less /scratch/bell/mathur20/MQU/ch3/revise/load/freq/$j.best66.$i.frq | \
#        sed '1d'| cut -f6 | cut -f2 -d ":" > /scratch/bell/mathur20/MQU/ch3/revise/load/freq/$j.best66.$i.MAF
#    done
#done

####### Load based on ROHs ########


cd /scratch/bell/mathur20/MQU/ch3/revise/load/rohs/freq/


for len in short medium long
do
    while read -a ind
    do
        vcftools --vcf /scratch/bell/mathur20/MQU/ch3/revise/vcfs/best66.old.auto.poly.nomiss.recode.vcf \
        --bed /scratch/bell/mathur20/MQU/ch3/revise/roh/ROHbed/${ind}.${len}ROHs.bed \
        --indv ${ind} \
        --freq \
        --out ${ind}.${len}
        done < /scratch/bell/mathur20/MQU/ch3/revise/lists/best66.samples.list
done


cd /scratch/bell/mathur20/MQU/ch3/revise/load/rohs/geno/

for len in short medium long
do
    while read -a ind
    do
        vcftools --vcf /scratch/bell/mathur20/MQU/ch3/revise/vcfs/best66.old.auto.poly.nomiss.recode.vcf \
        --bed /scratch/bell/mathur20/MQU/ch3/revise/roh/ROHbed/${ind}.${len}ROHs.bed \
        --indv ${ind} \
        --extract-FORMAT-info GT \
        --out ${ind}.${len}
    done < /scratch/bell/mathur20/MQU/ch3/revise/lists/best66.samples.list
done

# Remove header and extract MAF
cd /scratch/bell/mathur20/MQU/ch3/revise/load/rohs/freq/

for len in short medium long
do
    while read -a ind
    do
        less ${ind}.${len}.frq | \
        sed '1d'| cut -f6 | cut -f2 -d ":" > ${ind}.${len}.MAF
        done < /scratch/bell/mathur20/MQU/ch3/revise/lists/best66.samples.list
done


