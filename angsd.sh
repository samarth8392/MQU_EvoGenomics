##########################################################################
###                             Samarth Mathur                          ###
###                        PhD Candidate, DeWoody Lab                   ###
###                           Purdue University                         ###
###                                                                     ###
###     Date Created: 06/29/20                  Last Modified: 09/02/20 ###
###########################################################################
###########################################################################
###                     angsd.sh 			                  			###
###########################################################################

cd $SLURM_SUBMIT_DIR

module purge
module purge
module load bioinfo
module load samtools
module load R/3.4.2
module load bedops
module load vcftools

### 1D SFS ####
#AZ
/scratch/snyder/m/mathur20/MQU/2019/softwares/angsd/angsd -P 50 \
-GL 1 -doPost 3 -nInd 4 -minQ 20 -minMapQ 30 -doMajorMinor 1 -doMaf 1 \
-anc /scratch/snyder/m/mathur20/MQU/ch3/align/chick/ref/chicken_genome.fa \
-bam /scratch/snyder/m/mathur20/MQU/ch3/angsd/bamlist/az_mqu.list \
-pest /scratch/snyder/m/mathur20/MQU/ch3/angsd/sfs/az_chick.sfs \
-out /scratch/snyder/m/mathur20/MQU/ch3/angsd/sfs/az_chick2

/scratch/snyder/m/mathur20/MQU/2019/softwares/angsd/misc/realSFS -P 50 \
/scratch/snyder/m/mathur20/MQU/ch3/angsd/sfs/ -maxIter 1000 \
> /scratch/snyder/m/mathur20/MQU/ch3/angsd/sfs/az_chick2.sfs

#Bootstrap

/scratch/snyder/m/mathur20/MQU/2019/softwares/angsd/misc/realSFS -P 50 \
/scratch/snyder/m/mathur20/MQU/ch3/angsd/sfs/az_chick2.saf.idx \
-bootstrap 100 \
> /scratch/snyder/m/mathur20/MQU/ch3/angsd/sfs/az_chick_boostrap2.sfs

#MX

/scratch/snyder/m/mathur20/MQU/2019/softwares/angsd/angsd -P 50 \
-GL 1 -doPost 3 -nInd 4 -minQ 20 -minMapQ 30 -doMajorMinor 1 -doMaf 1 \
-anc /scratch/snyder/m/mathur20/MQU/ch3/align/chick/ref/chicken_genome.fa \
-bam /scratch/snyder/m/mathur20/MQU/ch3/angsd/bamlist/mx_mqu.list \
-pest /scratch/snyder/m/mathur20/MQU/ch3/angsd/sfs/mx_chick.sfs \
-out /scratch/snyder/m/mathur20/MQU/ch3/angsd/sfs/mx_chick2

/scratch/snyder/m/mathur20/MQU/2019/softwares/angsd/misc/realSFS -P 50 \
/scratch/snyder/m/mathur20/MQU/ch3/angsd/sfs/ -maxIter 1000 \
> /scratch/snyder/m/mathur20/MQU/ch3/angsd/sfs/mx_chick2.sfs

#Bootstrap

/scratch/snyder/m/mathur20/MQU/2019/softwares/angsd/misc/realSFS -P 50 \
/scratch/snyder/m/mathur20/MQU/ch3/angsd/sfs/mx_chick2.saf.idx \
-bootstrap 100 \
> /scratch/snyder/m/mathur20/MQU/ch3/angsd/sfs/mx_chick_boostrap2.sfs

#TX

/scratch/snyder/m/mathur20/MQU/2019/softwares/angsd/angsd -P 150 \
-GL 1 -doPost 3 -nInd 31 -minQ 20 -minMapQ 30 -doMajorMinor 1 -doMaf 1 \
-anc /scratch/snyder/m/mathur20/MQU/ch3/align/chick/ref/chicken_genome.fa \
-bam /scratch/snyder/m/mathur20/MQU/ch3/angsd/bamlist/tx_mqu.list \
-pest /scratch/snyder/m/mathur20/MQU/ch3/angsd/sfs/tx_chick.sfs \
-out /scratch/snyder/m/mathur20/MQU/ch3/angsd/sfs/tx_chick2

/scratch/snyder/m/mathur20/MQU/2019/softwares/angsd/misc/realSFS -P 150 \
/scratch/snyder/m/mathur20/MQU/ch3/angsd/sfs/tx_chick.saf.idx -maxIter 1000 \
> /scratch/snyder/m/mathur20/MQU/ch3/angsd/sfs/tx_chick2.sfs

#Bootstrap

/scratch/snyder/m/mathur20/MQU/2019/softwares/angsd/misc/realSFS -P 150 \
/scratch/snyder/m/mathur20/MQU/ch3/angsd/sfs/tx_chick2.saf.idx \
-bootstrap 100 \
> /scratch/snyder/m/mathur20/MQU/ch3/angsd/sfs/tx_chick_boostrap2.sfs

#########################

### Theta ###

#MX (Same for other populations i.e AZ and TX)

/scratch/snyder/m/mathur20/MQU/2019/softwares/angsd/misc/realSFS saf2theta -P 50 \
/scratch/snyder/m/mathur20/MQU/ch3/angsd/sfs/mx_chick.saf.idx \
-maxiter 1000 -bootstrap 1000 \
-sfs /scratch/snyder/m/mathur20/MQU/ch3/angsd/sfs/mx_chick.sfs \
-outname /scratch/snyder/m/mathur20/MQU/ch3/angsd/theta/mx_chick

/scratch/snyder/m/mathur20/MQU/2019/softwares/angsd/misc/thetaStat do_stat \
/scratch/snyder/m/mathur20/MQU/ch3/angsd/theta/mx_chick.thetas.idx \
-win 100000 -step 50000 \
-outnames /scratch/snyder/m/mathur20/MQU/ch3/angsd/theta/mx_chick_10050

/scratch/snyder/m/mathur20/MQU/2019/softwares/angsd/misc/thetaStat do_stat \
/scratch/snyder/m/mathur20/MQU/ch3/angsd/theta/mx_chick.thetas.idx \
-outnames /scratch/snyder/m/mathur20/MQU/ch3/angsd/theta/mx_chick


# PLINK FORMAT
/scratch/snyder/m/mathur20/MQU/2019/softwares/angsd/angsd -P 150 \
-doPlink 2 -doGeno -4 -doPost 1 -doMajorMinor 1 -GL 1 -doCounts 1 \
-doMaf 2 -postCutoff 0.99  -SNP_pval 1e-6 -geno_minDepth 5 \
-bam /scratch/snyder/m/mathur20/MQU/ch3/angsd/bamlist/mx_mqu.list \
-out /scratch/snyder/m/mathur20/MQU/ch3/angsd/plink/mx_chick

#########################

### Relatedness ###

/scratch/snyder/m/mathur20/MQU/2019/softwares/angsd/angsd -P 100 \
-GL 1 -doGlf 1 -nInd 66 -doHWE 1 -minMapQ 20 \
-doMajorMinor 1 -doMaf 3 -uniqueOnly 1 \
-setMinDepthInd 1  -minInd 5 \
-SNP_pval 1e-6 \
-skipTriallelic 1 -minMaf 0.05 -minHWEpval 0.01 \
-bam /scratch/snyder/m/mathur20/MQU/ch3/angsd/bamlist/all_mqu.list \
-out /scratch/snyder/m/mathur20/MQU/ch3/angsd/relatedness/all_chick

/scratch/snyder/m/mathur20/MQU/2019/softwares/angsd/misc/ibs \
-glf /scratch/snyder/m/mathur20/MQU/ch3/angsd/relatedness/all_chick.glf.gz \
-model 0 \
-nInd 66 -allpairs 1 \
-outFileName /scratch/snyder/m/mathur20/MQU/ch3/angsd/relatedness/all_chick

Rscript \
  -e "source('/scratch/snyder/m/mathur20/MQU/2019/softwares/angsd/R/read_IBS.R')" \
  -e "res = do_derived_stats(read_ibspair_model0('/scratch/snyder/m/mathur20/MQU/ch3/angsd/relatedness/all_chick.ibspair'))" \
  -e "print(res[,c('ind1', 'ind2', 'nSites', 'Kin', 'R0', 'R1') ])" \
  > /scratch/snyder/m/mathur20/MQU/ch3/angsd/relatedness/all_chick_relate.txt

#########################

### Admixture ###
#This will run for K=1 to 10, 10 times (using a nested loop)

/scratch/snyder/m/mathur20/MQU/2019/softwares/angsd/angsd -P 100 \
-nInd 66 -minQ 20 -minMapQ 30 -uniqueOnly 1 \
-GL 1 -doGlf 2 -doHWE 1 -doDepth 1 -doCounts 1 \
-doMajorMinor 1 -doMaf 3 \
-setMinDepthInd 1 -minInd 5 \
-SNP_pval 1e-6 -setMaxDepth 500 \
-skipTriallelic 1 -minMaf 0.05 -minHWEpval 0.01 \
-bam /scratch/snyder/m/mathur20/MQU/ch3/angsd/bamlist/all_mqu.list \
-out /scratch/snyder/m/mathur20/MQU/ch3/angsd/gl/all_chick

j="1"
while [ $j -lt 11 ]
do
	for i in 1 11
	do
		/scratch/snyder/m/mathur20/MQU/2019/softwares/angsd/misc/NGSadmix -P 100 \
		-K ${i} -minMaf 0.05 -maxiter 50000 -tol 1e-9 -tolLike50 1e-9 \
		-likes /scratch/snyder/m/mathur20/MQU/ch3/angsd/gl/all_chick.beagle.gz \
		-outfiles /scratch/snyder/m/mathur20/MQU/ch3/angsd/r${j}/all_chick_k${i} &> /scratch/snyder/m/mathur20/MQU/ch3/angsd/r${j}/all_chick_k${i}.log
 		i=$[$i+1]
	done
	j=$[$j+1]
done

#########################

### Genotype calling ###

/scratch/snyder/m/mathur20/MQU/ch3/softwares/angsd/angsd -P 400 \
-dovcf 1 -GL 2 -doCounts 1 -nInd 66 -doHWE 1 -doPost 1 -doMajorMinor 1 -doMaf 3 \
--ignore-RG 0 -dogeno 1 -geno_minDepth 5 -minMapQ 20 \
-skipTriallelic 1 -minQ 20  -doDepth 1 -SNP_pval 1e-6 \
-setMinDepthInd 5 -doGlf 2 \
-uniqueOnly 1 -dumpCounts 2 -doplink 2 -dosnpstat 1 \
-anc /scratch/snyder/m/mathur20/MQU/ch3/align/chick/ref/chicken_genome.fa \
-bam /scratch/snyder/m/mathur20/MQU/ch3/align/chick/final_bam/final_bam.list \
-out /scratch/snyder/m/mathur20/MQU/ch3/angsd/snps/all66

bcftools convert -O v -o /scratch/snyder/m/mathur20/MQU/ch3/angsd/snps/all66.vcf \
/scratch/snyder/m/mathur20/MQU/ch3/angsd/snps/all66.bcf

