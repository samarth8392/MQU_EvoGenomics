##########################################################################
###                             Samarth Mathur                          ###
###                        PhD Candidate, DeWoody Lab                   ###
###                           Purdue University                         ###
###                                                                     ###
###     Date Created: 09/01/20                  Last Modified: 09/09/20 ###
###########################################################################
###########################################################################
###                     roh_analysis.sh                        			###
###########################################################################

cd $SLURM_SUBMIT_DIR
module purge
module load bioinfo
#module load fastStructure
module load plink/1.90b6.4
#module load vcftools

Remove mito, and sex chr #
vcftools --vcf /scratch/snyder/m/mathur20/MQU/ch3/deep_cov/variants/deep_coverage_HardFilterSNPs.vcf \
--recode --recode-INFO-all \
--not-chr NC_006126.5 \
--not-chr NC_006127.5 \
--not-chr NC_040902.1 \
--out /scratch/snyder/m/mathur20/MQU/ch3/deep_cov/variants/deep_coverage_auto


# VCF to plink format

vcftools --vcf /scratch/snyder/m/mathur20/MQU/ch3/deep_cov/variants/deep_coverage_HardFilterSNPs.vcf \
--plink --out /scratch/snyder/m/mathur20/MQU/ch3/deep_cov/ROH/deep_coverage_ROH

plink --vcf /scratch/snyder/m/mathur20/MQU/ch3/impute/beagle/best68.auto.longchr.angsd.phased.vcf \
--allow-extra-chr -out /scratch/snyder/m/mathur20/MQU/ch3/roh/best68.auto.longchr.angsd.phased

# Test runs

for a in 10 20 50 100 1000 10000
do
	for b in 0.01 0.05 0.1 
	do
		for c in 5 50 100 500 1000
		do
			plink --bim best68.auto.longchr.angsd.phased.bim \
			--bed best68.auto.longchr.angsd.phased.bed \
			--fam best68.auto.longchr.angsd.phased.fam \
			--allow-extra-chr --no-sex --no-parents --no-pheno \
			--homozyg-window-het 3 \
			--homozyg-window-missing $c \
			--homozyg-window-snp $a \
			--homozyg-density $a \
			--homozyg-gap $a \
			--homozyg-window-threshold $b \
			--homozyg-snp $a \
			--homozyg-kb $a \
			--out tryruns/phase/best68.auto.longchr.angsd.phased.$a.$b.$c.roh
		done
	done
done

#cd /scratch/snyder/m/mathur20/MQU/ch3/roh/plink1.9/
cd /scratch/snyder/m/mathur20/MQU/ch3/roh/plink1.9/run2
plink --chr-set 33 no-xy no-mt  --make-bed \
-vcf /scratch/snyder/m/mathur20/MQU/ch3/variants/angsd/best66.auto.longch.noSing.rename2.vcf \
--out best66.auto.longchr.noSing

##Genome sequence length
cd /scratch/snyder/m/mathur20/MQU/ch3/align/chick/ref
cat chicken_genome.fa  | awk '$0 ~ ">" {if (NR > 1) {print c;} c=0;printf substr($0,2,100) "\t"; } $0 !~ ">" {c+=length($0);} END { print c; }' \
> chicken_chr_sequence.length 
cd /scratch/snyder/m/mathur20/MQU/ch3/roh/plink1.9/

#Final run

plink --bim best68.auto.longchr.angsd.phased.bim \
--bed best68.auto.longchr.angsd.phased.bed \
--fam best68.auto.longchr.angsd.phased.fam \
--homozyg --allow-extra-chr \
--homozyg-snp 50 \
--homozyg-kb 100 \
--homozyg-density 10 \
--homozyg-gap 100 \
--homozyg-window-snp 50 \
--homozyg-window-het 5 \
--out final
