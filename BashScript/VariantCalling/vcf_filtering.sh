##########################################################################
###                          Samarth Mathur, PhD                     	###
###                        The Ohio State University                 	###
###                                                                     ###
###########################################################################
###########################################################################
###                     vcf_filtering.sh                   				###
###########################################################################


cd $SLURM_SUBMIT_DIR
module load bioinfo
module load vcftools

# Filter singletons and private doubletons (N=74,709)

vcftools --vcf best66.recode.vcf \
--singletons \
--out /scratch/bell/mathur20/MQU/ch3/revise/vcfstats/best66

cut -f1,2 /scratch/bell/mathur20/MQU/ch3/revise/vcfstats/best66.singletons | \
> /scratch/bell/mathur20/MQU/ch3/revise/vcfstats/best66.singletons.SNPs

vcftools --vcf best66.old.auto.recode.vcf \
--recode --recode-INFO-all \
--exclude-positions /scratch/bell/mathur20/MQU/ch3/revise/vcfstats/best66.singletons.SNPs \
--out best66.noSing


# Filter VCF to retain only autosomes 
#i.e. no variants from sex chromosomes, mitogenome, and unplaced scaffolds.

cd /scratch/bell/mathur20/MQU/ch3/revise/vcfs/

# Only large autosomes
vcftools --vcf best66.noSing.recode.vcf \
--recode --recode-INFO-all --remove-indels --min-alleles 2 --max-alleles 2 \
--chr NC_006088.5 \
--chr NC_006089.5 \
--chr NC_006090.5 \
--chr NC_006091.5 \
--chr NC_006092.5 \
--chr NC_006093.5 \
--chr NC_006094.5 \
--chr NC_006095.5 \
--chr NC_006096.5 \
--chr NC_006097.5 \
--chr NC_006098.5 \
--chr NC_006099.5 \
--chr NC_006100.5 \
--chr NC_006101.5 \
--chr NC_006102.5 \
--chr NC_006103.5 \
--chr NC_006104.5 \
--chr NC_006105.5 \
--chr NC_006106.5 \
--chr NC_006107.5 \
--chr NC_006108.5 \
--chr NC_006109.5 \
--chr NC_006110.5 \
--chr NC_006111.5 \
--chr NC_006112.4 \
--chr NC_006113.5 \
--chr NC_006114.5 \
--chr NC_006115.5 \
--chr NC_028739.2 \
--chr NC_028740.2 \
--chr NC_006119.4 \
--chr NC_008465.4 \
--out  best66.noSing.auto


# Remove SNPs with >10% missingness 

vcftools --vcf best66.noSing.auto.recode.vcf \
--missing-site \
--out /scratch/bell/mathur20/MQU/ch3/revise/vcfstats/best66.noSing.auto

cat /scratch/bell/mathur20/MQU/ch3/revise/vcfstats/best66.noSing.auto.lmiss | \
awk '{if ($6 < 0.1) print $0;}' | cut -f1,2 \
> /scratch/bell/mathur20/MQU/ch3/revise/vcfstats/best66.noSing.auto.nomiss.SNPs

vcftools --vcf best66.noSing.auto.recode.vcf \
--recode --recode-INFO-all \
--positions /scratch/bell/mathur20/MQU/ch3/revise/vcfstats/best66.noSing.auto.nomiss.SNPs \
--out best66.noSing.auto.nomiss

