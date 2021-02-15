##########################################################################
###                             Samarth Mathur                          ###
###                        PhD Candidate, DeWoody Lab                   ###
###                           Purdue University                         ###
###                                                                     ###
###########################################################################
###########################################################################
###                     admixture.sh                        			###
###########################################################################

# Script to run admixture analysis from called genotypes #
# but first we need to change the name of chromosomes to 1-33


cd $SLURM_SUBMIT_DIR
module purge
#module load bioinfo
#module load plink

cd /scratch/snyder/m/mathur20/MQU/ch3/structure/admix/onlysyn/

#change chr to integer codes
cat /scratch/snyder/m/mathur20/MQU/ch3/synnonsyn/best66.syn.recode.vcf \
| sed -e 's/NC_006088.5/1/' \
-e 's/NC_006089.5/2/' \
-e 's/NC_006090.5/3/' \
-e 's/NC_006091.5/4/' \
-e 's/NC_006092.5/5/' \
-e 's/NC_006093.5/6/' \
-e 's/NC_006094.5/7/' \
-e 's/NC_006095.5/8/' \
-e 's/NC_006096.5/9/' \
-e 's/NC_006097.5/10/' \
-e 's/NC_006098.5/11/' \
-e 's/NC_006099.5/12/' \
-e 's/NC_006100.5/13/' \
-e 's/NC_006101.5/14/' \
-e 's/NC_006102.5/15/' \
-e 's/NC_006103.5/16/' \
-e 's/NC_006104.5/17/' \
-e 's/NC_006105.5/18/' \
-e 's/NC_006106.5/19/' \
-e 's/NC_006107.5/20/' \
-e 's/NC_006108.5/21/' \
-e 's/NC_006109.5/22/' \
-e 's/NC_006110.5/23/' \
-e 's/NC_006111.5/24/' \
-e 's/NC_006112.4/25/' \
-e 's/NC_006113.5/26/' \
-e 's/NC_006114.5/27/' \
-e 's/NC_006115.5/28/' \
-e 's/NC_028739.2/30/' \
-e 's/NC_028740.2/31/' \
-e 's/NC_006119.4/32/' \
-e 's/NC_008465.4/33/' \
> best66.syn.rename.vcf

# Select only autosomes

plink --make-bed --chr-set 33 no-xy no-mt \
--vcf best66.syn.rename.vcf \
--out best66.syn

for K in 1 2 3 4 5 6 7 8 9
do
	/scratch/snyder/m/mathur20/MQU/ch3/softwares/admixture_linux-1.3.0/admixture \
 	--cv best66.syn.bed $K -j2 
done
