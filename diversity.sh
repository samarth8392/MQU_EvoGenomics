##########################################################################
###                             Samarth Mathur                          ###
###                        PhD Candidate, DeWoody Lab                   ###
###                           Purdue University                         ###
###                                                                     ###
###     Date Created: 08/04/20                 Last Modified: 08/23/20 	###
###########################################################################
###########################################################################
###                     diversity.sh                   					###
###########################################################################

cd $SLURM_SUBMIT_DIR
module load bioinfo
module load vcftools
#module load samtools
#module load BEDTools/2.29.0
#module load bedops/2.4.29


# Get theta-pi from vcf in a sliding window

# Scan windows

cd /scratch/snyder/m/mathur20/MQU/ch3/diversity/run2/

for pop in az wtx etx mx
do
	vcftools --vcf /scratch/snyder/m/mathur20/MQU/ch3/variants/angsd/best66.auto.longch.noSing.rename2.vcf \
	--window-pi 100000 \
	--window-pi-step 50000 \
	--keep pops/$pop \
	--out pi/$pop.best66.100kb
done

