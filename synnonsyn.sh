##########################################################################
###                             Samarth Mathur                          ###
###                        PhD Candidate, DeWoody Lab                   ###
###                           Purdue University                         ###
###                                                                     ###
###     Date Created: 08/18/20                  Last Modified: 08/18/20 ###
###########################################################################
###########################################################################
###                     structure.sh                        				###
###########################################################################

# To get sites based on impact (deleterious, weak-del, and benign) for both pop AZ, TX
cd $SLURM_SUBMIT_DIR
module purge
module load bioinfo
module load vcftools
#module load plink

cd /scratch/snyder/m/mathur20/MQU/ch3/synnonsyn/final/

#Get syn sites (N=454,896)

~/ensembl-vep/filter_vep \
-i /scratch/snyder/m/mathur20/MQU/ch3/annotation/vep/wtx.best68.auto.longch.noSing.vep \
-o wtx.best66.poly.syn.txt \
--filter "Consequence is synonymous_variant"


#Get deleterious sites (N=156,864)
# DEL (for each pop Az, TX)
~/ensembl-vep/filter_vep \
-i /scratch/snyder/m/mathur20/MQU/ch3/annotation/vep/wtx.best68.auto.longch.noSing.vep \
-o  wtx.best66.poly.del.txt \
--filter "SIFT is deleterious"

#Get weakly-deleterious sites (for each pop Az, TX)
~/ensembl-vep/filter_vep \
-i /scratch/snyder/m/mathur20/MQU/ch3/annotation/vep/wtx.best68.auto.longch.noSing.vep \
-o  wtx.best66.poly.weakdel.txt \
--filter "SIFT is tolerated"

#Get benign sites (N=594,976) (for each pop Az, TX)
~/ensembl-vep/filter_vep \
-i /scratch/snyder/m/mathur20/MQU/ch3/annotation/vep/best68.auto.longch.noSing.vep \
-o /scratch/snyder/m/mathur20/MQU/ch3/synnonsyn/best66.benign.txt \
--filter "SIFT >= 0.1 or Consequence is synonymous_variant"


grep "synonymous_variant" best66.syn.txt | cut -f2 | awk '{gsub(":", " "); print $0}' \
> sites/syn.sites

grep "IMPACT=" best66.weakdel.txt | cut -f2 | awk '{gsub(":", " "); print $0}' | uniq \
> sites/weakdel.sites


vcftools --vcf /scratch/snyder/m/mathur20/MQU/ch3/variants/angsd/best66.auto.noSing.final.recode.vcf \
--recode --recode-INFO-all \
--positions sites/weakdel.sites \
--out weakbest66.del

cut -f2 best68.auto.longch.sny  \
| awk '{gsub(":", " "); print $0}' > best68.auto.longch.sny.bed


cd /scratch/snyder/m/mathur20/MQU/ch3/synnonsyn/final/share_priv/

vcftools --vcf /scratch/snyder/m/mathur20/MQU/ch3/variants/angsd/wtx.best66.auto.noSing.poly.final.recode.vcf \
--het \
--positions ../sites/del_sites_shared_azwtx.txt \
--out wtx.best66.poly.del.shared

vcftools --vcf /scratch/snyder/m/mathur20/MQU/ch3/variants/angsd/wtx.best66.auto.noSing.poly.final.recode.vcf \
--het \
--positions ../sites/del_sites_private_wtx.txt \
--out wtx.best66.poly.del.priv

vcftools --vcf /scratch/snyder/m/mathur20/MQU/ch3/variants/angsd/wtx.best66.auto.noSing.poly.final.recode.vcf \
--het \
--positions ../sites/weakdel_sites_shared_azwtx.txt \
--out wtx.best66.poly.weakdel.shared

vcftools --vcf /scratch/snyder/m/mathur20/MQU/ch3/variants/angsd/wtx.best66.auto.noSing.poly.final.recode.vcf \
--het \
--positions ../sites/weakdel_sites_private_wtx.txt \
--out wtx.best66.poly.weakdel.priv

vcftools --vcf /scratch/snyder/m/mathur20/MQU/ch3/variants/angsd/wtx.best66.auto.noSing.poly.final.recode.vcf \
--het \
--positions ../sites/syn_sites_shared_azwtx.txt \
--out wtx.best66.poly.syn.shared

vcftools --vcf /scratch/snyder/m/mathur20/MQU/ch3/variants/angsd/wtx.best66.auto.noSing.poly.final.recode.vcf \
--het \
--positions ../sites/syn_sites_private_wtx.txt \
--out wtx.best66.poly.syn.priv


#DIVERSITY
for i in site-pi het
do
	vcftools --vcf best68.auto.longch.sny.vcf \
	--$i \
	--out stats/best68.auto.longch.sny
done

#Get non-syn sites (N=1017729)
#non-syn = missense + inframe_insertion + inframe_deletion + stop_lost + frameshift_variant

# Has sift score

grep "SIFT=" ../annotation/vep/best68.auto.longch.angsd.mostsevere.vep > best68.auto.longch.nonsyn

cut -f2 best68.auto.longch.nonsyn | awk '{gsub(":", " "); print $0}' > best68.auto.longch.nosyn.bed

vcftools --vcf /scratch/snyder/m/mathur20/MQU/ch3/variants/angsd/best68.auto.longchr.angsd.final.recode.vcf \
--recode --recode-INFO-all \
--positions best68.auto.longch.nosyn.bed \
--out best68.auto.longch.nonsyn

# Deleterious, weakly deleterious, benign
#Del
cd sites
grep "deleterious(" best68.auto.longch.nonsyn |cut -f2 \
| awk '{gsub(":", " "); print $0}' > deleterious.bed

vcftools --vcf /scratch/snyder/m/mathur20/MQU/ch3/variants/angsd/best68.auto.longchr.angsd.final.recode.vcf \
--recode --recode-INFO-all --maf 0.008 \
--positions deleterious.bed \
--out ../best68.deleterious

#Weak Del
~/ensembl-vep/filter_vep \
-i /scratch/snyder/m/mathur20/MQU/ch3/annotation/vep/best68.auto.longch.angsd.vep \
-o /scratch/snyder/m/mathur20/MQU/ch3/synnonsyn//sites/best68.weaklydeleterious.txt \
--filter "SIFT >= 0.05 and SIFT < 0.1"

grep "IMPACT=" best68.weaklydeleterious.txt | cut -f2 \
| awk '{gsub(":", " "); print $0}' > weakly_deleterious.bed

vcftools --vcf /scratch/snyder/m/mathur20/MQU/ch3/variants/angsd/best68.auto.longchr.angsd.final.recode.vcf \
--recode --recode-INFO-all --maf 0.008 \
--positions weakly_deleterious.bed \
--out ../best68.weaklydeleterious

#Benign
~/ensembl-vep/filter_vep \
-i /scratch/snyder/m/mathur20/MQU/ch3/annotation/vep/best68.auto.longch.angsd.vep \
-o /scratch/snyder/m/mathur20/MQU/ch3/synnonsyn//sites/best68.benign.txt \
--filter "SIFT >= 0.1 or Consequence is synonymous_variant"
