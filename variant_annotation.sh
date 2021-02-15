##########################################################################
###                             Samarth Mathur                          ###
###                        PhD Candidate, DeWoody Lab                   ###
###                           Purdue University                         ###
###                                                                     ###
###     Date Created: 09/04/20                  Last Modified: 09/04/20 ###
###########################################################################
###########################################################################
###                     variant_annotation.sh                        	###
###########################################################################

# Get variant annotations from VEP and SNPeff and merge results

cd $SLURM_SUBMIT_DIR
module purge
module load bioinfo
module load perl

cd /scratch/snyder/m/mathur20/MQU/ch3/annotation/vep/

~/ensembl-vep/vep --everything --species gallus_gallus_merged \
-i /scratch/snyder/m/mathur20/MQU/ch3/variants/angsd/best68.auto.longchr.angsd.noSing.recode.vcf \
--cache --offline \
--canonical --allow_non_variant --symbol \
--output_file best68.auto.longch.noSing.vep \
--stats_file best68.auto.longch.noSing.vep.stats 


~/ensembl-vep/vep --everything --species gallus_gallus_merged \
-i /scratch/snyder/m/mathur20/MQU/ch3/variants/angsd/az.best66.auto.noSing.final.recode.vcf \
--cache --offline \
--canonical --allow_non_variant --symbol \
--output_file az.best68.auto.longch.noSing.vep \
--stats_file az.best68.auto.longch.noSing.vep.stats 

~/ensembl-vep/vep --everything --species gallus_gallus_merged \
-i /scratch/snyder/m/mathur20/MQU/ch3/variants/angsd/best66.auto.noSing.final.poly.recode.vcf \
--cache --offline \
--canonical --allow_non_variant --symbol \
--output_file best66.auto.longch.noSing.poly.vep \
--stats_file best66.auto.longch.noSing.poly.vep.stats 

~/ensembl-vep/vep --species gallus_gallus_merged \
-i /scratch/snyder/m/mathur20/MQU/ch3/variants/angsd/wtx.best66.auto.noSing.final.recode.vcf \
--cache --offline --most_severe \
--output_file wtx.best66.auto.longch.poly.mostsevere.vep \
--stats_file wtx.best66.auto.longch.poly.mostsevere.vep.stats

#SNPeff

module load snpEff/4.3
module load bedops

cd /scratch/snyder/m/mathur20/MQU/ch3/annotation/snpeff/
snpEff build -c snpEff.config -gff3 -v chick &> build.logfile.txt

snpEff ann -stats -c snpEff.config \
-no-downstream -no-intergenic -no-intron -no-upstream -no-utr -v \
chick /scratch/snyder/m/mathur20/MQU/ch3/variants/angsd/best68.auto.longchr.angsd.final.recode.vcf \
> best68.auto.longchr.angsd.final.snpeff.vcf

#remove warn and no annotation
grep -v "WARNING_" best68.auto.longchr.angsd.final.snpeff.vcf \
| grep "ANN=" > best68.auto.longchr.angsd.final.snpeff.nowarn.onlyANN.vcf

# Get SIFT scores

grep "SIFT=" best68.auto.longch.angsd.vep > best68.auto.longch.angsd.SIFT

#Merge SNPeff and VEP

cd /scratch/snyder/m/mathur20/MQU/ch3/annotation/snpeff

grep "ANN=" best68.auto.longchr.angsd.final.snpeff.nowarn.onlyANN.vcf | cut -f1,2,8 > snpeff.pos.good.eff

cd /scratch/snyder/m/mathur20/MQU/ch3/annotation/vep/
grep "IMPACT=" best68.auto.longch.angsd.vep | cut -f2,14 | sed 's/':'/'" "'/' > vep.pos.eff


cd /scratch/snyder/m/mathur20/MQU/ch3/annotation/vep/

grep "BIOTYPE=protein_coding" vep/vep.pos.eff > vep/vep.pcode.pos.eff
grep "protein_coding" snpeff/snpeff.pos.eff > snpeff/snpeff.pcode.pos.eff

rm vep/vep.pos.eff
rm snpeff/snpeff.pos.eff

cat vep/vep.pcode.pos.eff | awk '{gsub(";", " "); print $0}' | awk '{ print $1, $2, $3}' > vep/vep.impact

cut -f2 vep/vep.pcode.pos.eff | cut -f1 -d ';' | cut -f2 -d "=" > vep/vep.impact
cut -f3 snpeff/snpeff.pcode.pos.eff | cut -f3 -d '|' > snpeff/snpeff.impact

paste -d " " vep/vep.pcode.pos.eff vep/vep.impact snpeff/snpeff.pcode.pos.eff snpeff/snpeff.impact \
| cut -f1,2,4,5,6,8 > vep.snpeff.merged
