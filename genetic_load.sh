##########################################################################
###                             Samarth Mathur                          ###
###                        PhD Candidate, DeWoody Lab                   ###
###                           Purdue University                         ###
###                                                                     ###
###     Date Created: 08/18/20                  Last Modified: 08/18/20 ###
###########################################################################
###########################################################################
###                     genetic_load.sh                        			###
###########################################################################

cd $SLURM_SUBMIT_DIR
module purge
module load bioinfo
module load vcftools
#module load use.own
#module load conda-env/easySFS-py3.7.0
#module load R 

cd /scratch/snyder/m/mathur20/MQU/ch3/load/vcfs 

# get stats For del and weakdel sites

ls *vcf | while read -r line
do
	for i in depth site-depth site-pi het indv-freq-burden missing-site site-mean-depth missing-indv TsTv-summary 
	do
		vcftools --vcf best66.weakdel.recode.vcf \
		--$i \
		--out /scratch/snyder/m/mathur20/MQU/ch3/load/results/best66.weakdel.recode.vcf
	done
done



for pop in az wtx etx mx
do
	vcftools --vcf /scratch/snyder/m/mathur20/MQU/ch3/variants/angsd/best66.auto.noSing.final.recode.vcf \
	--recode --recode-INFO-all \
	--positions /scratch/snyder/m/mathur20/MQU/ch3/synnonsyn/sites/del.sites \
	--keep /scratch/snyder/m/mathur20/MQU/ch3/diversity/withSing/fst/$pop.best68.txt \
	--out $pop.best66.del

	vcftools --vcf /scratch/snyder/m/mathur20/MQU/ch3/variants/angsd/best66.auto.noSing.final.recode.vcf \
	--recode --recode-INFO-all \
	--positions /scratch/snyder/m/mathur20/MQU/ch3/synnonsyn/sites/weakdel.sites \
	--keep /scratch/snyder/m/mathur20/MQU/ch3/diversity/withSing/fst/$pop.best68.txt \
	--out $pop.best66.weakdel
	
	vcftools --vcf /scratch/snyder/m/mathur20/MQU/ch3/variants/angsd/best66.auto.noSing.final.recode.vcf \
	--recode --recode-INFO-all \
	--positions /scratch/snyder/m/mathur20/MQU/ch3/synnonsyn/sites/syn.sites \
	--keep /scratch/snyder/m/mathur20/MQU/ch3/diversity/withSing/fst/$pop.best68.txt \
	--out $pop.best68.syn
done


for pop in az wtx etx mx
do
	for i in freq2 counts2 
	do
		vcftools --vcf $pop.best66.del.recode.vcf \
		--$i \
		--out ../freq/$pop.best66.del

		vcftools --vcf $pop.best66.weakdel.recode.vcf \
		--$i \
		--out ../freq/$pop.best66.weakdel

		vcftools --vcf $pop.best66.syn.recode.vcf \
		--$i \
		--out ../freq/$pop.best66.syn
	done
done


for i in del weakdel syn
do
	for j in diff-site diff-site-discordance diff-discordance-matrix
	do
		vcftools --vcf az.best66.$i.recode.vcf \
		--diff wtx.best66.$i.recode.vcf \
		--$j \
		--out ../diff/az.wtx.best66.$i
	done
done



# AZ Del sites = 156,864


ls *vcf | while read -r line
do
	grep -v "#" $line > $line.bed
done

cd ../bed

ls *bed | while read -r line
do
	cut -f1,2,3,4,5,6,7,8,9 --complement $line > ../geno/$line
done
for i in del weakdel benign
do
	for pop in az wtx etx mx 
	do
		if [ $pop == "az" ]
		then
			ind=28
		fi
		if [ $pop == "wtx" ]
		then
			ind=31
		fi
		if [ $pop == "etx" ]
		then
			ind=3
		fi
		if [ $pop == "mx" ]
		then
			ind=4
		fi
		Rscript /scratch/snyder/m/mathur20/MQU/ch3/jobcodes/RScripts/genetic_load.R \
		-f /scratch/snyder/m/mathur20/MQU/ch3/load/geno/$pop.$i.geno \
		-o /scratch/snyder/m/mathur20/MQU/ch3/load/results/het/$pop.$i \
		-i $ind
	done
done
#



cd /scratch/snyder/m/mathur20/MQU/ch3/load/shared_private/

vcftools --vcf /scratch/snyder/m/mathur20/MQU/ch3/variants/angsd/best66.auto.noSing.final.recode.vcf \
--positions /scratch/snyder/m/mathur20/MQU/ch3/mut_age/sites/del_sites_shared_azwtx.txt \
--het \
--out del.shared

vcftools --vcf /scratch/snyder/m/mathur20/MQU/ch3/variants/angsd/best66.auto.noSing.final.recode.vcf \
--positions /scratch/snyder/m/mathur20/MQU/ch3/mut_age/sites/del_sites_private_az.txt \
--het \
--out del.private.az

vcftools --vcf /scratch/snyder/m/mathur20/MQU/ch3/variants/angsd/best66.auto.noSing.final.recode.vcf \
--positions /scratch/snyder/m/mathur20/MQU/ch3/mut_age/sites/del_sites_private_wtx.txt \
--het \
--out del.private.wtx

vcftools --vcf /scratch/snyder/m/mathur20/MQU/ch3/variants/angsd/best66.auto.noSing.final.recode.vcf \
--positions /scratch/snyder/m/mathur20/MQU/ch3/mut_age/sites/weakdel_sites_shared_azwtx.txt \
--het \
--out weakdel.shared

vcftools --vcf /scratch/snyder/m/mathur20/MQU/ch3/variants/angsd/best66.auto.noSing.final.recode.vcf \
--positions /scratch/snyder/m/mathur20/MQU/ch3/mut_age/sites/weakdel_sites_private_az.txt \
--het \
--out weakdel.private.az

vcftools --vcf /scratch/snyder/m/mathur20/MQU/ch3/variants/angsd/best66.auto.noSing.final.recode.vcf \
--positions /scratch/snyder/m/mathur20/MQU/ch3/mut_age/sites/weakdel_sites_private_wtx.txt \
--het \
--out weakdel.private.wtx

vcftools --vcf /scratch/snyder/m/mathur20/MQU/ch3/variants/angsd/best66.auto.noSing.final.recode.vcf \
--positions /scratch/snyder/m/mathur20/MQU/ch3/mut_age/sites/syn_sites_shared_azwtx.txt \
--het \
--out syn.shared

vcftools --vcf /scratch/snyder/m/mathur20/MQU/ch3/variants/angsd/best66.auto.noSing.final.recode.vcf \
--positions /scratch/snyder/m/mathur20/MQU/ch3/mut_age/sites/syn_sites_private_az.txt \
--het \
--out syn.private.az

vcftools --vcf /scratch/snyder/m/mathur20/MQU/ch3/variants/angsd/best66.auto.noSing.final.recode.vcf \
--positions /scratch/snyder/m/mathur20/MQU/ch3/mut_age/sites/syn_sites_private_wtx.txt \
--het \
--out syn.private.wtx

