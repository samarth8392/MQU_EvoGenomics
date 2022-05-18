##########################################################################
###                          Samarth Mathur, PhD                     	###
###                        The Ohio State University                 	###
###                                                                     ###
###########################################################################
###########################################################################
###                     snp_annotation.sh                           	###
###########################################################################

cd $SLURM_SUBMIT_DIR
module purge
module load bioinfo
module load perl

cd /scratch/bell/mathur20/MQU/ch3/revise/annotation/vep/

# All 66
~/ensembl-vep/vep --everything --species gallus_gallus_merged \
-i /scratch/bell/mathur20/MQU/ch3/revise/vcfs/best66.noSing.auto.nomiss.recode.vcf  \
--cache --offline --canonical --allow_non_variant --symbol --sift b \
--output_file best66.auto.noSing.nomiss.sift 


#Only AZ
~/ensembl-vep/vep --everything --species gallus_gallus_merged \
-i /scratch/bell/mathur20/MQU/ch3/revise/vcfs/onlyAZ.recode.vcf  \
--cache --offline --canonical --allow_non_variant --symbol --sift b \
--output_file onlyAZ.sift 

#Only WTX
~/ensembl-vep/vep --everything --species gallus_gallus_merged \
-i /scratch/bell/mathur20/MQU/ch3/revise/vcfs/onlyWTX.recode.vcf  \
--cache --offline --canonical --allow_non_variant --symbol --sift b \
--output_file onlyWTX.sift

# Filter VEP on consequences

cd /scratch/bell/mathur20/MQU/ch3/revise/annotation/vep/results/

for con in stop_gained stop_lost start_lost missense_variant synonymous_variant 5_prime_UTR_variant 3_prime_UTR_variant	intron_variant upstream_gene_variant downstream_gene_variant intergenic_variant
do
	~/ensembl-vep/filter_vep -i ../best66.auto.noSing.nomiss.sift \
	-o all66.${con}.txt \
	-filter "Consequence is ${con}"
done

for pop in onlyAZ onlyWTX
do
	for con in stop_gained stop_lost start_lost missense_variant synonymous_variant 5_prime_UTR_variant 3_prime_UTR_variant	intron_variant upstream_gene_variant downstream_gene_variant intergenic_variant
	do		
		~/ensembl-vep/filter_vep -i ${pop}.sift \
		-o ${pop}.${con}.txt \
		-filter "Consequence is ${con}"
	done
done

# FIlter SNPs on SIFT scores

for sift in deleterious tolerated
do
	~/ensembl-vep/filter_vep -i ../best66.auto.noSing.nomiss.sift \
	-o all66.${sift}.txt \
	-filter "SIFT is ${sift}"
done

~/ensembl-vep/filter_vep -i ../best66.auto.noSing.nomiss.sift \
-o all66.weakdel.txt \
-filter "SIFT < 0.1"

for pop in onlyAZ onlyWTX
do
	for sift in deleterious tolerated
	do
		~/ensembl-vep/filter_vep -i ${pop}.sift \
		-o ${pop}.${sift}.txt \
		-filter "SIFT is ${sift}"
	done
done

for pop in onlyAZ onlyWTX
do
	~/ensembl-vep/filter_vep -i ${pop}.sift \
	-o ${pop}.weakdel.txt \
	-filter "SIFT < 0.1"
done



for pop in onlyAZ onlyWTX
do
	for con in stop_gained stop_lost start_lost missense_variant synonymous_variant 5_prime_UTR_variant 3_prime_UTR_variant	intron_variant upstream_gene_variant downstream_gene_variant intergenic_variant
	do
		grep -v "#" ${pop}.${con}.txt | cut -f2 | uniq | sed -r 's/:/\t/g' > ../snps/${pop}.${con}.SNPs
	done
done
