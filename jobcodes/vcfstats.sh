#!/bin/sh -l
#SBATCH -A fnrquail
#SBATCH -N 1
#SBATCH -n 5
#SBATCH -t 14-00:00:00
#SBATCH --job-name=vcfstats
#SBATCH -e vcfstats
#SBATCH -o vcfstats


cd $SLURM_SUBMIT_DIR
module load bioinfo
module load vcftools


cd /scratch/bell/mathur20/MQU/ch3/revise/vcfstats/

#for i in depth het relatedness2 missing-indv missing-site
#do
#	vcftools --vcf /scratch/bell/mathur20/MQU/ch3/revise/vcfs/best66.old.auto.poly.nomiss.recode.vcf \
#	--$i \
#	--out best66.old.auto.poly.nomiss
#done

#vcftools --vcf /scratch/bell/mathur20/MQU/ch3/revise/vcfs/best66.old.auto.poly.nomiss.recode.vcf \
#--weir-fst-pop /scratch/bell/mathur20/MQU/ch3/revise/lists/az.best66.list \
#--weir-fst-pop /scratch/bell/mathur20/MQU/ch3/revise/lists/wtx.best66.list \
#--fst-window-size 100000 --fst-window-step 100000 \
#--out best66.old.auto.poly.nomiss

## Weir and Cockerham mean Fst estimate: 0.029349
## Weir and Cockerham weighted Fst estimate: 0.046356


for i in az wtx
do
	vcftools --vcf /scratch/bell/mathur20/MQU/ch3/revise/vcfs/best66.old.auto.noSing.nomiss.recode.vcf \
	--window-pi 1000 --window-pi-step 1000 \
	--keep /scratch/bell/mathur20/MQU/ch3/revise/lists/$i.best66.list \
	--out $i.best66.old.auto.noSing.nomiss_1kb

	vcftools --vcf /scratch/bell/mathur20/MQU/ch3/revise/vcfs/best66.old.auto.noSing.nomiss.recode.vcf \
	--window-pi 100000 --window-pi-step 100000 \
	--keep /scratch/bell/mathur20/MQU/ch3/revise/lists/$i.best66.list \
	--out $i.best66.old.auto.noSing.nomiss_100kb

done

