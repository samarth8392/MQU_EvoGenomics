#!/bin/sh -l
#SBATCH -A fnrdewoody
#SBATCH -N 1
#SBATCH -n 20
#SBATCH -t 10-00:00:00
#SBATCH --job-name=renameChr2
#SBATCH -e renameChr2
#SBATCH -o renameChr2


cd $SLURM_SUBMIT_DIR
module load bioinfo
module load vcftools

# Filter VCF to retain only autosomes 
#i.e. no variants from sex chromosomes, mitogenome, and unplaced scaffolds.

cd /scratch/bell/mathur20/MQU/ch3/revise/vcfs/

#change chr to integer codes
cat best66.old.auto.poly.nomiss.recode.vcf \
| sed -e 's/NC_006088.5/chr1/' \
-e 's/NC_006089.5/chr2/' \
-e 's/NC_006090.5/chr3/' \
-e 's/NC_006091.5/chr4/' \
-e 's/NC_006092.5/chr5/' \
-e 's/NC_006093.5/chr6/' \
-e 's/NC_006094.5/chr7/' \
-e 's/NC_006095.5/chr8/' \
-e 's/NC_006096.5/chr9/' \
-e 's/NC_006097.5/chr10/' \
-e 's/NC_006098.5/chr11/' \
-e 's/NC_006099.5/chr12/' \
-e 's/NC_006100.5/chr13/' \
-e 's/NC_006101.5/chr14/' \
-e 's/NC_006102.5/chr15/' \
-e 's/NC_006103.5/chr16/' \
-e 's/NC_006104.5/chr17/' \
-e 's/NC_006105.5/chr18/' \
-e 's/NC_006106.5/chr19/' \
-e 's/NC_006107.5/chr20/' \
-e 's/NC_006108.5/chr21/' \
-e 's/NC_006109.5/chr22/' \
-e 's/NC_006110.5/chr23/' \
-e 's/NC_006111.5/chr24/' \
-e 's/NC_006112.4/chr25/' \
-e 's/NC_006113.5/chr26/' \
-e 's/NC_006114.5/chr27/' \
-e 's/NC_006115.5/chr28/' \
-e 's/NC_028739.2/chr30/' \
-e 's/NC_028740.2/chr31/' \
-e 's/NC_006119.4/chr32/' \
-e 's/NC_008465.4/chr33/' \
> best66.old.auto.poly.nomiss.renameChr.vcf
