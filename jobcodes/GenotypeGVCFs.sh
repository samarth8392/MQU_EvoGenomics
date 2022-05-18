 #!/bin/sh -l
#SBATCH -A fnrquail
#SBATCH -N 1
#SBATCH -n 128
#SBATCH -t 14-00:00:00
#SBATCH --mem=220G
#SBATCH --job-name=GenotypeGVCFs_3
#SBATCH -e GenotypeGVCFs_3
#SBATCH -o GenotypeGVCFs_3

##########################################################################
###                             Samarth Mathur                          ###
###                        PhD Candidate, DeWoody Lab                   ###
###                           Purdue University                         ###
###                                                                     ###
###     Date Created: 06/06/20                  Last Modified: 12/03/21 ###
###########################################################################
###########################################################################
###                     vartiant_filtering.sh                   		###
###########################################################################

cd $SLURM_SUBMIT_DIR
module purge
module load bioinfo
#module load bwa
#module load picard-tools
#module load bedops
module load GATK/4.1.6.0
#module load GATK/3.6.0
#module load samtools

# best 66 #

cd /scratch/bell/mathur20/MQU/ch3/revise/gatk/vcf/

#gatk --java-options "-Xmx220g -XX:+UseParallelGC -XX:ParallelGCThreads=128" CombineGVCFs \
#-R /scratch/bell/mathur20/MQU/ch3/align/chick/ref/chicken_genome.fa \
#-V /scratch/bell/mathur20/MQU/ch3/revise/gatk/gvcf/E6536.raw.variants.g.vcf \
#-V /scratch/bell/mathur20/MQU/ch3/revise/gatk/gvcf/E6609.raw.variants.g.vcf \
#-V /scratch/bell/mathur20/MQU/ch3/revise/gatk/gvcf/E6628.raw.variants.g.vcf \
#-V /scratch/bell/mathur20/MQU/ch3/revise/gatk/gvcf/E6846.raw.variants.g.vcf \
#-V /scratch/bell/mathur20/MQU/ch3/revise/gatk/gvcf/E6877.raw.variants.g.vcf \
#-V /scratch/bell/mathur20/MQU/ch3/revise/gatk/gvcf/E7031.raw.variants.g.vcf \
#-V /scratch/bell/mathur20/MQU/ch3/revise/gatk/gvcf/E7125.raw.variants.g.vcf \
#-V /scratch/bell/mathur20/MQU/ch3/revise/gatk/gvcf/E7146.raw.variants.g.vcf \
#-V /scratch/bell/mathur20/MQU/ch3/revise/gatk/gvcf/E7208.raw.variants.g.vcf \
#-V /scratch/bell/mathur20/MQU/ch3/revise/gatk/gvcf/E7220.raw.variants.g.vcf \
#-V /scratch/bell/mathur20/MQU/ch3/revise/gatk/gvcf/E7563.raw.variants.g.vcf \
#-V /scratch/bell/mathur20/MQU/ch3/revise/gatk/gvcf/E7746.raw.variants.g.vcf \
#-V /scratch/bell/mathur20/MQU/ch3/revise/gatk/gvcf/E7747.raw.variants.g.vcf \
#-V /scratch/bell/mathur20/MQU/ch3/revise/gatk/gvcf/E7751.raw.variants.g.vcf \
#-V /scratch/bell/mathur20/MQU/ch3/revise/gatk/gvcf/E7752.raw.variants.g.vcf \
#-V /scratch/bell/mathur20/MQU/ch3/revise/gatk/gvcf/E7927.raw.variants.g.vcf \
#-V /scratch/bell/mathur20/MQU/ch3/revise/gatk/gvcf/E7932.raw.variants.g.vcf \
#-V /scratch/bell/mathur20/MQU/ch3/revise/gatk/gvcf/E7934.raw.variants.g.vcf \
#-V /scratch/bell/mathur20/MQU/ch3/revise/gatk/gvcf/E7946.raw.variants.g.vcf \
#-V /scratch/bell/mathur20/MQU/ch3/revise/gatk/gvcf/E7969.raw.variants.g.vcf \
#-V /scratch/bell/mathur20/MQU/ch3/revise/gatk/gvcf/E8013.raw.variants.g.vcf \
#-V /scratch/bell/mathur20/MQU/ch3/revise/gatk/gvcf/E8017.raw.variants.g.vcf \
#-V /scratch/bell/mathur20/MQU/ch3/revise/gatk/gvcf/E8024.raw.variants.g.vcf \
#-V /scratch/bell/mathur20/MQU/ch3/revise/gatk/gvcf/E8025.raw.variants.g.vcf \
#-V /scratch/bell/mathur20/MQU/ch3/revise/gatk/gvcf/E8030.raw.variants.g.vcf \
#-V /scratch/bell/mathur20/MQU/ch3/revise/gatk/gvcf/E8031.raw.variants.g.vcf \
#-V /scratch/bell/mathur20/MQU/ch3/revise/gatk/gvcf/E8032.raw.variants.g.vcf \
#-V /scratch/bell/mathur20/MQU/ch3/revise/gatk/gvcf/E8142.raw.variants.g.vcf \
#-V /scratch/bell/mathur20/MQU/ch3/revise/gatk/gvcf/E8946.raw.variants.g.vcf \
#-V /scratch/bell/mathur20/MQU/ch3/revise/gatk/gvcf/E8947.raw.variants.g.vcf \
#-V /scratch/bell/mathur20/MQU/ch3/revise/gatk/gvcf/E8948.raw.variants.g.vcf \
#-V /scratch/bell/mathur20/MQU/ch3/revise/gatk/gvcf/E8949.raw.variants.g.vcf \
#-V /scratch/bell/mathur20/MQU/ch3/revise/gatk/gvcf/E8954.raw.variants.g.vcf \
#-V /scratch/bell/mathur20/MQU/ch3/revise/gatk/gvcf/E9030.raw.variants.g.vcf \
#-V /scratch/bell/mathur20/MQU/ch3/revise/gatk/gvcf/E9031.raw.variants.g.vcf \
#-V /scratch/bell/mathur20/MQU/ch3/revise/gatk/gvcf/E9032.raw.variants.g.vcf \
#-V /scratch/bell/mathur20/MQU/ch3/revise/gatk/gvcf/E9033.raw.variants.g.vcf \
#-V /scratch/bell/mathur20/MQU/ch3/revise/gatk/gvcf/E9034.raw.variants.g.vcf \
#-V /scratch/bell/mathur20/MQU/ch3/revise/gatk/gvcf/E9035.raw.variants.g.vcf \
#-V /scratch/bell/mathur20/MQU/ch3/revise/gatk/gvcf/E9036.raw.variants.g.vcf \
#-V /scratch/bell/mathur20/MQU/ch3/revise/gatk/gvcf/E9037.raw.variants.g.vcf \
#-V /scratch/bell/mathur20/MQU/ch3/revise/gatk/gvcf/E9038.raw.variants.g.vcf \
#-V /scratch/bell/mathur20/MQU/ch3/revise/gatk/gvcf/E9039.raw.variants.g.vcf \
#-V /scratch/bell/mathur20/MQU/ch3/revise/gatk/gvcf/E9040.raw.variants.g.vcf \
#-V /scratch/bell/mathur20/MQU/ch3/revise/gatk/gvcf/E9041.raw.variants.g.vcf \
#-V /scratch/bell/mathur20/MQU/ch3/revise/gatk/gvcf/E9042.raw.variants.g.vcf \
#-V /scratch/bell/mathur20/MQU/ch3/revise/gatk/gvcf/E9043.raw.variants.g.vcf \
#-V /scratch/bell/mathur20/MQU/ch3/revise/gatk/gvcf/E9044.raw.variants.g.vcf \
#-V /scratch/bell/mathur20/MQU/ch3/revise/gatk/gvcf/E9045.raw.variants.g.vcf \
#-V /scratch/bell/mathur20/MQU/ch3/revise/gatk/gvcf/E9046.raw.variants.g.vcf \
#-V /scratch/bell/mathur20/MQU/ch3/revise/gatk/gvcf/E9047.raw.variants.g.vcf \
#-V /scratch/bell/mathur20/MQU/ch3/revise/gatk/gvcf/E9048.raw.variants.g.vcf \
#-V /scratch/bell/mathur20/MQU/ch3/revise/gatk/gvcf/E9049.raw.variants.g.vcf \
#-V /scratch/bell/mathur20/MQU/ch3/revise/gatk/gvcf/E9050.raw.variants.g.vcf \
#-V /scratch/bell/mathur20/MQU/ch3/revise/gatk/gvcf/E9051.raw.variants.g.vcf \
#-V /scratch/bell/mathur20/MQU/ch3/revise/gatk/gvcf/E9052.raw.variants.g.vcf \
#-V /scratch/bell/mathur20/MQU/ch3/revise/gatk/gvcf/E9053.raw.variants.g.vcf \
#-V /scratch/bell/mathur20/MQU/ch3/revise/gatk/gvcf/E9054.raw.variants.g.vcf \
#-V /scratch/bell/mathur20/MQU/ch3/revise/gatk/gvcf/E9055.raw.variants.g.vcf \
#-V /scratch/bell/mathur20/MQU/ch3/revise/gatk/gvcf/E9056.raw.variants.g.vcf \
#-V /scratch/bell/mathur20/MQU/ch3/revise/gatk/gvcf/E9057.raw.variants.g.vcf \
#-V /scratch/bell/mathur20/MQU/ch3/revise/gatk/gvcf/E9058.raw.variants.g.vcf \
#-V /scratch/bell/mathur20/MQU/ch3/revise/gatk/gvcf/E9059.raw.variants.g.vcf \
#-V /scratch/bell/mathur20/MQU/ch3/revise/gatk/gvcf/E9067.raw.variants.g.vcf \
#-V /scratch/bell/mathur20/MQU/ch3/revise/gatk/gvcf/E9569.raw.variants.g.vcf \
#-V /scratch/bell/mathur20/MQU/ch3/revise/gatk/gvcf/E9570.raw.variants.g.vcf \
#-O best66.raw.variants.g.vcf.gz


#PicardCommandLine CreateSequenceDictionary \
#reference=/scratch/bell/mathur20/MQU/ch3/align/chick/ref/chicken_genome.fa \
#output=/scratch/bell/mathur20/MQU/ch3/align/chick/ref/chicken_genome.dict


#gatk --java-options "-Xmx220g -XX:+UseParallelGC -XX:ParallelGCThreads=128" GenotypeGVCFs \
#-R /scratch/bell/mathur20/MQU/ch3/align/chick/ref/chicken_genome.fa \
#-V best66.raw.variants.g.vcf.gz \
#-O best66.raw.variants.vcf.gz

gatk --java-options "-Xmx220g -XX:+UseParallelGC -XX:ParallelGCThreads=64" SelectVariants \
-R /scratch/bell/mathur20/MQU/ch3/align/chick/ref/chicken_genome.fa \
-V best66.raw.variants.vcf.gz \
--select-type-to-include SNP \
-O best66.raw.SNPs.vcf.gz

gatk --java-options "-Xmx220g" VariantFiltration \
-R /scratch/bell/mathur20/MQU/ch3/align/chick/ref/chicken_genome.fa \
-V best66.raw.SNPs.vcf.gz \
-filter "QD < 2.0" --filter-name "QD2" \
-filter "MQ < 20.0" --filter-name "MQ20" \
-filter "MQRankSum < -3.0" --filter-name "MQRankSum-3.0" \
-filter "MQRankSum > 3.0" --filter-name "MQRankSum3.0" \
-filter "ReadPosRankSum < -3.0" --filter-name "ReadPosRankSum-3.0" \
-filter "ReadPosRankSum > 3.0" --filter-name "ReadPosRankSum3.0" \
-filter "SOR > 5.0" --filter-name "SOR5" \
-filter "FS > 40.0" --filter-name "FS40" \
-filter "AN < 105.0" --filter-name "AN105" \
-filter "AF < 0.05" --filter-name "MAF0.05" \
-O best66.filterflag.SNPs.vcf.gz

gatk SelectVariants \
-R /scratch/bell/mathur20/MQU/ch3/align/chick/ref/chicken_genome.fa \
-V best66.filterflag.SNPs.vcf.gz \
-select 'vc.isNotFiltered()' \
-O best66.final.filtered.SNPs.vcf.gz
